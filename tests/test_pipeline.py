"""
Unit tests for the micro-pKa screening pipeline.

This module contains comprehensive tests for all components of the pipeline
to ensure correctness and reproducibility.
"""

import pytest
import tempfile
import os
import pandas as pd
from pathlib import Path
import sys

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from micro_pka_pipeline import (
    MicroPKAPipeline, 
    DataPreparation, 
    TautomerSearch, 
    ReactionSiteEnumerator, 
    GraphBasedPredictor
)


class TestDataPreparation:
    """Test data preparation functionality."""
    
    def test_prepare_valid_molecule(self):
        """Test preparation of a valid molecule."""
        prep = DataPreparation()
        result = prep.prepare_molecule("CCO")  # Ethanol
        
        assert result['valid'] is True
        assert result['mol'] is not None
        assert result['num_atoms'] > 0
        assert result['molecular_weight'] > 0
    
    def test_prepare_invalid_smiles(self):
        """Test handling of invalid SMILES."""
        prep = DataPreparation()
        result = prep.prepare_molecule("invalid_smiles")
        
        assert result['valid'] is False
        assert 'error' in result
    
    def test_prepare_empty_smiles(self):
        """Test handling of empty SMILES."""
        prep = DataPreparation()
        result = prep.prepare_molecule("")
        
        assert result['valid'] is False
        assert 'error' in result


class TestTautomerSearch:
    """Test tautomer search functionality."""
    
    def test_find_tautomers(self):
        """Test tautomer search for a simple molecule."""
        search = TautomerSearch("/usr/bin/xtb")  # Mock path
        from rdkit import Chem
        
        mol = Chem.MolFromSmiles("CCO")
        result = search.find_tautomers(mol, "CCO")
        
        assert 'tautomers' in result
        assert 'num_tautomers' in result
        assert result['num_tautomers'] >= 1
        assert "CCO" in result['tautomers']


class TestReactionSiteEnumerator:
    """Test reaction site enumeration."""
    
    def test_enumerate_sites(self):
        """Test reaction site enumeration."""
        enumerator = ReactionSiteEnumerator()
        tautomers = ["CCO", "CCN"]
        
        result = enumerator.enumerate_sites(tautomers)
        
        assert 'reaction_sites' in result
        assert 'num_sites' in result
        assert result['num_sites'] >= 0


class TestGraphBasedPredictor:
    """Test graph-based predictor."""
    
    def test_predict_pka(self):
        """Test pKa prediction."""
        predictor = GraphBasedPredictor()
        reaction_sites = [
            {'atom_idx': 0, 'atomic_num': 8, 'tautomer': 'CCO'},
            {'atom_idx': 1, 'atomic_num': 1, 'tautomer': 'CCO'}
        ]
        
        result = predictor.predict_pka(reaction_sites)
        
        assert 'predictions' in result
        assert 'num_predictions' in result
        assert result['num_predictions'] == len(reaction_sites)
        assert all(isinstance(p, (int, float)) for p in result['predictions'])


class TestMicroPKAPipeline:
    """Test the main pipeline."""
    
    def setup_method(self):
        """Setup test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.cache_dir = os.path.join(self.temp_dir, "cache")
        
        # Create mock xtb path
        self.xtb_path = "/usr/bin/xtb"
    
    def teardown_method(self):
        """Cleanup test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_pipeline_initialization(self):
        """Test pipeline initialization."""
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=False,
            xtb_path=self.xtb_path
        )
        
        assert pipeline.cache_dir == Path(self.cache_dir)
        assert pipeline.n_processes == 1
        assert pipeline.use_caching is False
    
    def test_process_single_molecule(self):
        """Test processing a single molecule."""
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=False,
            xtb_path=self.xtb_path
        )
        
        result = pipeline.process_single_molecule("CCO")
        
        assert 'smiles' in result
        assert 'valid' in result
        assert result['smiles'] == "CCO"
    
    def test_process_invalid_molecule(self):
        """Test processing an invalid molecule."""
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=False,
            xtb_path=self.xtb_path
        )
        
        result = pipeline.process_single_molecule("invalid_smiles")
        
        assert result['valid'] is False
        assert 'error' in result
    
    def test_process_batch(self):
        """Test batch processing."""
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=False,
            xtb_path=self.xtb_path
        )
        
        smiles_list = ["CCO", "CCN", "invalid_smiles"]
        results = pipeline.process_batch(smiles_list)
        
        assert len(results) == len(smiles_list)
        assert results[0]['smiles'] == "CCO"
        assert results[1]['smiles'] == "CCN"
        assert results[2]['valid'] is False
    
    def test_run_screening(self):
        """Test complete screening workflow."""
        # Create test input file
        input_file = os.path.join(self.temp_dir, "input.csv")
        output_file = os.path.join(self.temp_dir, "output.csv")
        
        test_data = pd.DataFrame({
            'smiles': ['CCO', 'CCN', 'CCOO']
        })
        test_data.to_csv(input_file, index=False)
        
        # Initialize pipeline
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=False,
            xtb_path=self.xtb_path
        )
        
        # Run screening
        pipeline.run_screening(input_file, output_file)
        
        # Check output file exists and has correct format
        assert os.path.exists(output_file)
        
        output_data = pd.read_csv(output_file)
        assert 'smiles' in output_data.columns
        assert 'min_pka' in output_data.columns
        assert 'avg_pka' in output_data.columns
        assert 'max_pka' in output_data.columns
        assert 'num_pkas' in output_data.columns


class TestCaching:
    """Test caching functionality."""
    
    def setup_method(self):
        """Setup test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.cache_dir = os.path.join(self.temp_dir, "cache")
        self.xtb_path = "/usr/bin/xtb"
    
    def teardown_method(self):
        """Cleanup test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_caching_enabled(self):
        """Test caching when enabled."""
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=True,
            xtb_path=self.xtb_path
        )
        
        # Process molecule first time
        result1 = pipeline.process_single_molecule("CCO")
        
        # Process same molecule second time (should use cache)
        result2 = pipeline.process_single_molecule("CCO")
        
        # Results should be the same
        assert result1['smiles'] == result2['smiles']
        assert result1['valid'] == result2['valid']
    
    def test_caching_disabled(self):
        """Test caching when disabled."""
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=False,
            xtb_path=self.xtb_path
        )
        
        # Process molecule
        result = pipeline.process_single_molecule("CCO")
        
        # Cache directory should be empty
        cache_files = list(Path(self.cache_dir).glob("*.json"))
        assert len(cache_files) == 0


class TestErrorHandling:
    """Test error handling and edge cases."""
    
    def setup_method(self):
        """Setup test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.cache_dir = os.path.join(self.temp_dir, "cache")
        self.xtb_path = "/usr/bin/xtb"
    
    def teardown_method(self):
        """Cleanup test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_missing_xtb(self):
        """Test behavior when xtb is not found."""
        with pytest.raises(RuntimeError, match="xtb executable not found"):
            MicroPKAPipeline(
                cache_dir=self.cache_dir,
                n_processes=1,
                use_caching=False,
                xtb_path="/nonexistent/path/to/xtb"
            )
    
    def test_invalid_input_file(self):
        """Test handling of invalid input file."""
        pipeline = MicroPKAPipeline(
            cache_dir=self.cache_dir,
            n_processes=1,
            use_caching=False,
            xtb_path=self.xtb_path
        )
        
        # Create invalid input file
        input_file = os.path.join(self.temp_dir, "invalid.csv")
        output_file = os.path.join(self.temp_dir, "output.csv")
        
        with open(input_file, 'w') as f:
            f.write("invalid,data\n")
        
        with pytest.raises(Exception):
            pipeline.run_screening(input_file, output_file)


if __name__ == "__main__":
    pytest.main([__file__])
