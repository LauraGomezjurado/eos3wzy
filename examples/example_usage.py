#!/usr/bin/env python3
"""
Example usage of the micro-pKa screening pipeline.

This script demonstrates how to use the pipeline for drug discovery
screening with various molecular libraries.
"""

import os
import sys
import pandas as pd
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from micro_pka_pipeline import MicroPKAPipeline


def example_single_molecule():
    """Example: Process a single molecule."""
    print("=== Single Molecule Example ===")
    
    # Initialize pipeline
    pipeline = MicroPKAPipeline(
        cache_dir="./cache",
        n_processes=2,
        use_caching=True
    )
    
    # Process a single molecule
    smiles = "CCO"  # Ethanol
    result = pipeline.process_single_molecule(smiles)
    
    print(f"SMILES: {result['smiles']}")
    print(f"Valid: {result['valid']}")
    
    if result['valid']:
        stats = result['statistics']
        print(f"Min pKa: {stats['min_pka']:.2f}")
        print(f"Avg pKa: {stats['avg_pka']:.2f}")
        print(f"Max pKa: {stats['max_pka']:.2f}")
        print(f"Number of sites: {stats['num_sites']}")
        print(f"Processing time: {result['metadata']['processing_time']:.2f}s")
    else:
        print(f"Error: {result['error']}")


def example_batch_processing():
    """Example: Process a batch of molecules."""
    print("\n=== Batch Processing Example ===")
    
    # Initialize pipeline
    pipeline = MicroPKAPipeline(
        cache_dir="./cache",
        n_processes=4,
        use_caching=True
    )
    
    # Define test molecules
    test_molecules = [
        "CCO",           # Ethanol
        "CCN",           # Ethylamine
        "CCOO",          # Ethyl peroxide
        "CC(C)O",        # Isopropanol
        "CC(C)N",        # Isopropylamine
        "invalid_smiles" # Invalid input
    ]
    
    print(f"Processing {len(test_molecules)} molecules...")
    
    # Process batch
    results = pipeline.process_batch(test_molecules)
    
    # Analyze results
    valid_count = sum(1 for r in results if r['valid'])
    print(f"Successfully processed: {valid_count}/{len(results)} molecules")
    
    # Print results for valid molecules
    for result in results:
        if result['valid']:
            stats = result['statistics']
            print(f"{result['smiles']}: pKa = {stats['avg_pka']:.2f} "
                  f"(sites: {stats['num_sites']})")
        else:
            print(f"{result['smiles']}: ERROR - {result['error']}")


def example_file_processing():
    """Example: Process molecules from a file."""
    print("\n=== File Processing Example ===")
    
    # Create sample input file
    input_file = "sample_input.csv"
    output_file = "sample_output.csv"
    
    # Create sample data
    sample_data = pd.DataFrame({
        'smiles': [
            'CCO',           # Ethanol
            'CCN',           # Ethylamine
            'CCOO',          # Ethyl peroxide
            'CC(C)O',        # Isopropanol
            'CC(C)N',        # Isopropylamine
            'CC(C)(C)O',     # tert-Butanol
            'CC(C)(C)N',     # tert-Butylamine
            'CC(C)(C)CO',    # tert-Butyl alcohol
            'CC(C)(C)CN',    # tert-Butylamine
            'CC(C)(C)CCO'    # tert-Butyl alcohol
        ]
    })
    sample_data.to_csv(input_file, index=False)
    print(f"Created sample input file: {input_file}")
    
    # Initialize pipeline
    pipeline = MicroPKAPipeline(
        cache_dir="./cache",
        n_processes=4,
        use_caching=True
    )
    
    # Run screening
    print("Running screening...")
    pipeline.run_screening(input_file, output_file)
    
    # Load and display results
    results_df = pd.read_csv(output_file)
    print(f"\nResults saved to: {output_file}")
    print(f"Processed {len(results_df)} molecules")
    
    # Show summary statistics
    valid_results = results_df.dropna(subset=['avg_pka'])
    if len(valid_results) > 0:
        print(f"\nSummary Statistics:")
        print(f"Average pKa range: {valid_results['avg_pka'].min():.2f} - {valid_results['avg_pka'].max():.2f}")
        print(f"Mean pKa: {valid_results['avg_pka'].mean():.2f}")
        print(f"Standard deviation: {valid_results['avg_pka'].std():.2f}")
    
    # Clean up
    os.remove(input_file)
    print(f"\nCleaned up temporary files")


def example_drug_like_molecules():
    """Example: Process drug-like molecules."""
    print("\n=== Drug-like Molecules Example ===")
    
    # Drug-like molecules for testing
    drug_molecules = [
        "CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1",  # Chlorpromazine-like
        "COc1ccccc1N(C)C",                           # Diphenhydramine-like
        "CC(=O)Nc1ccc(O)cc1",                        # Acetaminophen-like
        "CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1",    # Another antipsychotic
        "COc1ccccc1N(C)C",                           # Another antihistamine
    ]
    
    # Initialize pipeline
    pipeline = MicroPKAPipeline(
        cache_dir="./cache",
        n_processes=4,
        use_caching=True
    )
    
    print(f"Processing {len(drug_molecules)} drug-like molecules...")
    
    # Process molecules
    results = pipeline.process_batch(drug_molecules)
    
    # Analyze results
    valid_results = [r for r in results if r['valid']]
    print(f"Successfully processed: {len(valid_results)}/{len(results)} molecules")
    
    # Group by pKa ranges
    pka_ranges = {
        'Very acidic (pKa < 4)': 0,
        'Acidic (4 ≤ pKa < 7)': 0,
        'Neutral (7 ≤ pKa < 9)': 0,
        'Basic (9 ≤ pKa < 12)': 0,
        'Very basic (pKa ≥ 12)': 0
    }
    
    for result in valid_results:
        avg_pka = result['statistics']['avg_pka']
        if avg_pka < 4:
            pka_ranges['Very acidic (pKa < 4)'] += 1
        elif avg_pka < 7:
            pka_ranges['Acidic (4 ≤ pKa < 7)'] += 1
        elif avg_pka < 9:
            pka_ranges['Neutral (7 ≤ pKa < 9)'] += 1
        elif avg_pka < 12:
            pka_ranges['Basic (9 ≤ pKa < 12)'] += 1
        else:
            pka_ranges['Very basic (pKa ≥ 12)'] += 1
    
    print("\nDistribution by pKa range:")
    for range_name, count in pka_ranges.items():
        print(f"  {range_name}: {count} molecules")


def example_performance_benchmark():
    """Example: Performance benchmarking."""
    print("\n=== Performance Benchmark Example ===")
    
    import time
    
    # Test molecules of varying complexity
    test_molecules = [
        "CCO",                    # Simple alcohol
        "CC(C)O",                 # Branched alcohol
        "CC(C)(C)O",             # Tertiary alcohol
        "CC(C)(C)(C)O",          # Quaternary alcohol
        "CC(C)(C)(C)(C)O",       # More complex
        "CC(C)(C)(C)(C)(C)O",    # Even more complex
    ]
    
    # Test different process counts
    process_counts = [1, 2, 4]
    
    for n_processes in process_counts:
        print(f"\nTesting with {n_processes} processes:")
        
        # Initialize pipeline
        pipeline = MicroPKAPipeline(
            cache_dir="./cache",
            n_processes=n_processes,
            use_caching=True
        )
        
        # Time the processing
        start_time = time.time()
        results = pipeline.process_batch(test_molecules)
        end_time = time.time()
        
        total_time = end_time - start_time
        molecules_per_second = len(test_molecules) / total_time
        
        print(f"  Total time: {total_time:.2f}s")
        print(f"  Molecules per second: {molecules_per_second:.2f}")
        print(f"  Average time per molecule: {total_time/len(test_molecules):.2f}s")


def main():
    """Run all examples."""
    print("Micro-pKa Screening Pipeline Examples")
    print("=" * 50)
    
    try:
        example_single_molecule()
        example_batch_processing()
        example_file_processing()
        example_drug_like_molecules()
        example_performance_benchmark()
        
        print("\n" + "=" * 50)
        print("All examples completed successfully!")
        
    except Exception as e:
        print(f"Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
