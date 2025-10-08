# API Documentation

## MicroPKAPipeline

Main pipeline class for micro-pKa screening.

### `__init__(cache_dir, n_processes, use_caching, xtb_path)`

Initialize the micro-pKa screening pipeline.

**Parameters:**
- `cache_dir` (str): Directory for caching intermediate results
- `n_processes` (int, optional): Number of processes for multiprocessing. If None, uses CPU count.
- `use_caching` (bool): Whether to use caching for intermediate results
- `xtb_path` (str, optional): Path to xtb executable. If None, tries to find in PATH.

### `process_single_molecule(smiles)`

Process a single molecule through the entire pipeline.

**Parameters:**
- `smiles` (str): SMILES string of the molecule

**Returns:**
- `Dict[str, Any]`: Dictionary containing pKa predictions and metadata

**Example:**
```python
pipeline = MicroPKAPipeline()
result = pipeline.process_single_molecule("CCO")
print(result['statistics']['avg_pka'])
```

### `process_batch(smiles_list)`

Process a batch of molecules using multiprocessing.

**Parameters:**
- `smiles_list` (List[str]): List of SMILES strings to process

**Returns:**
- `List[Dict[str, Any]]`: List of results for each molecule

**Example:**
```python
smiles_list = ["CCO", "CCN", "CCOO"]
results = pipeline.process_batch(smiles_list)
```

### `run_screening(input_file, output_file)`

Run a complete screening from input file to output file.

**Parameters:**
- `input_file` (str): Path to input CSV file with SMILES
- `output_file` (str): Path to output CSV file for results

**Example:**
```python
pipeline.run_screening("input.csv", "output.csv")
```

## DataPreparation

Handles molecular data preparation and validation.

### `prepare_molecule(smiles)`

Prepare and validate a molecule from SMILES.

**Parameters:**
- `smiles` (str): SMILES string

**Returns:**
- `Dict[str, Any]`: Preparation results including RDKit molecule object

## TautomerSearch

Handles tautomer search using GFN2-xTB quantum chemistry.

### `__init__(xtb_path)`

Initialize tautomer search with xtb path.

**Parameters:**
- `xtb_path` (str): Path to xtb executable

### `find_tautomers(mol, smiles)`

Find tautomers using GFN2-xTB calculations.

**Parameters:**
- `mol` (Chem.Mol): RDKit molecule object
- `smiles` (str): Original SMILES string

**Returns:**
- `Dict[str, Any]`: Tautomer search results

## ReactionSiteEnumerator

Enumerates potential reaction sites for pKa prediction.

### `enumerate_sites(tautomers)`

Enumerate reaction sites for pKa prediction.

**Parameters:**
- `tautomers` (List[str]): List of tautomer SMILES strings

**Returns:**
- `Dict[str, Any]`: Reaction site enumeration results

## GraphBasedPredictor

Graph-based pKa predictor using QupKake model.

### `predict_pka(reaction_sites)`

Predict pKa values for reaction sites.

**Parameters:**
- `reaction_sites` (List[Dict[str, Any]]): List of reaction sites

**Returns:**
- `Dict[str, Any]`: pKa predictions

## Result Format

### Single Molecule Result

```python
{
    'smiles': 'CCO',
    'valid': True,
    'pka_predictions': [8.5, 9.2],
    'statistics': {
        'min_pka': 8.5,
        'max_pka': 9.2,
        'avg_pka': 8.85,
        'num_sites': 2
    },
    'metadata': {
        'num_tautomers': 1,
        'num_reaction_sites': 2,
        'processing_time': 2.3
    }
}
```

### Batch Result

List of single molecule results, one for each input SMILES.

### Error Result

```python
{
    'smiles': 'invalid_smiles',
    'valid': False,
    'error': 'Invalid SMILES',
    'processing_time': 0.1
}
```

## Configuration

### Environment Variables

- `XTBPATH`: Path to xtb executable
- `CUDA_VISIBLE_DEVICES`: GPU device selection

### Pipeline Parameters

- `cache_dir`: Directory for caching (default: "./cache")
- `n_processes`: Number of processes (default: CPU count)
- `use_caching`: Enable caching (default: True)
- `xtb_path`: Path to xtb (default: auto-detect)

## Error Handling

The pipeline handles various error conditions gracefully:

- **Invalid SMILES**: Returns error result with validation message
- **Missing xtb**: Raises RuntimeError during initialization
- **File I/O errors**: Logs warnings and continues processing
- **Memory issues**: Reduces batch size automatically

## Performance Tips

1. **Use caching**: Enable caching for repeated calculations
2. **Optimize processes**: Set `n_processes` to CPU count
3. **Batch processing**: Use `process_batch()` for multiple molecules
4. **Memory management**: Process large datasets in chunks

## Examples

### Basic Usage

```python
from src.micro_pka_pipeline import MicroPKAPipeline

# Initialize pipeline
pipeline = MicroPKAPipeline()

# Process single molecule
result = pipeline.process_single_molecule("CCO")
if result['valid']:
    print(f"Average pKa: {result['statistics']['avg_pka']}")
```

### Batch Processing

```python
# Process multiple molecules
smiles_list = ["CCO", "CCN", "CCOO"]
results = pipeline.process_batch(smiles_list)

# Filter valid results
valid_results = [r for r in results if r['valid']]
print(f"Processed {len(valid_results)} molecules successfully")
```

### Custom Configuration

```python
# Custom pipeline configuration
pipeline = MicroPKAPipeline(
    cache_dir="/tmp/cache",
    n_processes=8,
    use_caching=True,
    xtb_path="/opt/xtb/bin/xtb"
)
```

### Error Handling

```python
result = pipeline.process_single_molecule("invalid_smiles")
if not result['valid']:
    print(f"Error: {result['error']}")
```
