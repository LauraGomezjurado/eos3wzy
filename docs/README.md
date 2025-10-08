# Micro-pKa Screening Pipeline

A comprehensive micro-pKa screening pipeline for drug discovery based on the QupKake model (Abarbanel & Hutchison, 2024). This pipeline combines quantum chemistry calculations with graph neural networks to predict micro-pKa values with high accuracy.

## Features

- **Data Preparation**: Molecular preprocessing and validation
- **Tautomer Search**: Uses GFN2-xTB quantum chemistry for tautomer enumeration
- **Reaction-Site Enumeration**: Identifies potential pKa sites
- **Graph-Based Prediction**: Leverages QupKake model for accurate pKa prediction
- **Multiprocessing**: Parallel processing for large-scale screening
- **Caching**: Intelligent caching to reduce computation time
- **Containerization**: Docker support for reproducible workflows
- **CLI Interface**: Easy-to-use command-line interface

## Installation

### Prerequisites

- Python 3.9+
- RDKit
- PyTorch
- xtb (GFN2-xTB quantum chemistry engine)

### Using Docker (Recommended)

```bash
# Build the Docker image
docker build -t micro-pka-pipeline .

# Run the pipeline
docker run -v $(pwd)/data:/data micro-pka-pipeline \
  python /repo/src/cli.py screen /data/input.csv /data/output.csv
```

### Manual Installation

```bash
# Clone the repository
git clone <repository-url>
cd micro-pka-pipeline

# Install dependencies
pip install -r requirements.txt

# Install xtb (GFN2-xTB)
# Download from: https://github.com/grimme-lab/xtb/releases
# Or use conda: conda install -c conda-forge xtb

# Install QupKake
pip install git+https://github.com/Shualdon/QupKake.git
```

## Usage

### Command Line Interface

```bash
# Show pipeline information
python src/cli.py info

# Create a sample input file
python src/cli.py create-sample sample_input.csv

# Run screening
python src/cli.py screen input.csv output.csv

# Run with custom settings
python src/cli.py screen input.csv output.csv \
  --n-processes 8 \
  --cache-dir ./cache \
  --verbose
```

### Python API

```python
from src.micro_pka_pipeline import MicroPKAPipeline

# Initialize pipeline
pipeline = MicroPKAPipeline(
    cache_dir="./cache",
    n_processes=4,
    use_caching=True
)

# Process single molecule
result = pipeline.process_single_molecule("CCO")

# Process batch
smiles_list = ["CCO", "CCN", "CCOO"]
results = pipeline.process_batch(smiles_list)

# Run complete screening
pipeline.run_screening("input.csv", "output.csv")
```

## Input Format

The pipeline expects a CSV file with molecular SMILES strings:

```csv
input
CCO
CCN
CCOO
```

## Output Format

The pipeline generates a CSV file with pKa predictions and statistics:

```csv
smiles,min_pka,avg_pka,max_pka,num_pkas,num_tautomers,processing_time
CCO,8.5,8.5,8.5,1,1,2.3
CCN,9.2,9.2,9.2,1,1,2.1
CCOO,7.8,7.8,7.8,1,1,2.4
```

## Performance

The pipeline is optimized for large-scale screening:

- **Multiprocessing**: Utilizes all available CPU cores
- **Caching**: Avoids redundant calculations
- **Memory Efficient**: Processes molecules in batches
- **Reproducible**: Deterministic results with same inputs

### Benchmark Results

Based on the QupKake paper validation:

- **Literature Set**: RMSE ≈ 0.54 pKa units
- **Novartis Set**: RMSE ≈ 0.79 pKa units
- **Processing Speed**: ~2-5 seconds per molecule (depending on complexity)

## Architecture

The pipeline consists of four main components:

1. **DataPreparation**: Validates and preprocesses molecular inputs
2. **TautomerSearch**: Uses GFN2-xTB to find tautomeric forms
3. **ReactionSiteEnumerator**: Identifies potential pKa sites
4. **GraphBasedPredictor**: Predicts pKa values using QupKake model

## Configuration

### Environment Variables

- `XTBPATH`: Path to xtb executable
- `CUDA_VISIBLE_DEVICES`: GPU device selection (if using GPU)

### Pipeline Parameters

- `cache_dir`: Directory for caching intermediate results
- `n_processes`: Number of parallel processes
- `use_caching`: Enable/disable caching
- `xtb_path`: Path to xtb executable

## Troubleshooting

### Common Issues

1. **xtb not found**: Ensure xtb is installed and in PATH
2. **Memory issues**: Reduce `n_processes` or use smaller batches
3. **CUDA errors**: Set `CUDA_VISIBLE_DEVICES=""` to use CPU only

### Debug Mode

```bash
python src/cli.py screen input.csv output.csv --verbose
```

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@article{abarbanel2024qupkake,
  title={QupKake: A Graph Neural Network Approach to Predicting Micro-pKa Values},
  author={Abarbanel, Samuel and Hutchison, Geoffrey R.},
  journal={Journal of Chemical Theory and Computation},
  year={2024},
  doi={10.1021/acs.jctc.4c00328}
}
```

## License

This project is licensed under the CC-BY-4.0 License.

## Contributing

Contributions are welcome! Please see the contributing guidelines for details.

## Support

For questions and support, please open an issue on the GitHub repository.
