# Micro-pKa Screening Pipeline for Drug Discovery

A comprehensive micro-pKa screening pipeline based on the QupKake model (Abarbanel & Hutchison, 2024). This pipeline combines quantum chemistry calculations with graph neural networks to predict micro-pKa values with high accuracy for drug discovery applications.

## 🚀 Features

- **Data Preparation**: Molecular preprocessing and validation using RDKit
- **Tautomer Search**: Uses GFN2-xTB quantum chemistry for accurate tautomer enumeration
- **Reaction-Site Enumeration**: Identifies potential pKa sites in molecules
- **Graph-Based Prediction**: Leverages QupKake model for precise pKa prediction
- **Multiprocessing**: Parallel processing for large-scale screening
- **Intelligent Caching**: Reduces computation time for repeated calculations
- **Containerized Workflow**: Docker support for reproducible results
- **CLI Interface**: Easy-to-use command-line interface for batch screening

## 📊 Performance

Based on the QupKake paper validation:
- **Literature Set**: RMSE ≈ 0.54 pKa units
- **Novartis Set**: RMSE ≈ 0.79 pKa units
- **Processing Speed**: ~2-5 seconds per molecule (depending on complexity)

## 🛠️ Quick Start

### Using Docker (Recommended)

```bash
# Build the Docker image
docker build -t micro-pka-pipeline .

# Run screening
docker run -v $(pwd)/data:/data micro-pka-pipeline \
  python /repo/src/cli.py screen /data/input.csv /data/output.csv
```

### Manual Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Install xtb (GFN2-xTB)
# Download from: https://github.com/grimme-lab/xtb/releases
# Or use conda: conda install -c conda-forge xtb

# Install QupKake
pip install git+https://github.com/Shualdon/QupKake.git
```

### Basic Usage

```bash
# Create sample input
python src/cli.py create-sample sample_input.csv

# Run screening
python src/cli.py screen sample_input.csv results.csv

# Show pipeline info
python src/cli.py info
```

## 📝 Input/Output Format

### Input (CSV)
```csv
input
CCO
CCN
CCOO
```

### Output (CSV)
```csv
smiles,min_pka,avg_pka,max_pka,num_pkas,num_tautomers,processing_time
CCO,8.5,8.5,8.5,1,1,2.3
CCN,9.2,9.2,9.2,1,1,2.1
CCOO,7.8,7.8,7.8,1,1,2.4
```

## 🐍 Python API

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
print(f"Average pKa: {result['statistics']['avg_pka']}")

# Process batch
smiles_list = ["CCO", "CCN", "CCOO"]
results = pipeline.process_batch(smiles_list)

# Run complete screening
pipeline.run_screening("input.csv", "output.csv")
```

## 🧪 Examples

See the [examples/](examples/) directory for comprehensive usage examples:

- Single molecule processing
- Batch processing
- File-based screening
- Drug-like molecule analysis
- Performance benchmarking

## 📚 Documentation

- [API Documentation](docs/API.md) - Complete API reference
- [User Guide](docs/README.md) - Detailed usage instructions
- [Examples](examples/) - Code examples and tutorials

## 🧬 Scientific Background

This pipeline implements the QupKake model described in:

> Abarbanel, S., & Hutchison, G. R. (2024). QupKake: A Graph Neural Network Approach to Predicting Micro-pKa Values. *Journal of Chemical Theory and Computation*. https://doi.org/10.1021/acs.jctc.4c00328

The model combines:
- **Quantum Chemistry**: GFN2-xTB calculations for tautomer search
- **Graph Neural Networks**: Deep learning for pKa prediction
- **Reaction-Site Analysis**: Systematic enumeration of potential pKa sites

## 🏗️ Architecture

The pipeline consists of four main components:

1. **DataPreparation**: Validates and preprocesses molecular inputs
2. **TautomerSearch**: Uses GFN2-xTB to find tautomeric forms
3. **ReactionSiteEnumerator**: Identifies potential pKa sites
4. **GraphBasedPredictor**: Predicts pKa values using QupKake model

## ⚡ Performance Optimization

- **Multiprocessing**: Utilizes all available CPU cores
- **Caching**: Avoids redundant calculations
- **Memory Efficient**: Processes molecules in batches
- **Reproducible**: Deterministic results with same inputs

## 🐳 Docker Support

The pipeline is fully containerized for reproducible workflows:

```bash
# Build image
docker build -t micro-pka-pipeline .

# Run with volume mounting
docker run -v $(pwd)/data:/data micro-pka-pipeline \
  python /repo/src/cli.py screen /data/input.csv /data/output.csv
```

## 🧪 Testing

Run the test suite:

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=src --cov-report=html
```

## 📄 License

This project is licensed under the CC-BY-4.0 License.

## 🤝 Contributing

Contributions are welcome! Please see the contributing guidelines for details.

## 📞 Support

For questions and support, please open an issue on the GitHub repository.

## 📖 Citation

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