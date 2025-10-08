#!/usr/bin/env python3
"""
Command-line interface for the micro-pKa screening pipeline.

This CLI provides easy access to the pipeline functionality for batch screening
of molecular libraries for pKa prediction.
"""

import argparse
import sys
import os
from pathlib import Path
import logging

# Add src to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from micro_pka_pipeline import MicroPKAPipeline


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('micro_pka_screening.log')
        ]
    )


def validate_input_file(input_file: str) -> bool:
    """Validate input file exists and has correct format."""
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        return False
    
    # Check if it's a CSV file
    if not input_file.endswith('.csv'):
        print(f"Error: Input file must be a CSV file.")
        return False
    
    return True


def create_sample_input(output_file: str):
    """Create a sample input file for testing."""
    sample_data = [
        "input",
        "Cc1[nH]c2ccccc2c1CCNCc1ccc(CCC(=O)N=O)cc1",
        "CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1",
        "COc1ccccc1N(C)C",
        "CC(=O)Nc1ccc(O)cc1",
        "CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1"
    ]
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(sample_data))
    
    print(f"Sample input file created: {output_file}")


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description='Micro-pKa Screening Pipeline for Drug Discovery',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run screening with default settings
  python cli.py screen input.csv output.csv
  
  # Run with custom settings
  python cli.py screen input.csv output.csv --n-processes 8 --cache-dir ./my_cache
  
  # Create sample input file
  python cli.py create-sample sample_input.csv
  
  # Validate input file
  python cli.py validate input.csv
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Screen command
    screen_parser = subparsers.add_parser('screen', help='Run micro-pKa screening')
    screen_parser.add_argument('input_file', help='Input CSV file with SMILES')
    screen_parser.add_argument('output_file', help='Output CSV file for results')
    screen_parser.add_argument('--cache-dir', default='./cache', 
                              help='Cache directory (default: ./cache)')
    screen_parser.add_argument('--n-processes', type=int, 
                              help='Number of processes (default: CPU count)')
    screen_parser.add_argument('--no-cache', action='store_true', 
                              help='Disable caching')
    screen_parser.add_argument('--xtb-path', help='Path to xtb executable')
    screen_parser.add_argument('--verbose', '-v', action='store_true', 
                              help='Verbose output')
    
    # Create sample command
    sample_parser = subparsers.add_parser('create-sample', 
                                         help='Create sample input file')
    sample_parser.add_argument('output_file', help='Output file name')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', 
                                          help='Validate input file')
    validate_parser.add_argument('input_file', help='Input file to validate')
    
    # Info command
    info_parser = subparsers.add_parser('info', help='Show pipeline information')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    if args.command == 'screen':
        # Setup logging
        setup_logging(args.verbose)
        
        # Validate input
        if not validate_input_file(args.input_file):
            sys.exit(1)
        
        # Create output directory if needed
        output_path = Path(args.output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            # Initialize pipeline
            pipeline = MicroPKAPipeline(
                cache_dir=args.cache_dir,
                n_processes=args.n_processes,
                use_caching=not args.no_cache,
                xtb_path=args.xtb_path
            )
            
            # Run screening
            print(f"Starting micro-pKa screening...")
            print(f"Input: {args.input_file}")
            print(f"Output: {args.output_file}")
            print(f"Processes: {pipeline.n_processes}")
            print(f"Cache: {args.cache_dir}")
            print()
            
            pipeline.run_screening(args.input_file, args.output_file)
            
            print(f"\nScreening completed successfully!")
            print(f"Results saved to: {args.output_file}")
            
        except Exception as e:
            print(f"Error during screening: {e}")
            sys.exit(1)
    
    elif args.command == 'create-sample':
        create_sample_input(args.output_file)
    
    elif args.command == 'validate':
        if validate_input_file(args.input_file):
            print(f"Input file '{args.input_file}' is valid.")
        else:
            sys.exit(1)
    
    elif args.command == 'info':
        print("Micro-pKa Screening Pipeline")
        print("=" * 40)
        print("Version: 1.0.0")
        print("Based on: QupKake model (Abarbanel & Hutchison, 2024)")
        print("Paper: https://doi.org/10.1021/acs.jctc.4c00328")
        print()
        print("Features:")
        print("- Data preparation and molecular preprocessing")
        print("- Tautomer search using GFN2-xTB quantum chemistry")
        print("- Reaction-site enumeration")
        print("- Graph-based pKa prediction")
        print("- Multiprocessing and caching for performance")
        print("- Containerized workflow for reproducibility")
        print()
        print("Requirements:")
        print("- Python 3.9+")
        print("- RDKit")
        print("- PyTorch")
        print("- xtb (GFN2-xTB)")
        print("- QupKake model")


if __name__ == "__main__":
    main()
