#!/usr/bin/env python3
"""
Test runner for the micro-pKa screening pipeline.

This script runs the complete test suite and provides detailed reporting.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def run_tests(test_dir="tests", verbose=False, coverage=False):
    """Run the test suite."""
    print("🧪 Running Micro-pKa Pipeline Tests")
    print("=" * 50)
    
    # Check if pytest is available
    try:
        import pytest
    except ImportError:
        print("❌ pytest not found. Installing...")
        subprocess.run([sys.executable, "-m", "pip", "install", "pytest", "pytest-cov"])
    
    # Build pytest command
    cmd = [sys.executable, "-m", "pytest", test_dir]
    
    if verbose:
        cmd.append("-v")
    
    if coverage:
        cmd.extend(["--cov=src", "--cov-report=html", "--cov-report=term"])
    
    # Run tests
    print(f"Running command: {' '.join(cmd)}")
    print()
    
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode == 0:
        print("\n✅ All tests passed!")
        if coverage:
            print("📊 Coverage report generated in htmlcov/index.html")
    else:
        print("\n❌ Some tests failed!")
        sys.exit(1)


def run_examples():
    """Run example scripts to verify functionality."""
    print("\n🚀 Running Examples")
    print("=" * 30)
    
    examples_dir = Path("examples")
    if not examples_dir.exists():
        print("❌ Examples directory not found")
        return
    
    example_files = list(examples_dir.glob("*.py"))
    if not example_files:
        print("❌ No example files found")
        return
    
    for example_file in example_files:
        print(f"Running {example_file.name}...")
        try:
            result = subprocess.run([sys.executable, str(example_file)], 
                                  capture_output=True, text=True, timeout=300)
            if result.returncode == 0:
                print(f"✅ {example_file.name} completed successfully")
            else:
                print(f"❌ {example_file.name} failed:")
                print(result.stderr)
        except subprocess.TimeoutExpired:
            print(f"⏰ {example_file.name} timed out")
        except Exception as e:
            print(f"❌ {example_file.name} error: {e}")


def check_imports():
    """Check if all required modules can be imported."""
    print("\n🔍 Checking Imports")
    print("=" * 20)
    
    required_modules = [
        "pandas",
        "numpy", 
        "rdkit",
        "torch",
        "torch_geometric"
    ]
    
    missing_modules = []
    
    for module in required_modules:
        try:
            __import__(module)
            print(f"✅ {module}")
        except ImportError:
            print(f"❌ {module}")
            missing_modules.append(module)
    
    if missing_modules:
        print(f"\n❌ Missing modules: {', '.join(missing_modules)}")
        print("Please install missing dependencies:")
        print("pip install -r requirements.txt")
        return False
    
    print("\n✅ All required modules available")
    return True


def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(description="Test runner for micro-pKa pipeline")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--coverage", "-c", action="store_true", help="Generate coverage report")
    parser.add_argument("--examples", "-e", action="store_true", help="Run examples")
    parser.add_argument("--imports", "-i", action="store_true", help="Check imports only")
    parser.add_argument("--all", "-a", action="store_true", help="Run all checks")
    
    args = parser.parse_args()
    
    if args.all:
        args.imports = True
        args.examples = True
        args.coverage = True
    
    # Check imports first
    if args.imports or args.all:
        if not check_imports():
            sys.exit(1)
    
    # Run tests
    if not args.imports or args.all:
        run_tests(verbose=args.verbose, coverage=args.coverage)
    
    # Run examples
    if args.examples or args.all:
        run_examples()
    
    print("\n🎉 Test suite completed!")


if __name__ == "__main__":
    main()
