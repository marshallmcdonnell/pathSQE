#!/usr/bin/env python
"""
Entry point wrapper for pathSQE driver.
This can be called as: python -m pathSQE.pathSQE_driver
"""

import sys
import os

# Ensure the modules are importable
sys.path.insert(0, os.path.dirname(__file__))


def main():
    """Main entry point for pathSQE driver."""
    # Import and run pathSQE - this will be customized based on your needs
    from . import core
    print(f"pathSQE version {core.__dict__.get('__version__', '0.1.0')}")
    print("To use the driver, run the pathSQE_driver.py or configure your script here.")


if __name__ == "__main__":
    main()
