#!/bin/bash
set -e

# Compile the modified einverted and install it into $PREFIX/bin using the conda-set compiler
$CC -O2 -o $PREFIX/bin/einverted einverted.c -lemboss -lm

# Install the Python package
$PYTHON -m pip install . --no-deps -vv