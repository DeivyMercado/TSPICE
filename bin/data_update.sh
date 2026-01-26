#!/bin/bash

# Script to update/check data files
# Currently a placeholder as tspice includes kernels directly in package data

echo "Checking data files..."
DATA_DIR="src/tspice/data"

if [ -d "$DATA_DIR" ]; then
    echo "Data directory exists at $DATA_DIR"
    ls -l "$DATA_DIR"
else
    echo "Warning: Data directory not found at $DATA_DIR"
fi

echo "Data check complete."
