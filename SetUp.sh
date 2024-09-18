#!/bin/bash

CONFIG_FILE="custom_config"

echo "Setting up FitCheck..."

read -p "Enter MATLAB path [/usr/local/MATLAB/R2021a/bin/matlab]: " MATLAB_PATH
MATLAB_PATH=${MATLAB_PATH:-/usr/local/MATLAB/R2021a/bin/matlab}

read -p "Enter LCModel path [/usr/local/lcmodel/bin/lcmodel]: " LCMODEL_PATH
LCMODEL_PATH=${LCMODEL_PATH:-/usr/local/lcmodel/bin/lcmodel}

read -p "Enter MRS basis set path [/path/to/basis_set]: " BASIS_SET_PATH
BASIS_SET_PATH=${BASIS_SET_PATH:-/path/to/basis_set}

# Write to custom_config
cat <<EOL > $CONFIG_FILE
# custom_config

MATLAB_PATH="$MATLAB_PATH"
LCMODEL_PATH="$LCMODEL_PATH"
BASIS_SET_PATH="$BASIS_SET_PATH"
EOL

echo "Configuration saved to $CONFIG_FILE."

