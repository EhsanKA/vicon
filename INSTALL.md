# VICON Installation Guide

## Prerequisites

### 1. Install Miniconda

#### For Intel Macs:
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

#### For Apple Silicon (M1/M2) Macs:
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
```

#### Installation Steps:
```bash
# Make the installer executable
chmod +x Miniconda3-latest-MacOSX-*.sh

# Run the installer
./Miniconda3-latest-MacOSX-*.sh

# Follow the prompts:
# 1. Press Enter to review the license
# 2. Type 'yes' to accept
# 3. Press Enter to confirm the install location (default is ~/miniconda3)
# 4. Type 'yes' to initialize Miniconda

# Restart your shell
source ~/.bash_profile  # or ~/.zshrc if you use Zsh

# Verify installation
conda --version
```

## Installing VICON

### 1. Clone the Repository
```bash
git clone https://github.com/EhsanKA/vicon.git
cd vicon
```

### 2. Create and Activate Environment
```bash
# Create the environment
conda env create -f environment.yaml

# Activate the environment
conda activate vicon
```

### 3. Build and Install VICON
```bash
# Set environment variables
export PYTHONIOENCODING=utf-8
export CONDA_BLD_SKIP_XATTRS=True

# Create build directory
mkdir -p ~/tmp-conda-build

# Build the package
CONDA_BLD_SKIP_XATTRS=True conda build conda/ -c bioconda -c conda-forge --croot ~/tmp-conda-build

# Create index
conda index ~/tmp-conda-build

# Install the package
conda install -c file://$HOME/tmp-conda-build vicon

# Make Python executable
chmod +x "$CONDA_PREFIX/bin/python"
```

## Installing from Anaconda Cloud

If you want to install directly from Anaconda Cloud:

```bash
# Create a new environment (recommended)
conda create -n vicon python=3.11
conda activate vicon

# Install VICON
conda install -c conda-forge -c bioconda -c eka97 vicon
```

## Running VICON

1. Review the example configuration file:
```bash
cat configs/config_orov_s.yaml
```

2. Run the pipeline:
```bash
python run_pipeline.py --config configs/config_orov_s.yaml
```

## Using ViralMSA

The ViralMSA tool is included in the VICON package and can be used directly from the command line:

```bash
# Check if viralmsa is available
which viralmsa

# Run viralmsa
viralmsa --help
```

For more information about using ViralMSA, please refer to the [ViralMSA documentation](https://github.com/niemasd/ViralMSA).

## Troubleshooting

If you encounter any issues:
1. Ensure all environment variables are set correctly
2. Verify that the conda environment is activated
3. Check that all dependencies are installed correctly
4. Make sure you have the correct permissions for the build directory

## Support

For issues or questions, please:
1. Check the [GitHub repository](https://github.com/EhsanKA/vicon)
2. Open an issue if you encounter any problems
3. Review the documentation in the repository 