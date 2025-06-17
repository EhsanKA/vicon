# VICON Test Procedure

This document outlines the procedure for testing a VICON installation.

## Prerequisites

1. Ensure you have conda installed
2. Have access to the VICON repository
3. Have the OROV-S test data available

## Installation Test

1. Activate the conda environment:
```bash
conda activate vicon
```

2. Make the executables accessible:
```bash
chmod +x "$CONDA_PREFIX/bin/vicon-run"
chmod +x "$CONDA_PREFIX/bin/viralmsa"
```

3. Run the installation test script:
```bash
python tests/test_installation.py
```

The test script will verify:
- Python version
- All VICON modules are importable
- Required commands are available (`vicon-run`, `viralmsa`, `vsearch`)

Expected output should show:
```
✓ VICON installation appears to be working correctly!
```

## Pipeline Test

1. Run the pipeline with OROV-S configuration:
```bash
vicon-run --config configs/config_orov_s.yaml
```

2. Verify the output files in `results/orov_s_test/`:
- `kmer1.csv` and `kmer2.csv`: Kmer information
- `kmers1_mutations.png` and `kmers2_mutations.png`: Mutation visualizations
- `Histogram_of_the_kmer_coverage_orov_s_test_150.png`: Coverage histogram
- `logs/` directory: Detailed run logs

## Expected Results

A successful test run should:
1. Complete without errors
2. Generate all expected output files
3. Show appropriate progress information during execution
4. Create visualizations in the results directory

## Building and Uploading to Anaconda Cloud

### 1. Prepare the Environment
```bash
# Activate the environment
conda activate vicon

# Install and configure Anaconda client
conda install anaconda-client
conda config --set anaconda_upload yes

# Set environment variables
export PYTHONIOENCODING=utf-8
export CONDA_BLD_SKIP_XATTRS=True
```

### 2. Clean Previous Builds
```bash
# Remove previous build artifacts
rm -rf ~/tmp-conda-build/*
```

### 3. Build the Package
```bash
# Create build directory
mkdir -p ~/tmp-conda-build

# Build the package
CONDA_BLD_SKIP_XATTRS=True conda build conda/ -c bioconda -c conda-forge --croot ~/tmp-conda-build

# Create index
conda index ~/tmp-conda-build
```

### 4. Upload to Anaconda Cloud
```bash
# Upload the package
anaconda upload ~/tmp-conda-build/noarch/vicon-0.2.6-*.tar.bz2
```

### 5. Install from Anaconda Cloud
```bash
# Create a new environment (recommended)
conda create -n vicon python=3.11
conda activate vicon

# Install VICON
conda install -c conda-forge -c bioconda -c eka97 vicon
```

### Important Notes
1. Make sure the version in `conda/meta.yaml` matches the version in `pyproject.toml`
2. Increment the build number in `meta.yaml` for each new upload
3. Ensure all dependencies are correctly specified in `meta.yaml`
4. Test the package in a clean environment after uploading

## Usage Examples

### Basic Usage
```bash
# Run VICON with a configuration file
vicon-run --config configs/config_orov_s.yaml
```

### Configuration File Structure
```yaml
# Example configuration file (config_orov_s.yaml)
input:
  fasta: data/orov/orov_s/samples/OROV-S.fasta
  metadata: data/orov/orov_s/samples/metadata.csv

output:
  directory: results/orov_s_test
  prefix: orov_s_test

parameters:
  dereplication:
    identity: 0.99
    threads: 1
  alignment:
    method: viralmsa
    threads: 1
  processing:
    min_coverage: 0.8
    min_identity: 0.8
```

### Output Files
The pipeline generates several output files in the specified output directory:

1. Kmer Information:
   - `kmer1.csv`: Information about the first kmer
   - `kmer2.csv`: Information about the second kmer

2. Visualizations:
   - `kmers1_mutations.png`: Mutations in the first kmer
   - `kmers2_mutations.png`: Mutations in the second kmer
   - `Histogram_of_the_kmer_coverage_*.png`: Coverage distribution

3. Logs:
   - Detailed logs in the `logs/` directory

### Common Parameters

1. Dereplication:
   - `identity`: Sequence identity threshold (default: 0.99)
   - `threads`: Number of threads to use

2. Alignment:
   - `method`: Alignment method (default: viralmsa)
   - `threads`: Number of threads to use

3. Processing:
   - `min_coverage`: Minimum coverage threshold (default: 0.8)
   - `min_identity`: Minimum identity threshold (default: 0.8)

### Example Workflow

1. Prepare your input files:
   - FASTA file containing viral sequences
   - Metadata CSV file with sample information

2. Create a configuration file:
   - Copy the example configuration
   - Update paths and parameters as needed

3. Run the pipeline:
```bash
vicon-run --config your_config.yaml
```

4. Check the results:
```bash
# List output files
ls results/your_output_directory/

# View kmer information
head -n 5 results/your_output_directory/kmer1.csv

# Check logs
cat results/your_output_directory/logs/pipeline.log
```

### Using ViralMSA Directly
You can also use the ViralMSA tool directly:
```bash
# Get help
viralmsa --help

# Run ViralMSA
viralmsa -s your_sequences.fasta -r reference.fasta -o output_directory
```

## Troubleshooting

If you encounter issues:

1. Permission errors:
   - Ensure executables have proper permissions
   - Run the chmod commands listed above

2. Module import errors:
   - Verify the conda environment is activated
   - Check that all dependencies are installed

3. Command not found errors:
   - Verify the conda environment is activated
   - Check that the executables are in the correct location
   - Ensure proper permissions are set

4. Pipeline errors:
   - Check the input data exists
   - Verify the configuration file is correct
   - Check the logs directory for detailed error messages

5. Build errors:
   - Check version numbers match in all configuration files
   - Verify all dependencies are available in specified channels
   - Check for syntax errors in meta.yaml
   - Ensure build directory has proper permissions

6. Upload errors:
   - Verify Anaconda Cloud credentials
   - Check if package version already exists
   - Ensure build artifacts are present
   - Check network connection 