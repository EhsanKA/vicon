# VICON - Viral Sequence Analysis Toolkit

VICON is a Python package for processing and analyzing viral sequence data, with specialized tools for viral genome coverage analysis and sequence alignment.

## Features

- Viral sequence alignment and coverage analysis
- K-mer analysis and sliding window coverage calculations
<!-- - Support for segmented viral genomes (rotavirus, influenza, etc.) -->
- Visualization tools for coverage plots
- Wrapper scripts for vsearch and viralmsa
<!-- - Support for multiple input formats (FASTA, WIG) -->

## Installation

### Local Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/EhsanKA/vicon.git
   cd vicon
   ```

2. Create and activate a conda environment:
   ```bash
   conda env create -f environment.yaml
   conda activate vicon
   ```

3. Install the package:
   ```bash
   pip install -e .
   ```

### Docker Installation

1. Build the Docker image:
   ```bash
   docker build -t vicon .
   ```

2. Run the container:
   ```bash
   docker run -it vicon
   ```

## Usage

### Main Scripts
#### todo 
1. Process viral sequences:
   ```bash
   python -m vicon processing/coverage_analysis.py --input input.fasta
   ```

2. Run k-mer analysis:
   ```bash
   python -m vicon processing/kmer_analysis.py --input input.fasta --kmer 21
   ```

## Wrapper Scripts

### run_vsearch
Wrapper for vsearch sequence clustering and analysis.

Usage:
```bash
run_vsearch -i input.fasta -o output.fasta [options]
```

Options:
- `-i/--input`: Input FASTA file (required)
- `-o/--output`: Output FASTA file (required)
- `-t/--threads`: Number of threads (default: 4)
- `-id/--identity`: Minimum identity threshold (default: 0.97)
- `-m/--memory`: Memory limit in MB (default: 4000)

### run_viralmsa
Wrapper for ViralMSA multiple sequence alignment.

Usage:
```bash
run_viralmsa -i input.fasta -r reference.fasta -o output_dir [options]
```

Options:
- `-i/--input`: Input FASTA file (required)
- `-r/--reference`: Reference genome FASTA (required)
- `-o/--output`: Output directory (required)
- `-t/--threads`: Number of threads (default: 4)
- `-a/--algorithm`: Alignment algorithm (default: "mafft")

## License
This project is licensed under the terms of the MIT license.
