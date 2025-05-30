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

3. Dependencies:
   - ViralMSA:
      ```bash
      mkdir -p scripts && cd scripts
      wget "https://github.com/niemasd/ViralMSA/releases/latest/download/ViralMSA.py"
      chmod a+x ViralMSA.py
      cd ../
      ```
   - minimap2:
      ```bash
      cd scripts/
      curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
      cd ../
      ```
      Consider running the following command before the execution each time or add it to your `.bashrc`:
      ```bash
      export PATH="$PWD/scripts/minimap2-2.28_x64-linux:$PATH"
      ```
   - Vsearch:
      ```bash
      conda install bioconda::vsearch
      ```

3. Install the package:
   ```bash
   pip install -e .
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
