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
   - Depending on your os version, download the miniconda from:
   ```
   https://www.anaconda.com/docs/getting-started/miniconda/install#macos-linux-installation
   ```
   - ViralMSA:
      ```bash
      mkdir -p scripts && cd scripts
      wget "https://github.com/niemasd/ViralMSA/releases/latest/download/ViralMSA.py"
      chmod a+x ViralMSA.py
      cd ../
      ```

3. Install the package:
   ```bash
   pip install -e .
   ```

## Usage

```bash
   python run_pipeline.py --config configs/config_rsva.yaml
   ```

## License
This project is licensed under the terms of the MIT license.
