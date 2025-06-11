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

To run the VICON pipeline, use the following command:

```bash
vicon-run --config path/to/your/config.yaml
```

### Example Configuration

Here's an example of what your configuration file (`config.yaml`) should look like:

```yaml
project_path: "project_path"
virus_name: "orov"
input_sample: "data/orov/samples/samples.fasta"
input_reference: "data/orov/reference/reference.fasta"
email: "email@address.com"
kmer_size: 150
threshold: 147 # shows a tolerance of 150-147 =3 degenerations
l_gene_start: 8000
l_gene_end: 16000
coverage_ratio: 0.5
min_year: 2020
threshold_ratio: 0.01
drop_old_samples: false
drop_mischar_samples: true
```

## License
This project is licensed under the terms of the MIT license.
