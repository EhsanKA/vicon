# VICON - Viral Sequence Analysis Toolkit

VICON is a Python package for processing and analyzing viral sequence data, with specialized tools for viral genome coverage analysis and sequence alignment.

## Features

- Viral sequence alignment and coverage analysis
- K-mer analysis and sliding window coverage calculations
- Visualization tools for coverage plots
- Wrapper scripts for vsearch and ViralMSA

## Installation

### Option 1: Conda (Recommended)

The easiest way to install VICON with all dependencies:

```bash
# Create and activate a new environment
conda create -n vicon python=3.11
conda activate vicon

# Install VICON and all dependencies
conda install -c conda-forge -c bioconda -c eka97 vicon

# Set required permissions
chmod +x "$CONDA_PREFIX/bin/vicon-run"
chmod +x "$CONDA_PREFIX/bin/viralmsa"
chmod +x "$CONDA_PREFIX/bin/minimap2"
```

### Option 2: PyPI (pip)

Install from [PyPI](https://pypi.org/project/vicon/):

```bash
pip install vicon
```

> **Note:** When installing via pip, you must manually install these external dependencies:
> - **minimap2** (≥2.30)
> - **vsearch**
> - **ViralMSA**

#### Installing External Dependencies

**Ubuntu / Debian:**
```bash
sudo apt-get update
sudo apt-get install -y minimap2 vsearch
```

**macOS (Homebrew):**
```bash
brew install minimap2 vsearch
```

**ViralMSA:**
```bash
mkdir -p ~/bin && cd ~/bin
wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py"
chmod +x ViralMSA.py
ln -sf "$PWD/ViralMSA.py" ~/.local/bin/viralmsa
```

## Usage

Run the VICON pipeline with:

```bash
vicon-run --config path/to/your/config.yaml
```

### Input FASTA Preprocessing

> **Note:**  
> VICON automatically preprocesses your input FASTA files (both sample and reference) before analysis:
> - Converts all sequences to uppercase
> - Cleans and standardizes FASTA headers
> - Replaces any non-ATCG characters in sequences with 'N'
>
> You do not need to manually edit or check your FASTA files for these issues.

### Example Configuration

Create a configuration file (`config.yaml`):

```yaml
project_path: "project_path"
virus_name: "orov"
input_sample: "data/orov/samples/samples.fasta"
input_reference: "data/orov/reference/reference.fasta"
email: "email@address.com"
kmer_size: 150
threshold: 147 # shows a tolerance of 150-147 = 3 degenerations
l_gene_start: 8000
l_gene_end: 16000
coverage_ratio: 0.5
min_year: 2020
threshold_ratio: 0.01
drop_old_samples: false
drop_mischar_samples: true
```

### FASTA Header Year Extraction

The pipeline automatically extracts years from FASTA headers using a two-step approach:

1. **Priority extraction**: Years following separators (`|`, `_`, `/`, `-`)
2. **Fallback extraction**: Any standalone 4-digit number between 1850-2030

| Header Example           | Year Extracted? | Extracted Year | Reason                          |
|-------------------------|:--------------:|:--------------:|---------------------------------|
| `>sample|2021`           | ✅ Yes          | 2021           | After pipe separator            |
| `>sample_2020`           | ✅ Yes          | 2020           | After underscore separator      |
| `>sample/2019/data`      | ✅ Yes          | 2019           | After slash separator           |
| `>sample-2022-final`     | ✅ Yes          | 2022           | After dash separator            |
| `>data 2021 sequence`    | ✅ Yes          | 2021           | Standalone 4-digit number       |
| `>sample.2020.version`   | ✅ Yes          | 2020           | Standalone 4-digit number       |
| `>test2021extra`         | ✅ Yes          | 2021           | Standalone 4-digit number       |
| `>sample|202`            | ❌ No           | -              | Not 4 digits                    |
| `>sample_1800_old`       | ❌ No           | -              | Outside valid range (1850-2030) |
| `>sample20213long`       | ❌ No           | -              | 5 consecutive digits            |

> **Best Practice:** Use `|YYYY`, `_YYYY`, `/YYYY`, or `-YYYY` patterns for reliable year extraction.

## License

This project is licensed under the terms of the MIT license.
