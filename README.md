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

### Standard Installation

1. Create and activate a conda environment:
   ```bash
   conda create -n vicon python=3.11
   conda activate vicon
   ```

2. Install VICON and its dependencies:
   ```bash
   conda install -c conda-forge -c bioconda -c eka97 vicon
   ```

3. Set required permissions:
   ```bash
   chmod +x "$CONDA_PREFIX/bin/vicon-run"
   chmod +x "$CONDA_PREFIX/bin/viralmsa"
   chmod +x "$CONDA_PREFIX/bin/minimap2"
   ```

<!-- ### Development Installation

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
   - Install vsearch:
     ```bash
     conda install -c bioconda vsearch -y
     ```
   - ViralMSA:
      ```bash
      mkdir -p scripts && cd scripts
      wget "https://github.com/niemasd/ViralMSA/releases/latest/download/ViralMSA.py"
      chmod a+x ViralMSA.py
      cd ../
      ```

4. Install VICON in development mode:
   ```bash
   pip install -e .
   ```

5. Set required permissions:
   ```bash
   chmod +x "$CONDA_PREFIX/bin/vicon-run"
   chmod +x "$CONDA_PREFIX/bin/viralmsa"
   ``` -->

## Usage

To run the VICON pipeline, use the following command:

```bash
vicon-run --config path/to/your/config.yaml
```

### Input FASTA Preprocessing

> **Note:**  
> When you run the pipeline, VICON will automatically preprocess your input FASTA files (both sample and reference) before any analysis.  
> This step:
> - Converts all sequences to uppercase
> - Cleans and standardizes FASTA headers
> - Replaces any non-ATCG characters in sequences with 'N'
>
> The cleaned files are used for all downstream analysis, so you do not need to manually edit or check your FASTA files for these issues.

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

### FASTA Header Year Extraction: Supported Formats

The pipeline will extract the year from your FASTA headers if they match one of these patterns:

| Header Example           | Year Extracted? | Extracted Year | Reason                |
|-------------------------|:--------------:|:--------------:|-----------------------|
| `>sample|2021`           | Yes            | 2021           | After pipe            |
| `>sample_2020`           | Yes            | 2020           | After underscore      |
| `>sample|15-JAN-2019`    | Yes            | 2019           | Date format           |
| `>sample_23-DEC-2022`    | Yes            | 2022           | Date format           |
| `>sample_ABC2021`        | Yes            | 2021           | Ends with 4 digits    |
| `>sample_2024_04_15`     | Yes            | 2024           | First _YYYY           |
| `>sample|2024_04_15`     | Yes            | 2024           | First |YYYY            |
| `>sample2021extra`       | No             | -              | Not at end/after sep  |
| `>sample|202`            | No             | -              | Not 4 digits          |

> **Tip:** For best results, use `|YYYY` or `_YYYY` at the end of your FASTA header.

## License
This project is licensed under the terms of the MIT license.
