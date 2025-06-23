# VICON Installation and Usage Guide (Mac)

## 1. Installation Guideline (One-Time Setup)

### A. Install Miniconda (if not already installed)

1. Go to: https://docs.conda.io/en/latest/miniconda.html  
2. Download the installer for MacOS (choose Apple Silicon or Intel, as appropriate).
3. Open Terminal and run the downloaded installer, following the prompts.

### B. Create and Activate a New Environment

```bash
conda create -n vicon python=3.11 -y
conda activate vicon
```

### C. Install VICON

Run the following command to install VICON and all dependencies:

```bash
conda install -c conda-forge -c bioconda -c eka97 vicon
```

### D. Set Permissions (if needed)

```bash
chmod +x "$CONDA_PREFIX/bin/vicon-run"
chmod +x "$CONDA_PREFIX/bin/viralmsa"
chmod +x "$CONDA_PREFIX/bin/minimap2"
```

---

## 2. Usage Guideline (Every Time You Want to Run VICON)

### A. Open Terminal

- You can find Terminal in Applications > Utilities > Terminal.

### B. Activate the VICON Environment

```bash
conda activate vicon
```

### C. Prepare Your Input Files

- Make sure your input FASTA files and config file are in the correct locations as specified in your config YAML.
- If you just downloaded or edited files, double-check their paths.

### D. Run the Pipeline

```bash
vicon-run --config path/to/your/config.yaml
```

- The pipeline will automatically preprocess your FASTA files (uppercase, clean headers, replace non-ATCG with N).
- Results and logs will be saved in the output directory specified in your config.

### E. When Finished

To leave the environment:

```bash
conda deactivate
```

---

## 3. Sample Config File

Save this as `config.yaml` and edit the paths as needed:

```yaml
project_path: "/Users/yourusername/vicon_project"
virus_name: "orov"
input_sample: "data/orov/samples/samples.fasta"
input_reference: "data/orov/reference/reference.fasta"
email: "your.email@address.com"
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

### Directory Structure Example

```
/Users/yourusername/vicon_project/
├── data/
│   └── orov/
│       ├── samples/
│       │   └── samples.fasta
│       └── reference/
│           └── reference.fasta
├── config.yaml
``` 