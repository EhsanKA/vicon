# VICON Installation and Usage Guide

This guide explains how to install and run VICON version 1.0.0 using Conda, including setting up necessary permissions and executing the workflow with a sample configuration.

## 1. Environment Activation

Activate the dedicated Conda environment for VICON:

```bash
conda activate vicon
```


## 2. VICON Installation

Install VICON version 1.0.0 along with its dependencies from the specified channels:

```bash
conda install -c conda-forge -c bioconda -c eka97 vicon=1.0.0 -y
```

- **conda-forge** and **bioconda** provide bioinformatics packages.
- **eka97** is a custom channel for VICON-specific releases.


## 3. Set Executable Permissions

Ensure the main scripts are executable by running:

```bash
chmod +x "$CONDA_PREFIX/bin/vicon-run"
chmod +x "$CONDA_PREFIX/bin/viralmsa"
chmod +x "$CONDA_PREFIX/bin/minimap2"
```

- This step is crucial for Unix-based systems to allow script execution.


## 4. Running VICON

Execute the VICON pipeline with your desired configuration:

```bash
vicon-run --config configs/config_orov_s.yaml
```

- Replace `configs/config_orov_s.yaml` with your specific configuration file if needed.
- The configuration file should define parameters for your analysis (e.g., reference data, input files, workflow settings).


## Notes

- Ensure you have the required permissions to install packages and modify files in your Conda environment.
- If you encounter errors regarding missing dependencies or permissions, verify your Conda environment is activated and that you have write access to the `$CONDA_PREFIX/bin` directory.
- For detailed configuration options, refer to the documentation included with the VICON package or the sample YAML configuration files.


## Troubleshooting

- **Command not found:** Make sure the Conda environment is activated and the installation completed without errors.
- **Permission denied:** Double-check the `chmod` commands and your user permissions.
- **Dependency issues:** Ensure all specified channels are accessible and up-to-date.

This guide should help you set up and run VICON for your bioinformatics workflows efficiently.

