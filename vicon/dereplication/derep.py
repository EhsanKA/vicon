import subprocess
from pathlib import Path

import subprocess

def run_vsearch(input_fasta, output_fasta, cluster_output, sizeout=True):
    """
    Runs the `vsearch` command to process a fasta file and displays logs in real-time in Jupyter Notebook.
    
    Args:
        input_fasta (str): Path to the input fasta file.
        output_fasta (str): Path to save the dereplicated fasta file.
        cluster_output (str): Path to save the cluster file.
        sizeout (bool): Whether to include the `--sizeout` flag.
    """
    cmd = [
        "vsearch",
        "--fastx_uniques", input_fasta,
        "--fastaout", output_fasta,
        "--uc", cluster_output
    ]
    if sizeout:
        cmd.append("--sizeout")

    # Use Popen to stream logs to Jupyter in real-time
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    
    # Stream logs line by line
    try:
        for line in process.stdout:
            print(line, end='')  # Print to Jupyter Notebook output
    except Exception as e:
        print(f"Error reading process output: {e}")
    
    # Wait for the process to complete
    process.wait()
    
    if process.returncode != 0:
        raise RuntimeError(f"Vsearch failed with return code {process.returncode}.")

