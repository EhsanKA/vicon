import subprocess
from pathlib import Path

def run_vsearch(input_fasta, output_fasta, cluster_output, sizeout=True, logger=None):
    """
    Runs the `vsearch` command to process a fasta file and displays logs in real-time in Jupyter Notebook.
    
    Args:
        input_fasta (str): Path to the input fasta file.
        output_fasta (str): Path to save the dereplicated fasta file.
        cluster_output (str): Path to save the cluster file.
        sizeout (bool): Whether to include the `--sizeout` flag.
        logger (logging.Logger, optional): Logger to use for logging output. If None, prints to stdout.
    """
    cmd = [
        "vsearch",
        "--fastx_uniques", input_fasta,
        "--fastaout", output_fasta,
        "--uc", cluster_output
    ]
    if sizeout:
        cmd.append("--sizeout")

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    
    # Stream logs line by line
    try:
        for line in process.stdout:
            if logger:
                logger.info(line.strip())
            else:
                print(line, end='')
    except Exception as e:
        if logger:
            logger.error(f"Error reading process output: {e}")
        else:
            print(f"Error reading process output: {e}")
    
    # Wait for the process to complete
    process.wait()
    
    if process.returncode != 0:
        error_msg = f"Vsearch failed with return code {process.returncode}."
        if logger:
            logger.error(error_msg)
        raise RuntimeError(error_msg)

