import subprocess
from pathlib import Path

def run_viralmsa(sample_fasta, output_dir, reference_fasta, email='email@address.com'):
    """
    Runs the ViralMSA script for multiple sequence alignment.
    
    Args:
        email (str): Email address for the tool.
        sample_fasta (str): Path to the sample fasta file.
        output_dir (str): Directory to save the alignment.
        reference_fasta (str): Path to the reference fasta file.
    """

    cmd = [
        "viralmsa",
        "-e", email,
        "-s", sample_fasta,
        "-o", output_dir,
        "-r", reference_fasta,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"ViralMSA failed:\n{result.stderr}")
    return result
