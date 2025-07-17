import numpy as np
import pandas as pd
import yaml
import os



def compare_sequences(reference, sequence):
    """
    Compares a subsequence to the reference sequence and returns an alignment string.
    Matching positions are replaced with '-', and non-matching positions retain the sequence's character.

    Parameters:
    - reference: The reference subsequence (most frequent).
    - sequence: The sequence to compare.

    Returns:
    - comparison: A string showing the comparison result.
    """
    # If the sequence is the same as the reference, return the sequence itself (no dashes)
    if reference == sequence:
        return reference

    return ''.join('-' if ref_char == seq_char else seq_char for ref_char, seq_char in zip(reference, sequence))


def count_changes(reference, sequence):
    """
    Counts the number of changes between the reference sequence and the given sequence,
    ignoring dashes ('-') in either sequence.

    Parameters:
    - reference: The reference subsequence (most frequent).
    - sequence: The sequence to compare.

    Returns:
    - changes: The number of differences between the sequences where both have valid characters (ignoring dashes).
    """
    changes = 0
    for ref_char, seq_char in zip(reference, sequence):
        if ref_char != '-' and seq_char != '-' and ref_char != seq_char:
            changes += 1
    return changes


def count_non_gap_characters_from_dataframe(df, sequence_column='alignment'):
    """
    Counts the number of non-gap ('-') characters at each position across all sequences in a DataFrame column.

    Parameters:
    df (pd.DataFrame): The DataFrame containing sequences.
    sequence_column (str): The name of the column containing sequences.

    Returns:
    pd.DataFrame: A DataFrame with positions as the index and counts as the values.
    """
    if sequence_column not in df.columns:
        raise ValueError(f"Column '{sequence_column}' does not exist in the DataFrame.")
    
    sequences = df[sequence_column].tolist()
    
    if not sequences:
        raise ValueError("The list of sequences is empty.")
    
    sequence_length = len(sequences[0])
    
    # Ensure all sequences are of the same length
    for seq in sequences:
        if len(seq) != sequence_length:
            raise ValueError("All sequences must be of the same length.")
    
    # Convert sequences to a 2D NumPy array
    seq_array = np.array([list(seq) for seq in sequences])
    
    # Create a boolean array where True indicates a non-gap character
    non_gap_array = seq_array != '-'
    
    # Sum over the sequences to get counts at each position
    counts = non_gap_array.sum(axis=0)
    
    # Create a DataFrame from the counts
    df_counts = pd.DataFrame({
        'Position': range(1, sequence_length + 1),
        'NonGapCount': counts
    })
    
    df_counts.set_index('Position', inplace=True)
    return df_counts

def find_min_coverage_threshold(df, coverage_ratio=0.5, logger=None):
    """
    Calculates the minimum coverage threshold based on the coverage ratio.
    """
    min_coverage = int(df.shape[0] * coverage_ratio)
    if logger:
        logger.info(f"Minimum coverage threshold set to {min_coverage} based on coverage ratio {coverage_ratio}")
    else:
        print(f"Minimum coverage threshold set to {min_coverage} based on coverage ratio {coverage_ratio}")
    return min_coverage


import tempfile
import shutil

def replace_hyphen_with_n(input_fasta, output_fasta):
    """
    Replaces all '-' characters with 'N' in sequences from a FASTA file.
    The processed content is first saved to a temporary file, 
    then moved to overwrite the original file.
    
    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to save the processed FASTA file.
    """
    try:
        # Create a temporary file to store modified content
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_file:
            temp_file_name = temp_file.name
            
            with open(input_fasta, "r") as infile:
                for line in infile:
                    if line.startswith(">"):
                        # Write header lines as-is
                        temp_file.write(line)
                    else:
                        # Replace '-' with 'N' in sequence lines and write
                        modified_sequence = line.replace("-", "N")
                        temp_file.write(modified_sequence)

        # Overwrite the original file or save to the output path
        shutil.move(temp_file_name, output_fasta)
        print(f"Processed FASTA file saved to: {output_fasta}")
        
    except Exception as e:
        print(f"An error occurred: {e}")



from Bio import SeqIO
import sys

def check_fasta_non_atcg(fasta_path):
    valid_chars = set("ATCG")
    
    with open(fasta_path, "r") as fasta_file:
        lines = fasta_file.readlines()
        for i in range(0, len(lines), 2):  # Process header and sequence pairs
            header = lines[i].strip()
            sequence = lines[i+1].strip().upper() if i+1 < len(lines) else ""
            invalid_chars = set(sequence) - valid_chars
            if invalid_chars:
                print(f"{header} contains non-ATCG characters: {invalid_chars}")
                break
            else:
                print(f"{header} is valid.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_fasta.py <fasta_file>")
    else:
        check_fasta_non_atcg(sys.argv[1])


import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

def process_fasta_content(content):
    """
    Process FASTA content:
    - Convert sequences to uppercase
    - Replace spaces and tabs in headers with underscores
    
    Args:
        content (str): Raw FASTA content
    Returns:
        str: Processed FASTA content
    """
    lines = content.split('\n')
    processed_lines = []
    
    for line in lines:
        if line.startswith('>'):  # Header line
            # Replace spaces and tabs with underscores
            processed_lines.append(line.replace(' ', '_').replace('\t', '_'))
        else:  # Sequence line
            # Convert sequence to uppercase
            processed_lines.append(line.upper())
    
    return '\n'.join(processed_lines)

def combine_fasta_files(input_dir, output_file):
    """
    Combine multiple FASTA files from a directory into a single FASTA file.
    
    Args:
        input_dir (str): Path to directory containing FASTA files
        output_file (str): Path to output combined FASTA file
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Get all .fasta files from input directory
    fasta_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.fasta')])
    
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}")
        return
    
    print(f"Found {len(fasta_files)} FASTA files")
    
    # Combine files
    with open(output_file, 'w') as outfile:
        for fasta_file in fasta_files:
            input_path = os.path.join(input_dir, fasta_file)
            print(f"Processing: {fasta_file}")
            
            with open(input_path, 'r') as infile:
                # Read content, process it, and write to output file
                content = infile.read()
                processed_content = process_fasta_content(content)
                outfile.write(processed_content)
                # Add a newline between files if needed
                outfile.write('\n')
    
    print(f"Combined FASTA file saved to: {output_file}")


from Bio import SeqIO
from Bio.Seq import Seq

def process_fasta_file(input_file, output_file):
    """
    Process a FASTA file:
    - Sequence: uppercased, non-ATCG replaced with N
    - Header: preserve the entire original header line
    """
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            # Preserve the full header (including spaces)
            record.id = record.description  # This keeps the full header after '>'
            record.id = record.id.replace(" ", "_")
            record.description = ""         # Prevents Biopython from appending description again

            # Clean the sequence
            record.seq = record.seq.replace("\n", "")  # Remove any newlines
            seq_str = ''.join(['N' if c.upper() not in 'ATCG' else c.upper() for c in str(record.seq)])
            record.seq = Seq(seq_str)

            SeqIO.write(record, out_handle, "fasta")

def filter_by_most_common_kmers(df):
    kmer1_most = df['kmer1'].mode()[0]
    kmer2_most = df['kmer2'].mode()[0]
    kmer1_count = df['kmer1'].value_counts().get(kmer1_most, 0)
    kmer2_count = df['kmer2'].value_counts().get(kmer2_most, 0)
    filtered_df = df[(df['kmer1'] == kmer1_most) | (df['kmer2'] == kmer2_most)].copy()
    return filtered_df, kmer1_most, kmer2_most, kmer1_count, kmer2_count

def create_results_dir(config):
    """
    Creates the results directory structure for the pipeline.
    
    Args:
        config (dict): Configuration dictionary containing project settings
        
    Returns:
        str: Path to the results directory
    """
    base_path = config["project_path"]
    virus = config["virus_name"]
    results_dir = os.path.join(base_path, "results", virus)
    
    # Create the results directory if it doesn't exist
    os.makedirs(results_dir, exist_ok=True)
    
    # Create subdirectories
    os.makedirs(os.path.join(results_dir, "logs"), exist_ok=True)
    
    return results_dir

def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
