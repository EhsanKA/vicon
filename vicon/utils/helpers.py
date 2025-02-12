import numpy as np
import pandas as pd



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

def find_min_coverage_threshold(df, coverage_ratio=0.5):
    """
    Calculates the minimum coverage threshold based on the coverage ratio.
    """
    min_coverage = int(df.shape[0] * coverage_ratio)
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
