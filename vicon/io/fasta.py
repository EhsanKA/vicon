import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def read_fasta(file_path):
    """Reads a FASTA file and returns the sequence as a string."""
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record.seq)

def read_fasta_to_dataframe(fasta_file):
    """Reads a FASTA file and returns a DataFrame with sequence IDs and sequences."""
    sequences = []
    sequence_id = ''
    sequence = ''
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append([sequence_id, sequence])
                    sequence = ''
                sequence_id = line[1:]
            else:
                sequence += line
        if sequence:
            sequences.append([sequence_id, sequence])
    return pd.DataFrame(sequences, columns=['ID', 'Sequence'])


# def create_folders_and_save_sequences(fasta_file, new_address= 'a'):
#     with open(fasta_file, 'r') as file:
#         current_folder = None
#         for line in file:
#             if line.startswith('>'):
#                 header = line.strip().lstrip('>')
#                 folder_name = header.split('|')[0].replace("/", "_")
#                 current_folder = new_address+"/"+folder_name
#                 os.makedirs(current_folder, exist_ok=True)
#                 sequence_file = os.path.join(current_folder, "genome.fasta")
#                 with open(sequence_file, 'w') as seq_file:
#                     seq_file.write(line)
#             else:
#                 with open(sequence_file, 'a') as seq_file:
#                     seq_file.write(line)

import os

def create_folders_and_save_sequences(fasta_file, new_address='a'):
    with open(fasta_file, 'r') as file:
        current_folder = None
        sequence_file = None

        for line in file:
            if line.startswith('>'):
                # Process header line
                header = line.strip().lstrip('>')
                folder_name = header.split('|')[0].replace("/", "_")
                current_folder = os.path.join(new_address, folder_name)
                os.makedirs(current_folder, exist_ok=True)

                # Initialize sequence file for the current folder
                sequence_file = os.path.join(current_folder, "genome.fasta")
                with open(sequence_file, 'w') as seq_file:
                    seq_file.write(line)  # Write the header line
            else:
                # Write the sequence line to the last initialized sequence_file
                if sequence_file:
                    with open(sequence_file, 'a') as seq_file:
                        seq_file.write(line)



def remove_first_record(input_fasta, output_fasta, logger=None):
    """
    Removes the first record (header and sequence) from a FASTA file.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to save the output FASTA file without the first record.

    Returns:
        None
    """
    # Read all records from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    if len(records) <= 1:
        raise ValueError("The FASTA file must contain at least two records to remove the first one.")
    
    # Write the remaining records to the output FASTA file
    SeqIO.write(records[1:], output_fasta, "fasta",)
    if logger:
        logger.info(f"The first record has been removed. Updated FASTA saved to: {output_fasta}")
    else:
        print(f"The first record has been removed. Updated FASTA saved to: {output_fasta}")


def generate_remaining_fasta(df3, df_samples, kmer1, kmer2, fasta_file_path, output_path, logger=None):
    """
    Generate a FASTA file containing sequences that are not covered by the degenerate kmers.
    
    Args:
        df3 (pd.DataFrame): Binary matrix with kmer coverage data
        df_samples (pd.DataFrame): DataFrame containing sample information with 'ID' column
        kmer1 (int): Position of first kmer
        kmer2 (int): Position of second kmer  
        fasta_file_path (str): Path to the original FASTA file to extract sequences from
        output_path (str): Path where the remaining sequences FASTA file will be saved
        logger (logging.Logger, optional): Logger for logging information
    
    Returns:
        int: Number of remaining sequences written to file
    """
    if kmer1 not in df3.columns or kmer2 not in df3.columns:
        if logger:
            logger.warning(f"Kmer1 ({kmer1}) or kmer2 ({kmer2}) not found in df3 columns")
        return 0
    
    # Get common indices between df3 and df_samples
    common_idx = df3.index.intersection(df_samples['ID'])
    deg_df = df3.loc[common_idx, [kmer1, kmer2]]
    
    if logger:
        logger.info(f"Total samples in df3: {len(df3)}")
        logger.info(f"Total samples in df_samples: {len(df_samples)}")
        logger.info(f"Common samples between df3 and df_samples: {len(common_idx)}")
        logger.info(f"Samples covered by kmer1: {deg_df[kmer1].sum()}")
        logger.info(f"Samples covered by kmer2: {deg_df[kmer2].sum()}")
    
    # Find samples that are NOT covered by either kmer (sum == 0)
    not_covered_mask = (deg_df.sum(axis=1) == 0)
    not_covered_ids = deg_df[not_covered_mask].index.tolist()
    
    if logger:
        logger.info(f"Found {len(not_covered_ids)} sequences not covered by degenerate kmers")
        if len(not_covered_ids) > 0:
            logger.info(f"Not covered IDs: {not_covered_ids[:5]}...")  # Show first 5
    
    # Read the original FASTA file
    original_fasta_df = read_fasta_to_dataframe(fasta_file_path)
    
    # Filter sequences that are not covered
    remaining_sequences = original_fasta_df[original_fasta_df['ID'].isin(not_covered_ids)]
    
    # Write remaining sequences to new FASTA file manually
    with open(output_path, 'w') as f:
        for _, row in remaining_sequences.iterrows():
            f.write(f">{row['ID']}\n")
            f.write(f"{row['Sequence']}\n")
    
    if logger:
        logger.info(f"Remaining sequences saved to: {output_path}")
    
    return len(remaining_sequences)
