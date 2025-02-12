import pandas as pd
from Bio import SeqIO
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



def remove_first_record(input_fasta, output_fasta):
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
    SeqIO.write(records[1:], output_fasta, "fasta")
    print(f"The first record has been removed. Updated FASTA saved to: {output_fasta}")
