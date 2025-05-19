#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

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

def main():
    parser = argparse.ArgumentParser(description='Combine multiple FASTA files into one.')
    parser.add_argument('input_dir', help='Directory containing FASTA files')
    parser.add_argument('output_file', help='Path to output combined FASTA file')
    
    args = parser.parse_args()
    
    combine_fasta_files(args.input_dir, args.output_file)

if __name__ == '__main__':
    main() 