#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import os
import re
from pathlib import Path

def extract_segment_number(filename):
    """Extract segment number from filename."""
    match = re.search(r'segment\s*(\d+)', filename, re.IGNORECASE)
    return int(match.group(1)) if match else None

def parse_collection_date(date_str):
    """Parse and standardize collection date."""
    if not date_str:
        return "NA"
    # Remove any text in parentheses
    date_str = re.sub(r'\([^)]*\)', '', date_str).strip()
    return date_str

def extract_sequence_from_origin(record):
    """Extract sequence from ORIGIN section of GenBank record."""
    if not hasattr(record, '_raw') or not record._raw:
        return None
    
    # Find ORIGIN section
    origin_match = re.search(r'ORIGIN\s*(.*?)\/\/', record._raw, re.DOTALL)
    if not origin_match:
        return None
    
    # Extract sequence
    origin_text = origin_match.group(1)
    # Remove numbers and whitespace, join all lines
    sequence = ''.join(re.sub(r'[\d\s]', '', line) for line in origin_text.split('\n'))
    return sequence.upper()

def extract_record_info(record):
    """Extract required information from a GenBank record."""
    # Initialize default values
    info = {
        'accession': 'NA',
        'organism': 'NA',
        'geo_loc_name': 'NA',
        'collection_date': 'NA'
    }
    
    # Get accession (without version)
    if hasattr(record, 'annotations'):
        if 'accessions' in record.annotations and record.annotations['accessions']:
            # Use the first accession number from the annotations
            info['accession'] = record.annotations['accessions'][0]
        else:
            # Fall back to parsing from record.id (removing version number)
            info['accession'] = record.id.split('.')[0]
    
    # Get organism
    if hasattr(record, 'annotations') and 'organism' in record.annotations:
        info['organism'] = record.annotations['organism']
    
    # Extract features from source
    for feature in record.features:
        if feature.type == 'source':
            qualifiers = feature.qualifiers
            if 'country' in qualifiers:
                info['geo_loc_name'] = qualifiers['country'][0]
            if 'collection_date' in qualifiers:
                info['collection_date'] = parse_collection_date(qualifiers['collection_date'][0])
    
    return info

def process_genbank_files(input_dir):
    """Process all GenBank files and organize sequences by strain."""
    # Dictionary to store sequences by accession
    sequences_by_strain = defaultdict(dict)
    
    # Process each GenBank file
    for gb_file in sorted(Path(input_dir).glob('*.gb')):
        segment_num = extract_segment_number(gb_file.name)
        if segment_num is None:
            print(f"Warning: Could not extract segment number from {gb_file.name}")
            continue
            
        print(f"Processing segment {segment_num} from {gb_file.name}")
        
        # Parse the GenBank file
        for record in SeqIO.parse(gb_file, "genbank"):
            info = extract_record_info(record)
            
            # Try to get sequence from record.seq first
            try:
                sequence = str(record.seq).upper()
            except:
                # If that fails, try to extract from ORIGIN section
                sequence = extract_sequence_from_origin(record)
                if not sequence:
                    print(f"Warning: Could not extract sequence for {record.id} in segment {segment_num}")
                    continue
            
            # Create a unique identifier for this strain
            strain_key = f"{info['accession']}"
            
            # Store the sequence and information
            sequences_by_strain[strain_key][segment_num] = {
                'sequence': sequence,
                'info': info
            }
    
    return sequences_by_strain

def write_concatenated_sequences(sequences_by_strain, output_file):
    """Write concatenated sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for strain_key, segments in sequences_by_strain.items():
            # Only process strains that have all segments
            if len(segments) != 11:  # Rotavirus A has 11 segments
                print(f"Warning: Strain {strain_key} has only {len(segments)} segments, skipping")
                continue
            
            # Get info from any segment (they should all be the same)
            info = next(iter(segments.values()))['info']
            
            # Create header
            header = f">{info['accession']}|{','.join(str(i) for i in range(1,12))}|{info['organism']}|{info['geo_loc_name']}|{info['collection_date']}"
            
            # Concatenate sequences in order of segments
            concatenated_seq = ''.join(segments[i]['sequence'] for i in range(1, 12))
            
            # Write to file
            f.write(f"{header}\n")
            
            # Write sequence in lines of 80 characters
            for i in range(0, len(concatenated_seq), 80):
                f.write(concatenated_seq[i:i+80] + '\n')

def main():
    # Set input and output paths
    input_dir = "/fast/AG_Ohler/ekarimi/projects/vicon/data/rotav/Human_Rotavirus_A_1"
    output_dir = os.path.join(os.path.dirname(input_dir), "concatenated")
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, "concatenated_sequences.fasta")
    
    print("Processing GenBank files...")
    sequences_by_strain = process_genbank_files(input_dir)
    
    print("Writing concatenated sequences...")
    write_concatenated_sequences(sequences_by_strain, output_file)
    
    print(f"Done! Output written to {output_file}")

if __name__ == "__main__":
    main() 