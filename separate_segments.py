#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import os
import re
from pathlib import Path
import warnings
from Bio import BiopythonWarning

# Filter out Biopython warnings
warnings.filterwarnings('ignore', category=BiopythonWarning)
warnings.filterwarnings('ignore', category=UserWarning)

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
    
    try:
        # Find ORIGIN section
        origin_match = re.search(r'ORIGIN\s*(.*?)\/\/', record._raw, re.DOTALL)
        if not origin_match:
            return None
        
        # Extract sequence
        origin_text = origin_match.group(1)
        # Remove numbers and whitespace, join all lines
        sequence = ''.join(re.sub(r'[\d\s]', '', line) for line in origin_text.split('\n'))
        return sequence.upper() if sequence else None
    except Exception as e:
        print(f"Error extracting sequence from ORIGIN for {record.id}: {str(e)}")
        return None

def extract_record_info(record):
    """Extract required information from a GenBank record."""
    # Initialize default values
    info = {
        'accession': 'NA',
        'organism': 'NA',
        'geo_loc_name': 'NA',
        'collection_date': 'NA'
    }
    
    try:
        # Get accession (without version)
        if hasattr(record, 'annotations'):
            if 'accessions' in record.annotations and record.annotations['accessions']:
                info['accession'] = record.annotations['accessions'][0]
            else:
                info['accession'] = record.id.split('.')[0]
        
        # Get organism
        if hasattr(record, 'annotations') and 'organism' in record.annotations:
            info['organism'] = record.annotations['organism']
        
        # Extract features from source
        for feature in record.features:
            if feature.type == 'source':
                qualifiers = feature.qualifiers
                
                # First check for geo_loc_name directly
                if 'geo_loc_name' in qualifiers:
                    info['geo_loc_name'] = qualifiers['geo_loc_name'][0]
                # Fallback to other location fields if geo_loc_name is not present
                elif 'country' in qualifiers:
                    info['geo_loc_name'] = qualifiers['country'][0]
                elif 'isolation_source' in qualifiers:
                    info['geo_loc_name'] = qualifiers['isolation_source'][0]
                elif 'lat_lon' in qualifiers:
                    info['geo_loc_name'] = qualifiers['lat_lon'][0]
                elif 'note' in qualifiers and 'location' in qualifiers['note'][0].lower():
                    info['geo_loc_name'] = qualifiers['note'][0]
                
                if 'collection_date' in qualifiers:
                    info['collection_date'] = parse_collection_date(qualifiers['collection_date'][0])
    except Exception as e:
        print(f"Error extracting info for {record.id}: {str(e)}")
    
    return info

def process_genbank_files(input_dir, base_output_dir):
    """Process GenBank files and create separate FASTA files for each segment."""
    # Dictionary to store file handles and statistics
    segment_files = {}
    stats = defaultdict(lambda: {'processed': 0, 'failed': 0})
    
    try:
        # Process each GenBank file
        for gb_file in sorted(Path(input_dir).glob('*.gb')):
            segment_num = extract_segment_number(gb_file.name)
            if segment_num is None:
                print(f"Warning: Could not extract segment number from {gb_file.name}")
                continue
                
            print(f"\nProcessing segment {segment_num} from {gb_file.name}")
            
            # Create segment-specific directory and file
            if segment_num not in segment_files:
                # Create directory for this segment
                segment_dir = os.path.join(base_output_dir, f"segment_{segment_num}")
                os.makedirs(segment_dir, exist_ok=True)
                
                # Create file handle in the segment directory
                output_file = os.path.join(segment_dir, "sequences.fasta")
                segment_files[segment_num] = open(output_file, 'w')
            
            # Parse the GenBank file
            for record in SeqIO.parse(gb_file, "genbank"):
                stats[segment_num]['processed'] += 1
                info = extract_record_info(record)
                
                # Try to get sequence from record.seq first
                sequence = None
                try:
                    sequence = str(record.seq).upper()
                except:
                    # If that fails, try to extract from ORIGIN section
                    sequence = extract_sequence_from_origin(record)
                
                if not sequence:
                    print(f"Warning: Could not extract sequence for {record.id} in segment {segment_num}")
                    stats[segment_num]['failed'] += 1
                    continue
                
                # Create header
                header = f">{info['accession']}|{segment_num}|{info['organism']}|{info['geo_loc_name']}|{info['collection_date']}"
                
                # Write to appropriate segment file
                segment_files[segment_num].write(f"{header}\n")
                
                # Write sequence in lines of 80 characters
                for i in range(0, len(sequence), 80):
                    segment_files[segment_num].write(sequence[i:i+80] + '\n')
    
    finally:
        # Close all file handles
        for fh in segment_files.values():
            fh.close()
    
    return len(segment_files), stats

def main():
    # Set input and output paths
    input_dir = "/fast/AG_Ohler/ekarimi/projects/vicon/data/rotav/Human_Rotavirus_A_1"
    # Create a base directory for all segments
    base_output_dir = os.path.join(input_dir, "segments")
    
    print(f"Processing GenBank files...")
    print(f"Base output directory: {base_output_dir}")
    num_segments, stats = process_genbank_files(input_dir, base_output_dir)
    
    print(f"\nProcessing Summary:")
    print("-" * 50)
    for segment_num in sorted(stats.keys()):
        processed = stats[segment_num]['processed']
        failed = stats[segment_num]['failed']
        success = processed - failed
        print(f"Segment {segment_num}:")
        print(f"  - Processed: {processed}")
        print(f"  - Successful: {success}")
        print(f"  - Failed: {failed}")
    
    print(f"\nCreated {num_segments} segment directories in {base_output_dir}")
    print("\nFiles created:")
    total_size = 0
    for i in range(1, 12):
        segment_dir = os.path.join(base_output_dir, f"segment_{i}")
        fasta_file = os.path.join(segment_dir, "sequences.fasta")
        if os.path.exists(fasta_file):
            size = os.path.getsize(fasta_file)
            total_size += size
            print(f"segment_{i}/sequences.fasta: {size/1024:.2f} KB")
    print(f"\nTotal size of all FASTA files: {total_size/1024/1024:.2f} MB")

if __name__ == "__main__":
    main() 