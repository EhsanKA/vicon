#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict, Counter
import os
import re
from pathlib import Path

def extract_segment_number(filename):
    """Extract segment number from filename."""
    match = re.search(r'segment\s*(\d+)', filename, re.IGNORECASE)
    return int(match.group(1)) if match else None

def get_accession(record):
    """Get accession number without version."""
    if hasattr(record, 'annotations'):
        if 'accessions' in record.annotations and record.annotations['accessions']:
            return record.annotations['accessions'][0]
    return record.id.split('.')[0]

def analyze_segment_distribution(input_dir):
    """Analyze how many accessions have how many segments."""
    # Dictionary to store segments for each accession
    accession_segments = defaultdict(set)
    
    # Process each GenBank file
    for gb_file in sorted(Path(input_dir).glob('*.gb')):
        segment_num = extract_segment_number(gb_file.name)
        if segment_num is None:
            print(f"Warning: Could not extract segment number from {gb_file.name}")
            continue
            
        print(f"Processing segment {segment_num} from {gb_file.name}")
        
        # Parse the GenBank file
        for record in SeqIO.parse(gb_file, "genbank"):
            accession = get_accession(record)
            accession_segments[accession].add(segment_num)
    
    # Count how many accessions have each number of segments
    segment_counts = Counter(len(segments) for segments in accession_segments.values())
    
    # Print results
    print("\nSegment distribution:")
    print("Number of segments | Number of accessions")
    print("-" * 40)
    for num_segments in sorted(segment_counts.keys()):
        print(f"{num_segments:^16} | {segment_counts[num_segments]:^19}")
    
    # Print accessions with all 11 segments
    complete_accessions = [acc for acc, segments in accession_segments.items() 
                         if len(segments) == 11]
    
    print(f"\nFound {len(complete_accessions)} accessions with all 11 segments")
    
    # Print accessions with only one segment
    single_segment_accessions = [acc for acc, segments in accession_segments.items() 
                               if len(segments) == 1]
    
    print(f"Found {len(single_segment_accessions)} accessions with only 1 segment")
    
    # Save detailed results to a file
    output_dir = os.path.join(os.path.dirname(input_dir), "analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    with open(os.path.join(output_dir, "segment_distribution.txt"), 'w') as f:
        f.write("Accession\tNumber of segments\tSegment numbers\n")
        for accession, segments in sorted(accession_segments.items()):
            segments_str = ','.join(str(s) for s in sorted(segments))
            f.write(f"{accession}\t{len(segments)}\t{segments_str}\n")
    
    return accession_segments, segment_counts

def main():
    # Set input path
    input_dir = "/fast/AG_Ohler/ekarimi/projects/vicon/data/rotav/Human_Rotavirus_A_1"
    
    print("Analyzing segment distribution...")
    accession_segments, segment_counts = analyze_segment_distribution(input_dir)
    
    print("\nDetailed results have been saved to 'segment_distribution.txt'")

if __name__ == "__main__":
    main() 