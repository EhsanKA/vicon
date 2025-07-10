import pandas as pd
import numpy as np
from collections import Counter
from vicon.utils.year_extraction import get_subsequences_with_years, get_sorted_subsequences, simulate_year_bars
from vicon.utils.helpers import compare_sequences, count_changes



def compare_kmers(kmer, reference):
    """
    Compares a kmer with a reference sequence and returns a comparison string.
    Matching positions are replaced with '-', and non-matching positions retain the kmer's character.
    """
    result = []
    for i in range(len(kmer)):
        if kmer[i] == reference[i]:
            result.append('-')
        else:
            result.append(kmer[i])
    return ''.join(result)

def filter_kmers_by_changes(df_samples, kmer_name='kmer1', num_changes_threshold=3):
    """
    Filters the DataFrame to return rows where the subsequence has more than a specified number of changes compared to the
    most frequent subsequence and adds a new column showing the number of changes.
    
    Parameters:
    - df_samples: pandas DataFrame with columns ['ID', 'kmer1'].
    - kmer_name: Name of the kmer column (default is 'kmer1').
    - num_changes_threshold: The minimum number of changes to filter by (default is 3).
    
    Returns:
    - DataFrame with rows having more than `num_changes_threshold` changes and a new column 'num_changes' showing the count.
    """
    # Drop rows where 'Subsequence' (kmer1) is NaN
    df_samples_clean = df_samples[['ID', kmer_name]].dropna()

    # Count frequency of each subsequence
    subseq_freq = Counter(df_samples_clean[kmer_name])

    # Sort subsequences by frequency in descending order
    sorted_subsequences = sorted(subseq_freq.items(), key=lambda x: x[1], reverse=True)
    
    # Get the most frequent subsequence
    most_frequent_subseq = sorted_subsequences[0][0]

    # Create a list to hold subsequences with more than the threshold number of changes
    subseq_with_changes = []

    # Iterate through the subsequences and calculate the number of changes
    for idx, row in df_samples_clean.iterrows():
        subseq = row[kmer_name]
        comparison = compare_kmers(subseq, most_frequent_subseq)
        num_changes = len([char for char in comparison if char != '-'])  # Count only the changed characters
        
        # If the number of changes is greater than the threshold, add it to the list
        if num_changes > num_changes_threshold:
            subseq_with_changes.append({
                'ID': row['ID'],
                kmer_name: subseq,
                'num_changes': num_changes
            })

    # Convert the list of results into a DataFrame
    df_with_changes = pd.DataFrame(subseq_with_changes)

    return df_with_changes

def process_kmers_to_dataframe(df_samples, kmer_column):
    """
    Processes kmers and creates a DataFrame with relevant information about each subsequence.

    Parameters:
    - df_samples: Pandas DataFrame containing the sequence data.
    - kmer_column: Name of the column containing kmers.

    Returns:
    - df_kmers: A Pandas DataFrame with Seq_ID, year range, percentage, alignment, and changes.
    """
    if df_samples.empty:
        raise ValueError("Input DataFrame is empty")
    
    if kmer_column not in df_samples.columns:
        raise ValueError(f"Column '{kmer_column}' not found in DataFrame")
    
    subseq_to_years, subseq_freq = get_subsequences_with_years(df_samples, kmer_column)
    sorted_subsequences = get_sorted_subsequences(subseq_freq)
    
    if not sorted_subsequences:
        raise ValueError(f"No valid subsequences found in column '{kmer_column}'. Check if the data contains valid sequences and years.")
    
    most_frequent_subseq = sorted_subsequences[0][0]
    total_subsequences = sum(subseq_freq.values())

    return create_kmers_dataframe(sorted_subsequences, subseq_to_years, most_frequent_subseq, total_subsequences)


def create_kmers_dataframe(sorted_subsequences, subseq_to_years, most_frequent_subseq, total_subsequences):
    """
    Creates a DataFrame with Seq_ID, year range, percentage, alignment, and changes.

    Parameters:
    - sorted_subsequences: List of sorted subsequences with their frequencies.
    - subseq_to_years: Dictionary mapping subsequences to years.
    - most_frequent_subseq: The most frequent subsequence (reference).
    - total_subsequences: Total number of subsequences (for percentage calculation).

    Returns:
    - df_kmers: A Pandas DataFrame containing the processed subsequences and associated data.
    """
    data = []
    for i, (subseq, freq) in enumerate(sorted_subsequences, start=1):
        percentage = (freq / total_subsequences) * 100 
        # percentage = freq # just for testing
        years = subseq_to_years[subseq]
        s_year = min(years) if years else None
        e_year = max(years) if years else None
        # print(f"Subsequence {i}: {subseq} ({percentage:.2f}%)")
        # print(f"Start Year: {s_year}, End Year: {e_year}")
        # print("*******************************")
        year_range = simulate_year_bars(years)

        # Handle reference sequence alignment and changes
        if subseq == most_frequent_subseq:
            alignment = subseq  # The reference sequence is shown as is (no dashes)
            changes = 0  # No changes for the reference sequence
        else:
            alignment = compare_sequences(most_frequent_subseq, subseq)
            changes = count_changes(most_frequent_subseq, subseq)

        data.append({
            "Seq_ID": i,
            "s_year": s_year,
            "e_year": e_year,
            "year_range": year_range,
            "count": freq,
            "percentage": f"{percentage:.2f}%",
            "changes": changes,
            "alignment": alignment,
        })

    return pd.DataFrame(data)

def mask_kmers_with_reference(df_samples, most_freq_seq, kmer_col):
    """
    Returns a copy of df_samples where the specified kmer column is masked:
    - Any character in the kmer column that matches most_freq_seq is replaced with '-'
    - Rows where the kmer column is equal to most_freq_seq are removed
    - Drops the 'Sequence' column if it exists
    """
    df_masked = df_samples.copy()
    # Drop the 'Sequence' column if it exists
    if 'Sequence' in df_masked.columns:
        df_masked = df_masked.drop(columns=['Sequence'])

    # Remove rows where kmer_col == most_freq_seq
    mask = df_masked[kmer_col] != most_freq_seq
    df_masked = df_masked[mask].copy()

    def mask_seq(seq, ref):
        return ''.join('-' if a == b else a for a, b in zip(seq, ref)) if pd.notnull(seq) and pd.notnull(ref) else seq
    df_masked[kmer_col] = df_masked[kmer_col].apply(lambda x: mask_seq(x, most_freq_seq))

    return df_masked

