from ..utils.helpers import find_min_coverage_threshold
from ..io.fasta import read_fasta_to_dataframe
import itertools
import numpy as np
import pandas as pd
import time
from multiprocessing import Pool, cpu_count
from collections import Counter

def abundant_kmers(df):
    """Finds abundant kmers and the samples covering them."""
    print("Finding abundant kmers...")
    data = df.copy()
    output = dict()
    samples = dict()
    while True:
        pos = data.sum().idxmax()
        rep = data.sum().max()
        output[pos] = int(rep)
        print(f"Kmer start pos: {pos} has appeared in {int(rep)} samples / {data.shape[0]} samples")
        samples[pos] = data[data[pos]==1].index.values
        data = data[data[pos]!=1]
        total_remaining = data.sum().sum()
        if total_remaining <1:
            print("No more kmers left to process.")
            break
        else:
            print(f"Remaining kmers to process: {total_remaining}")
    return output, samples

def crop_df(df, start, end, coverage_ratio=0.5):
    """Crops the DataFrame to the specified gene region and coverage threshold."""
    ldf = limit_to_l_gene(df, start, end)
    min_coverage = find_min_coverage_threshold(ldf, coverage_ratio)
    ldf = ldf.loc[:, ldf.sum() > min_coverage]
    print(f"DataFrame cropped to {ldf.shape[1]} columns with coverage above threshold.")
    return ldf

def limit_to_l_gene(df, start, end):
    """Limits the DataFrame to the specified gene region."""
    print(f"Limiting DataFrame to gene region from position {start} to {end}")
    ldf = df.loc[:, (df.columns >= start) & (df.columns <= end)]
    return ldf

def build_coverage_table(ldf):
    """Builds a coverage table for combinations of kmers."""
    print("Building coverage table for kmers...")
    data = ldf.values
    column_names = ldf.columns.tolist()
    num_kmers = len(column_names)
    coverage_list = []

    for i, j in itertools.combinations(range(num_kmers), 2):
        k1, k2 = column_names[i], column_names[j]
        union_sum = np.sum((data[:, i] + data[:, j]) > 0)
        k1_cov = data[:, i].sum()
        k2_cov = data[:, j].sum()
        total_cov = k1_cov + k2_cov
        coverage_list.append({
            'k1': k1,
            'k2': k2,
            'union_sum': union_sum,
            'k1_cov': k1_cov,
            'k2_cov': k2_cov,
            'k1_cov_plus_k2_cov': total_cov
        })
    coverage_df = pd.DataFrame(coverage_list)
    coverage_df = coverage_df.loc[(coverage_df[['k1_cov', 'k2_cov', 'union_sum']] != 0).any(axis=1)]
    coverage_df = coverage_df.sort_values(by=['union_sum', 'k1_cov_plus_k2_cov'], ascending=False)
    print(f"Coverage table built with {len(coverage_df)} combinations.")
    return coverage_df

def top_kmers_df(cov_df):
    """Returns the top kmers based on maximum union_sum and total coverage."""
    max_union_sum = cov_df['union_sum'].max()
    sub_df = cov_df[cov_df['union_sum'] == max_union_sum]
    max_total_cov = sub_df['k1_cov_plus_k2_cov'].max()
    sub_df = sub_df[sub_df['k1_cov_plus_k2_cov'] == max_total_cov]
    print(f"Top kmers found covering {max_union_sum} samples with total coverage {max_total_cov}")
    return sub_df

def find_most_frequent_and_calculate_mismatches(sequences):
    """Finds the most frequent sequence and calculates mismatches."""
    sequence_counts = Counter(sequences)
    most_frequent = sequence_counts.most_common(1)[0][0]
    total_mismatches = sum(1 for seq in sequences 
                          for a, b in zip(seq, most_frequent) if a != b)
    return most_frequent, total_mismatches

def get_i_th_kmers(fasta_file, i, mask, window_size=150):
    df = read_fasta_to_dataframe(fasta_file)
    df = df.iloc[mask[:, i].astype(bool)] 
    df['kmer'] = df['Sequence'].str.slice(i, i + window_size)
    return df['kmer'].values, df['ID'].values

def count_sequences_with_max_mismatches(sequences, ids, most_frequent, max_mismatches=3):
    """Counts sequences with <= max_mismatches compared to most_frequent."""
    count = 0
    seq_indices = []
    for seq, i in zip(sequences, ids):
        mismatches = sum(1 for a, b in zip(seq, most_frequent) if a != b)
        if mismatches <= max_mismatches:
            count += 1
            seq_indices.append(i)
    return count, seq_indices

def count_seq_coverage(kmer_index, fasta_file, mask, window_size=150):
    seqs, ids = get_i_th_kmers(fasta_file, kmer_index, mask, window_size)
    most_freq, min_value = find_most_frequent_and_calculate_mismatches(seqs)
    coverage, seq_indices = count_sequences_with_max_mismatches(seqs, ids, most_freq, 3)
    return coverage, seq_indices, min_value

def process_kmer_column(args):
    """Process single kmer column in parallel."""
    c, fasta_file, mask, window_size = args
    coverage, seq_indices, min_value = count_seq_coverage(c, fasta_file, mask, window_size)
    return (c, {'coverage': coverage, 'indices': seq_indices, 'mismatches': min_value})

def process_coverage_chunk(args):
    """Process chunk of coverage matrix in parallel."""
    chunk_indices, columns, kmer_indices = args
    results = []
    for i, j in chunk_indices:
        c1 = columns[i]
        c2 = columns[j]
        results.append((i, j, len(set(kmer_indices[c1]).union(kmer_indices[c2]))))
    return results

def find_best_pair_kmer(ldf, fasta_file, mask, window_size=150, sort_by_mismatches=True ,n_processes=None):
    """
    Finds the best pair of kmers from a DataFrame using parallel processing.
    This function evaluates all possible pairs of kmer columns in the input DataFrame (`ldf`)
    to identify the pair that maximizes coverage (and optionally minimizes mismatches).
    It leverages multiprocessing to speed up both the per-kmer and per-pair computations.
    Args:
        ldf (pd.DataFrame): DataFrame where each column represents a kmer.
        fasta_file (str): Path to the FASTA file used for coverage analysis.
        mask (np.ndarray or similar): Mask to apply during kmer processing.
        window_size (int, optional): Window size for kmer analysis. Defaults to 150.
        sort_by_mismatches (bool, optional): If True, sorts best pairs by coverage and then mismatches. Defaults to True.
        n_processes (int or None, optional): Number of processes to use for parallelization. Defaults to all available CPUs.
    Returns:
        tuple: The names of the two best kmer columns (kmer1, kmer2) as strings.
    Notes:
        - Requires the helper functions `process_kmer_column` and `process_coverage_chunk` to be defined.
        - The function prints timing and progress information to stdout.
    """

    start_time = time.time()
    
    if n_processes is None:
        n_processes = cpu_count()
        print(f"Using {n_processes} processes")

    # Parallel process kmer_dict
    with Pool(n_processes) as pool:
        args = [(c, fasta_file, mask, window_size) for c in ldf.columns]
        results = pool.map(process_kmer_column, args)
    
    kmer_dict = dict(results)
    print(f"Kmer dict computed in {time.time()-start_time:.2f}s")
    
    # Extract just the indices for parallel processing
    kmer_indices = {k: v['indices'] for k, v in kmer_dict.items()}
    columns = ldf.columns
    n = len(columns)
    cov = np.zeros((n, n))
    
    # Split work into chunks
    indices = [(i, j) for i in range(n) for j in range(i, n)]
    chunk_size = max(1, len(indices) // (n_processes * 4))  # 4 chunks per core
    chunks = [indices[i:i+chunk_size] for i in range(0, len(indices), chunk_size)]
    
    # Process chunks in parallel
    with Pool(n_processes) as pool:
        results = pool.map(process_coverage_chunk, 
                         [(chunk, columns, kmer_indices) for chunk in chunks])
    
    # Combine results
    for chunk_results in results:
        for i, j, value in chunk_results:
            cov[i, j] = value
    
    print(f"Coverage matrix computed in {time.time()-start_time:.2f}s")
    
    # Get all pairs and their unique coverage
    pair_indices = [(i, j) for i in range(n) for j in range(i, n)]
    data = []
    for i, j in pair_indices:
        unique_cov = int(cov[i, j])
        data.append([
            columns[i], columns[j],
            kmer_dict[columns[i]]['coverage'],
            kmer_dict[columns[j]]['coverage'],
            kmer_dict[columns[i]]['mismatches'],
            kmer_dict[columns[j]]['mismatches'],
            unique_cov
        ])

    df_best = pd.DataFrame(data, columns=["kmer1", "kmer2", "cov1", "cov2", "mism1", "mism2", "unique_cov"])
    df_best['sum_cov'] = df_best['cov1'] + df_best['cov2']
    df_best['sum_mism'] = df_best['mism1'] + df_best['mism2']

    # Sort by unique coverage first, then sum_cov, then mismatches
    df_best = df_best.sort_values(['unique_cov', 'sum_cov', 'sum_mism'], ascending=[False, False, True])

    # Take the top 1000 pairs
    df_best = df_best.head(1000)

    print(f"[INFO] Degenerate Kmer1 Coverage: {df_best.iloc[0]['cov1']}")
    print(f"[INFO] Degenerate Kmer2 Coverage: {df_best.iloc[0]['cov2']}")
    print(f"[INFO] Overall Degenerate Coverage: {df_best.iloc[0]['sum_cov']}")
    
    return df_best.iloc[0]['kmer1'], df_best.iloc[0]['kmer2']


def select_best_kmers(fasta_file, mask, kmer_set, set2, window_size=150):

    def find_min_mismatch_and_max_coverage(kmer_set, fasta_file, mask, window_size=150):
        min_kmer_set = 1e6
        best_kmer_set = None
        max_coverage = 0

        for i in kmer_set:
            seqs, ids = get_i_th_kmers(fasta_file, i ,mask, window_size=window_size)
            most_freq, min_value = find_most_frequent_and_calculate_mismatches(seqs)
            coverage, seq_indices = count_sequences_with_max_mismatches(seqs, most_freq, max_mismatches=3)
            if coverage > max_coverage:
                max_coverage = coverage
                min_kmer_set = min_value
                best_kmer_set = i
            elif coverage == max_coverage:
                if min_value < min_kmer_set:
                    min_kmer_set = min_value
                    best_kmer_set = i
        return best_kmer_set, min_kmer_set, max_coverage
    
    best_kmer1, min_kmer1, max_coverage1 = find_min_mismatch_and_max_coverage(kmer_set, fasta_file, mask, window_size=window_size)
    print(f"Best kmer1: {best_kmer1}, min_kmer1: {min_kmer1}, max_coverage1: {max_coverage1}")
    best_kmer2, min_kmer2, max_coverage2 = find_min_mismatch_and_max_coverage(set2, fasta_file, mask, window_size=window_size)
    print(f"Best kmer2: {best_kmer2}, min_kmer2: {min_kmer2}, max_coverage2: {max_coverage2}")
    return best_kmer1, best_kmer2


def calculate_kmer_coverage(ldf, fasta_file, mask, kmer1, kmer2, window_size=150):
    """
    Calculate coverage for a pair of kmers using the same method as find_best_pair_kmer.
    
    Parameters:
    -----------
    ldf : pandas.DataFrame
        DataFrame containing sequence data
    fasta_file : str
        Path to the FASTA file
    mask : numpy.ndarray
        Mask array for filtering sequences
    kmer1 : int
        Position of first kmer
    kmer2 : int
        Position of second kmer
    window_size : int, optional
        Size of the kmers to analyze (default: 150)
    
    Returns:
    --------
    int
        Number of sequences that have both kmers at their respective positions
    """
    # Get coverage and indices for each kmer using the same method as find_best_pair_kmer
    coverage1, indices1, _ = count_seq_coverage(kmer1, fasta_file, mask, window_size=window_size)
    coverage2, indices2, _ = count_seq_coverage(kmer2, fasta_file, mask, window_size=window_size)
    
    # Calculate union of indices (same as in find_best_pair_kmer)
    union_indices = set(indices1).union(set(indices2))
    
    return len(union_indices)

def find_kmer_position(df_ref, kmer_sequence, window_size=150):
    """
    Find the position of a kmer sequence in the reference sequence.
    
    Parameters:
    -----------
    df_ref : pandas.DataFrame
        DataFrame containing the reference sequence
    kmer_sequence : str
        The kmer sequence to find
    window_size : int, optional
        Size of the kmers (default: 150)
    
    Returns:
    --------
    int or None
        Position of the kmer sequence in the reference sequence, or None if not found
    """
    # Read reference sequence

    ref_seq = df_ref['Sequence'].values[0]
    
    # Search for the kmer sequence
    for i in range(len(ref_seq) - window_size + 1):
        if ref_seq[i:i+window_size] == kmer_sequence:
            return i
    
    return None
    