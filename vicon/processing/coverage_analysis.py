from ..utils.helpers import find_min_coverage_threshold
from ..io.fasta import read_fasta_to_dataframe

import itertools
import numpy as np
import pandas as pd


def abundant_kmers(df):
    """
    Finds abundant kmers and the samples covering them.
    """
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
    """
    Crops the DataFrame to the specified gene region and coverage threshold.
    """
    ldf = limit_to_l_gene(df, start, end)
    min_coverage = find_min_coverage_threshold(ldf, coverage_ratio)
    ldf = ldf.loc[:, ldf.sum() > min_coverage]
    print(f"DataFrame cropped to {ldf.shape[1]} columns with coverage above threshold.")
    return ldf

def limit_to_l_gene(df, start, end):
    """
    Limits the DataFrame to the specified gene region.
    """
    print(f"Limiting DataFrame to gene region from position {start} to {end}")
    ldf = df.loc[:, (df.columns >= start) & (df.columns <= end)]
    return ldf

def build_coverage_table(ldf):
    """
    Builds a coverage table for combinations of kmers.
    """
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
    """
    Returns the top kmers based on maximum union_sum and total coverage.
    """
    max_union_sum = cov_df['union_sum'].max()
    sub_df = cov_df[cov_df['union_sum'] == max_union_sum]
    max_total_cov = sub_df['k1_cov_plus_k2_cov'].max()
    sub_df = sub_df[sub_df['k1_cov_plus_k2_cov'] == max_total_cov]
    print(f"Top kmers found covering {max_union_sum} samples with total coverage {max_total_cov}")
    return sub_df


from collections import Counter

def find_most_frequent_and_calculate_mismatches(sequences):
    """
    Finds the most frequent sequence and calculates the total sum of mismatches
    of all other sequences from the most frequent sequence.

    Parameters:
        sequences (list of str): A list of DNA sequences.

    Returns:
        tuple: A tuple containing:
            - most_frequent (str): The most frequent sequence.
            - total_mismatches (int): The total sum of mismatches of other sequences from the most frequent one.
    """    
    
    # Find the most frequent sequence
    sequence_counts = Counter(sequences)
    most_frequent = sequence_counts.most_common(1)[0][0]
    
    # Calculate mismatches
    total_mismatches = 0
    for seq in sequences:
        if len(seq) != len(most_frequent):
            raise ValueError("All sequences must have the same length.")
        mismatches = sum(1 for a, b in zip(seq, most_frequent) if a != b)
        total_mismatches += mismatches

    # print(f"total_mismatches: {total_mismatches}")

    return most_frequent, total_mismatches


def get_i_th_kmers(fasta_file, i, mask, window_size=150):
    df = read_fasta_to_dataframe(fasta_file)
    df = df.iloc[mask[:, i].astype(bool)] 
    df['kmer'] = df['Sequence'].str.slice(i, i + window_size)

    return  df['kmer'].values, df['ID'].values


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


def count_sequences_with_max_mismatches(sequences, ids, most_frequent, max_mismatches=3):
    """
    Counts the number of sequences that have at most a specified number of mismatches
    compared to the most frequent sequence.

    Parameters:
        sequences (list of str): A list of DNA sequences.
        most_frequent (str): The most frequent sequence in the list.
        max_mismatches (int): The maximum number of allowed mismatches.

    Returns:
        int: The count of sequences with mismatches <= max_mismatches.
    """
    count = 0
    seq_indices = []
    for seq, i in zip(sequences, ids):
        # Count mismatches using zip
        mismatches = sum(1 for a, b in zip(seq, most_frequent) if a != b)
        if mismatches <= max_mismatches:
            count += 1
            seq_indices.append(i)

    return count, seq_indices


def count_seq_coverage(kmer_index, fasta_file, mask, window_size=150):
    seqs, ids = get_i_th_kmers(fasta_file, kmer_index, mask, window_size=window_size)
    most_freq, min_value = find_most_frequent_and_calculate_mismatches(seqs)
    coverage, seq_indices = count_sequences_with_max_mismatches(seqs, ids, most_freq, max_mismatches=3)
    return coverage, seq_indices, min_value


def find_best_pair_kmer(ldf, fasta_file, mask, window_size=150):
    """
    Finds the best pair of kmers using the most frequent sequences in each position
      based on coverage and mismatches.
    """

    # Find abundant kmers in ldf after the coverage threshold
    kmer_dict = dict()
    for c in ldf.columns:
        kmer_dict[c] = dict()
        coverage, seq_indices, min_value = count_seq_coverage(c, fasta_file, mask, window_size=window_size)
        kmer_dict[c]['coverage'] = coverage
        kmer_dict[c]['indices'] = seq_indices
        kmer_dict[c]['mismatches'] = min_value


    # Build the coverage table
    cov = np.zeros((len(ldf.columns), len(ldf.columns)))
    for i, c in enumerate(ldf.columns):
        for j, c2 in enumerate(ldf.columns):
            if i <= j:
                cov[i,j] = list(set(kmer_dict[c]['indices']).union(kmer_dict[c2]['indices'])).__len__()

    # Find the pair with the maximum coverage
    pairs = np.argwhere(cov == cov.max())
    print(f"Max coverage with 2 kmers: {cov.max()/ldf.shape[0]}")

    # Create the data for the DataFrame of the best kmers
    data = []
    for pair in pairs:
        kmer1 = ldf.columns[pair[0]]
        kmer2 = ldf.columns[pair[1]]

        # Extract values from kmer_dict
        coverage_k1 = kmer_dict[kmer1]["coverage"]
        coverage_k2 = kmer_dict[kmer2]["coverage"]
        mismatches_k1 = kmer_dict[kmer1]["mismatches"]
        mismatches_k2 = kmer_dict[kmer2]["mismatches"]

        # Append row to data
        data.append([kmer1, kmer2, coverage_k1, coverage_k2, mismatches_k1, mismatches_k2])

    # Create the DataFrame
    df_best_kmers = pd.DataFrame(data, columns=["kmer1", "kmer2", "coverage_k1", "coverage_k2", "mismatches_k1", "mismatches_k2"])
    df_best_kmers['sum_coverage'] = df_best_kmers['coverage_k1'] + df_best_kmers['coverage_k2']
    df_best_kmers['sum_mismatches'] = df_best_kmers['mismatches_k1'] + df_best_kmers['mismatches_k2']

    df_best_kmers = df_best_kmers.sort_values(by=["sum_coverage", "sum_mismatches"], ascending=[False, True])

    return df_best_kmers.iloc[0]['kmer1'], df_best_kmers.iloc[0]['kmer2']
