import os
import numpy as np
import pandas as pd
from .sliding_window import sliding_window_matches, calculate_validity_mask
from .kmer_analysis import process_kmers_to_dataframe, filter_kmers_by_changes
from ..io.fasta import read_fasta, read_fasta_to_dataframe
from vicon.utils.year_extraction import extract_year, filter_group_by_year
from vicon.processing.kmer_analysis import process_kmers_to_dataframe

def process_all_samples(reference_path,
                        samples_fasta,
                        output_dir,
                        window_size=150,
                        threshold=145,
                        only_valid_kmers=True,
                        ):
    """
    Processes all samples against the reference and aggregates binary results.

    Parameters:
        reference_path (str): Path to the reference FASTA file.
        samples_dir (str): Directory containing sample directories with genome.fasta files.
        output_dir (str): Directory to save individual binary results.
        window_size (int): Size of the sliding window.
        threshold (int): Threshold for binary match results.
        only_valid_kmers (bool): If True, applies validity mask to the results.

    Returns:
        pd.DataFrame: Aggregated binary results as a DataFrame.
    """
    # Read the reference sequence
    reference_seq = read_fasta(reference_path)

    sample_seqs = read_fasta_to_dataframe(samples_fasta)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Prepare to aggregate results
    aggregated_data = []
    sample_ids = []
    mask = []

    for index, row in sample_seqs.iterrows():
        sample_seq = row['Sequence']
        binary_result = sliding_window_matches(reference_seq, sample_seq, window_size, threshold)
        if only_valid_kmers:
            valid_mask = calculate_validity_mask(reference_seq, sample_seq, window_size)
            mask.append(valid_mask)
            binary_result = binary_result & valid_mask

        # Append binary result to aggregation
        aggregated_data.append(binary_result)
        sample_ids.append(row['ID'])

    # Aggregate all results into a DataFrame
    aggregated_array = np.vstack(aggregated_data)
    df = pd.DataFrame(aggregated_array, index=sample_ids)
    df.index.name = "SampleID"

    if only_valid_kmers:
        mask_array = np.vstack(mask)
    else:
        mask_array = None

    return df, mask_array


def pipeline_results_cleaner(sample_address='../data/rsva/samples/aligned/derep.fasta.aln',
                             kmer1=11000,
                             kmer2=12000,
                             drop_old_samples=True,
                             drop_mutants=False,
                             kmer_size= 150,
                             min_year=2020,
                             threshold_ratio=0.01,
                             drop_mischar_samples=False,
                             return_droped_samples=False,
                             logger=None
                             ):
    
    df_samples = read_fasta_to_dataframe(sample_address)
    if logger:
        logger.info(f"Read {df_samples.shape[0]} samples from {sample_address}.")
    else:
        print(f"Read {df_samples.shape[0]} samples from {sample_address}.")

    df_samples['kmer1'] = df_samples['Sequence'].str.slice(kmer1, kmer1+kmer_size)
    df_samples['kmer2'] = df_samples['Sequence'].str.slice(kmer2, kmer2+kmer_size)

    df_samples['year'] = df_samples['ID'].apply(extract_year)

    # Convert the 'year' column to numeric, coercing errors to NaN
    df_samples['year'] = pd.to_numeric(df_samples['year'], errors='coerce')

    # Drop rows where 'year' is NaN (i.e., non-numerical entries were converted to NaN)
    if logger:
        logger.info(f"df_samples shape before dropping NaN years: {df_samples.shape}")
    else:
        print(f"df_samples shape before dropping NaN years: {df_samples.shape}")
    df_samples_without_year = df_samples[df_samples['year'].isna()]
    df_samples = df_samples.dropna(subset=['year'])
    if logger:
        logger.info(f"df_samples shape after dropping NaN years: {df_samples.shape}")
    else:
        print(f"df_samples shape after dropping NaN years: {df_samples.shape}")

    ## Drop the samples containing Non ATCG chars in their kmer1 or kmer2
    # print(df_samples['kmer2'])
    if drop_mischar_samples:
        df_sample_s_with_mischar = df_samples[df_samples['kmer1'].str.contains('[^ATCG]', regex=True) |
                                              df_samples['kmer2'].str.contains('[^ATCG]', regex=True)]
        df_samples = df_samples[~df_samples['kmer1'].str.contains('[^ATCG]', regex=True) & 
                                ~df_samples['kmer2'].str.contains('[^ATCG]', regex=True)]
        if logger:
            logger.info(f"df_samples shape after dropping samples with non ATCG chars: {df_samples.shape}")
        else:
            print(f"df_samples shape after dropping samples with non ATCG chars: {df_samples.shape}")


    # drop old samples 
    if drop_old_samples:
        # Calculate the threshold for 1% of the total number of rows
        threshold = len(df_samples) * threshold_ratio

        # Group by 'kmer1'
        gk1 = df_samples.groupby('kmer1')
        gk2 = df_samples.groupby('kmer2')

        # Apply the filter and create a new DataFrame with only the groups to be kept
        filtered_df1 = gk1.filter(filter_group_by_year, min_year=min_year, threshold=threshold)
        filtered_df2 = gk2.filter(filter_group_by_year, min_year=min_year, threshold=threshold)

        filtered_ids1 = filtered_df1['ID']
        filtered_ids2 = filtered_df2['ID']

        ids = list(set(filtered_ids1.values).intersection(filtered_ids2.values))
        df_samples = df_samples[df_samples['ID'].isin(ids)]

    # Drop kmers which have more than 3 changes to the most frequent kmer
    if drop_mutants==True:
        df_changes1 = filter_kmers_by_changes(df_samples, kmer_name='kmer1', num_changes_threshold=3)
        df_changes2 = filter_kmers_by_changes(df_samples, kmer_name='kmer2', num_changes_threshold=3)

        ids_to_drop = list(set(df_changes1.ID.tolist()).union(df_changes2.ID.tolist()))
        df_samples = df_samples.loc[~df_samples['ID'].isin(ids_to_drop)]

    df_kmers1 = process_kmers_to_dataframe(df_samples, 'kmer1')
    df_kmers1 = df_kmers1.sort_values(by=["e_year", "s_year"], ascending=[False, False])

    df_kmers2 = process_kmers_to_dataframe(df_samples, 'kmer2')
    df_kmers2 = df_kmers2.sort_values(by=["e_year", "s_year"], ascending=[False, False])

    if return_droped_samples:
        return df_kmers1, df_kmers2, df_samples_without_year, df_sample_s_with_mischar
    else:
        return df_kmers1, df_kmers2, df_samples


