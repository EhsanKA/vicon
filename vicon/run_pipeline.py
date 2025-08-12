## How to run??
#### python run_pipeline.py --config configs/config_rsva.yaml


import os
import shutil
import pandas as pd
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vicon.dereplication.derep import run_vsearch
from vicon.alignment.ref_align import run_viralmsa
from vicon.processing.sample_processing import process_all_samples, pipeline_results_cleaner
from vicon.visualization.plots import plot_non_gap_counts, plot_rel_cons
from vicon.processing.coverage_analysis import crop_df, find_best_pair_kmer
from vicon.io.fasta import read_fasta_to_dataframe, remove_first_record, generate_remaining_fasta
from vicon.processing.kmer_analysis import mask_kmers_with_reference
from vicon.utils.helpers import (
    load_config,
    count_non_gap_characters_from_dataframe,
    filter_by_most_common_kmers,
    process_fasta_file
)


def setup_alignment_directory(aligned_dir):
    if os.path.exists(aligned_dir):
        shutil.rmtree(aligned_dir)
    # os.makedirs(aligned_dir)


def extract_kmer_sequences(reference_path, kmer1, kmer2, kmer_size):
    df_ref = read_fasta_to_dataframe(reference_path)
    ref_seq = df_ref['Sequence'].values[0]
    return ref_seq[kmer1:kmer1+kmer_size], ref_seq[kmer2:kmer2+kmer_size]


def parse_args():
    parser = argparse.ArgumentParser(description="Run the VICON pipeline.")
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to the YAML config file (default: config.yaml)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)

    # Paths and Constants
    base_path = config.get("base_path", "")  
    virus = config["virus_name"]

    input_sample = os.path.join(base_path, config["input_sample"])
    input_reference = os.path.join(base_path, config["input_reference"])
    output_dir = os.path.join(base_path, "results", virus)
    # output_dir = os.path.join(base_path, "results/reproduce_old_results", virus)

    # --- FASTA cleaning step: process and replace input files with cleaned versions ---
    input_sample_upper = input_sample.replace('.fasta', '_upper.fasta')
    input_reference_upper = input_reference.replace('.fasta', '_upper.fasta')
    process_fasta_file(input_sample, input_sample_upper)
    process_fasta_file(input_reference, input_reference_upper)
    input_sample = input_sample_upper
    input_reference = input_reference_upper
    # -------------------------------------------------------------------------------

    sample_dir = os.path.dirname(input_sample)
    aligned_dir = os.path.join(sample_dir, "aligned")
    derep_fasta = os.path.join(sample_dir, "derep.fasta")
    clusters_uc = os.path.join(sample_dir, "clusters.uc")
    derep_fasta_aln = os.path.join(aligned_dir, "derep.fasta.aln")
    kmer1_path = os.path.join(output_dir, "kmer1.csv")
    kmer2_path = os.path.join(output_dir, "kmer2.csv")
    log_dir = os.path.join(output_dir, "logs")

    email = config.get("email", "example@example.com")
    kmer_size = config.get("kmer_size", 150)
    threshold = config.get("threshold", 147)
    only_valid_kmers = config.get("only_valid_kmers", True)
    l_gene_start = config.get("l_gene_start", -1)
    l_gene_end = config.get("l_gene_end", 40000)
    coverage_ratio = config.get("coverage_ratio", 0.5)
    sort_by_mismatches = config.get("sort_by_mismatches", True)
    
    drop_old_samples = config.get("drop_old_samples", False)
    min_year = config.get("min_year", 2020)
    threshold_ratio = config.get("threshold_ratio", 0.1)
    drop_mischar_samples = config.get("drop_mischar_samples", True)

    # Define a logger
    logger = logging.getLogger('pipeline_logger')
    logger.setLevel(logging.INFO)

    # Create a file handler with dynamic log file name
    os.makedirs(log_dir, exist_ok=True)  # <-- This must come BEFORE FileHandler
    log_file_path = os.path.join(log_dir, 'pipeline.log')
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)

    # Create a logging format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    # Add the handler to the logger
    logger.addHandler(file_handler)

    # Ensure the log file is created fresh for each run
    if os.path.exists(log_file_path):
        os.remove(log_file_path)
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setFormatter(formatter)  # Reapply formatter
    logger.addHandler(file_handler)  # Re-add handler to logger

    # Log configuration parameters
    logger.info("Pipeline Configuration Parameters:")
    for key, value in config.items():
        logger.info(f"{key}: {value}")


    # Pass logger to coverage_analysis
    from vicon.processing.coverage_analysis import crop_df, find_best_pair_kmer
    coverage_analysis_logger = logger

    # Setup
    logger.info(f"[INFO] Using base path: {base_path}")
    logger.info(f"aligned_dir: {aligned_dir}")
    setup_alignment_directory(aligned_dir)

    # Dereplication and Alignment
    run_vsearch(input_sample, derep_fasta, clusters_uc, logger=logger)

    run_viralmsa(
        email=email,
        sample_fasta=derep_fasta,
        output_dir=aligned_dir,
        reference_fasta=input_reference,
    )

    remove_first_record(derep_fasta_aln, derep_fasta_aln, logger=logger)

    # Process Samples
    df3, mask3 = process_all_samples(
        input_reference, derep_fasta_aln, log_dir,
        window_size=kmer_size, threshold=threshold, only_valid_kmers=only_valid_kmers
    )
    df3.columns = df3.columns.astype(int)

    print(f"max kmer present in the samples at the position: {df3.sum(axis=0).idxmax()}, max count is: {df3.sum(axis=0).max()}")

    plot_rel_cons(df3, kmer_size=kmer_size, threshold=kmer_size-threshold, save_path=output_dir, sample_name=virus)

    # L-gene region crop and kmer detection
    ldf = crop_df(df3, l_gene_start, l_gene_end, coverage_ratio=coverage_ratio, logger=logger)
    kmer1, kmer2 = find_best_pair_kmer(
        ldf, derep_fasta_aln, mask3,
        sort_by_mismatches=sort_by_mismatches, window_size=kmer_size, logger=logger
    )
    # kmer1, kmer2 = int(config['kmer1']), int(config['kmer2'])

    kmer1_seq, kmer2_seq = extract_kmer_sequences(input_reference, kmer1, kmer2, kmer_size)
    logger.info(f"[INFO] Degenerate Kmer1 sequence (from reference) (position {kmer1})")
    logger.info(f"[INFO] Degenerate Kmer2 sequence (from reference) (position {kmer2})")

    # Clean results
    df_kmers1, df_kmers2, df_samples = pipeline_results_cleaner(
        sample_address=derep_fasta_aln,
        kmer1=kmer1,
        kmer2=kmer2,
        drop_old_samples=drop_old_samples,
        kmer_size=kmer_size,
        min_year=min_year,
        threshold_ratio=threshold_ratio,
        drop_mischar_samples=drop_mischar_samples,
        logger=logger
    )

    df_kmers1.to_csv(kmer1_path)
    df_kmers2.to_csv(kmer2_path)

    # Non-gap mutation plots
    for df_kmers, label in [(df_kmers1, "kmers1"), (df_kmers2, "kmers2")]:
        df_counts = count_non_gap_characters_from_dataframe(df_kmers, sequence_column='alignment') - 1
        plot_non_gap_counts(
            df_counts,
            title=f'{virus} - {label} Non-Gap Counts',
            save=os.path.join(output_dir, f"{label}_mutations.png"),
            logger=logger
        )

    # Filter by common kmers
    filtered_df, kmer1_most, kmer2_most, kmer1_count, kmer2_count = filter_by_most_common_kmers(df_samples)
    logger.info("*"*100)
    logger.info("output summary: [cleaned samples , native kmers] ")
    logger.info("*"*100)
    logger.info(f"[INFO] Native kmer1 sequence: \n {kmer1_most}")
    logger.info(f"[INFO] Native kmer2 sequence: \n {kmer2_most}")
    logger.info(f"[INFO] Native Kmer1 count: {kmer1_count}")
    logger.info(f"[INFO] Native Kmer2 count: {kmer2_count}")
    logger.info(f"[INFO] Overall Native Coverage : {filtered_df.shape[0]} out of {df_samples.shape[0]}")

    # Mask and save kmer1
    sample_id_with_kmer1_masked = mask_kmers_with_reference(df_samples, kmer1_most, 'kmer1')
    masked_kmer1_path = os.path.join(output_dir, "sample_id_with_kmer1_masked.csv")
    sample_id_with_kmer1_masked = sample_id_with_kmer1_masked.drop(columns=['kmer2', 'year'], errors='ignore')  # Drop kmer2 column if it exists
    sample_id_with_kmer1_masked.to_csv(masked_kmer1_path, index=False)
    logger.info(f"Saved sample_id_with_kmer1_masked.csv to {masked_kmer1_path}")
    logger.info(f"This file contains the sample IDs along with their corresponding kmer1 sequences (showing differences from the native kmer1) and years.")

    # Mask and save kmer2
    sample_id_with_kmer2_masked = mask_kmers_with_reference(df_samples, kmer2_most, 'kmer2')
    masked_kmer2_path = os.path.join(output_dir, "sample_id_with_kmer2_masked.csv")
    sample_id_with_kmer2_masked = sample_id_with_kmer2_masked.drop(columns=['kmer1', 'year'], errors='ignore')  # Drop kmer1 column if it exists
    sample_id_with_kmer2_masked.to_csv(masked_kmer2_path, index=False)
    logger.info(f"Saved sample_id_with_kmer2_masked.csv to {masked_kmer2_path}")
    logger.info(f"This file contains the sample IDs along with their corresponding kmer2 sequences (showing differences from the native kmer2) and years.")

    # Calculate overall degenerate coverage using the binary matrix for cleaned samples
    if kmer1 in df3.columns and kmer2 in df3.columns:
        # Only keep rows (samples) present in both df3 and df_samples
        common_idx = df3.index.intersection(df_samples['ID'])
        deg_df = df3.loc[common_idx, [kmer1, kmer2]]
        # For each sample, check if either kmer is present
        deg_covered = (deg_df.sum(axis=1) > 0).sum()
        # Print kmer1 and kmer2 sequences and their counts from df3 among cleaned samples
        logger.info("*"*100)
        logger.info("output summary: [cleaned samples , degenerate kmers] ")
        logger.info("*"*100)
        logger.info(f"[INFO] Cleaned Degenerate Kmer1 sequence (from reference): \n {kmer1_seq}")
        logger.info(f"[INFO] Cleaned Degenerate Kmer2 sequence (from reference): \n {kmer2_seq}")
        logger.info(f"[INFO] Cleaned Degenerate Kmer1 count (from binary matrix): {deg_df[kmer1].sum()}")
        logger.info(f"[INFO] Cleaned Degenerate Kmer2 count (from binary matrix): {deg_df[kmer2].sum()}")
        logger.info(f"[INFO] Cleaned Overall Degenerate Coverage (from binary matrix): {deg_covered} out of {deg_df.shape[0]}")
    
    #     # Generate remaining sequences FASTA file
    #     remaining_fasta_path = os.path.join(output_dir, "remaining_sequences.fasta")
    #     remaining_count = generate_remaining_fasta(
    #         df3, df_samples, kmer1, kmer2, 
    #         derep_fasta_aln, remaining_fasta_path, logger
    #     )
    #     logger.info(f"[INFO] Generated remaining sequences file with {remaining_count} sequences: {remaining_fasta_path}")
    else:
        logger.warning("[WARN] kmer1 or kmer2 not found in df3 columns for degenerate coverage calculation.")

    print("*"*100)
    print("- vicon pipeline finished")
    print(f"- Check the {log_file_path} for more details.")
    print(f"- Logs are being acuumulated in {log_file_path}! to have a clean log, please remove the log file before running the pipeline again.")
    print("*"*100)


if __name__ == "__main__":
    main()
