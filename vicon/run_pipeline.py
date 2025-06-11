## How to run??
#### python run_pipeline.py --config configs/config_rsva.yaml


import os
import shutil
import pandas as pd
import argparse

from vicon.dereplication.derep import run_vsearch
from vicon.alignment.ref_align import run_viralmsa
from vicon.processing.sample_processing import process_all_samples, pipeline_results_cleaner
from vicon.visualization.plots import plot_non_gap_counts, plot_rel_cons
from vicon.processing.coverage_analysis import crop_df, find_best_pair_kmer
from vicon.io.fasta import read_fasta_to_dataframe, remove_first_record
from vicon.utils.helpers import (
    load_config,
    count_non_gap_characters_from_dataframe,
    filter_by_most_common_kmers
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
    # base_path = os.path.join(config["project_path"], "vicon")
    base_path = config["project_path"]
    virus = config["virus_name"]

    input_sample = os.path.join(base_path, config["input_sample"])
    input_reference = os.path.join(base_path, config["input_reference"])
    output_dir = os.path.join(base_path, "results", virus)

    sample_dir = os.path.dirname(input_sample)
    aligned_dir = os.path.join(sample_dir, "aligned")
    derep_fasta = os.path.join(sample_dir, "derep.fasta")
    clusters_uc = os.path.join(sample_dir, "clusters.uc")
    derep_fasta_aln = os.path.join(aligned_dir, "derep.fasta.aln")
    kmer1_path = os.path.join(output_dir, "kmer1.csv")
    kmer2_path = os.path.join(output_dir, "kmer2.csv")
    log_dir = os.path.join(output_dir, "logs")

    email = config["email"]
    kmer_size = config["kmer_size"]
    threshold = config["threshold"]
    l_gene_start = config["l_gene_start"]
    l_gene_end = config["l_gene_end"]
    coverage_ratio = config["coverage_ratio"]

    # Setup
    print(f"[INFO] Using base path: {base_path}")
    print(f"aligned_dir: {aligned_dir}")
    setup_alignment_directory(aligned_dir)

    # Dereplication and Alignment
    run_vsearch(input_sample, derep_fasta, clusters_uc)

    run_viralmsa(
        email=email,
        sample_fasta=derep_fasta,
        output_dir=aligned_dir,
        reference_fasta=input_reference,
    )

    remove_first_record(derep_fasta_aln, derep_fasta_aln)

    # Process Samples
    df3, mask3 = process_all_samples(
        input_reference, derep_fasta_aln, log_dir,
        window_size=kmer_size, threshold=threshold, only_valid_kmers=True
    )
    df3.columns = df3.columns.astype(int)

    plot_rel_cons(df3, kmer_size=kmer_size, threshold=kmer_size-threshold, save_path=output_dir, sample_name=virus)

    # L-gene region crop and kmer detection
    ldf = crop_df(df3, l_gene_start, l_gene_end, coverage_ratio=coverage_ratio)
    kmer1, kmer2 = find_best_pair_kmer(
        ldf, derep_fasta_aln, mask3,
        sort_by_mismatches=False, window_size=kmer_size
    )

    # print(f"[INFO] Kmer1: {kmer1}, Kmer2: {kmer2}")

    kmer1_seq, kmer2_seq = extract_kmer_sequences(input_reference, kmer1, kmer2, kmer_size)
    print(f"[INFO] Degenerate Kmer1 sequence (from reference) (position {kmer1}):\n{kmer1_seq}")
    print(f"[INFO] Degenerate Kmer2 sequence (from reference) (position {kmer2}):\n{kmer2_seq}")

    # Clean results
    df_kmers1, df_kmers2, df_samples = pipeline_results_cleaner(
        sample_address=derep_fasta_aln,
        kmer1=kmer1,
        kmer2=kmer2,
        drop_old_samples=config["drop_old_samples"],
        kmer_size=kmer_size,
        min_year=config["min_year"],
        threshold_ratio=config["threshold_ratio"],
        drop_mischar_samples=config["drop_mischar_samples"]
    )

    df_kmers1.to_csv(kmer1_path)
    df_kmers2.to_csv(kmer2_path)

    # Non-gap mutation plots
    for df_kmers, label in [(df_kmers1, "kmers1"), (df_kmers2, "kmers2")]:
        df_counts = count_non_gap_characters_from_dataframe(df_kmers, sequence_column='alignment') - 1
        plot_non_gap_counts(
            df_counts,
            title=f'{virus} – {label} Non-Gap Counts',
            save=os.path.join(output_dir, f"{label}_mutations.png")
        )

    # Filter by common kmers
    filtered_df, kmer1_most, kmer2_most, kmer1_count, kmer2_count = filter_by_most_common_kmers(df_samples)
    print(f"[INFO] Native kmer1 sequence: {kmer1_most}")
    print(f"[INFO] Native Kmer1 count: {kmer1_count}")
    print(f"[INFO] Native kmer2 sequence: {kmer2_most}")
    print(f"[INFO] Native Kmer2 count: {kmer2_count}")
    print(f"[INFO] Overall Native Coverage : {filtered_df.shape[0]} out of {df_samples.shape[0]}")


if __name__ == "__main__":
    main()
