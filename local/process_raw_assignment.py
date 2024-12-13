#!/usr/bin/env python3
"""
Author: Bhargav Pulugundla
Last Updated: 13-Dec-2024

This script references code from the GitHub repository, https://github.com/GCA-VH-lab/RiboSeqPy
"""

import os
import pandas as pd
import argparse
from pathlib import Path
import pysam
from collections import Counter, defaultdict

import warnings
import tables

warnings.filterwarnings("ignore", category=tables.NaturalNameWarning)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Script for assigning raw sequencing reads directly to genomic features according to the read lengths."
    )
    parser.add_argument(
        "--rawdir", type=str, help="Directory containing raw FASTQ files."
    )
    parser.add_argument("--names", type=str, help="Sample names separated by space.")
    parser.add_argument(
        "--mapped_twice",
        type=int,
        default=0,
        choices=[0, 1],
        help="Include reads mapped twice.",
    )
    parser.add_argument(
        "--mapping",
        type=str,
        default="5",
        choices=["5", "3"],
        help="Map reads based on 5' or 3' end.",
    )
    parser.add_argument(
        "--read_len_min", type=int, default=25, help="Minimum read length."
    )
    parser.add_argument(
        "--read_len_max", type=int, default=35, help="Maximum read length."
    )
    return parser.parse_args()


def yeast_chromosomes():
    """Return an ordered list of yeast chromosomes."""
    return [
        "I",
        "II",
        "III",
        "IV",
        "V",
        "VI",
        "VII",
        "VIII",
        "IX",
        "X",
        "XI",
        "XII",
        "XIII",
        "XIV",
        "XV",
        "XVI",
        "Mito",
    ]


def update_dataframe(df, chromosome, strand):
    """
    Update the DataFrame to include chromosome, position, and strand information.
    """
    df.fillna(0, inplace=True)
    df["sum"] = df.sum(axis=1)
    columns = list(df.columns)
    columns = ["Chromosome", "Position", "Strand"] + columns
    df["Chromosome"] = chromosome
    df["Strand"] = strand
    df["Position"] = df.index
    return df[columns]


def save_to_csv(dataframe, output_file):
    """Save DataFrame to a CSV file."""
    dataframe.to_csv(output_file, sep="\t", header=True, index=True)


def process_raw_assignment(args):
    """
    Perform raw assignment of reads for the given samples and save results to HDF5 and optionally CSV.
    """
    base_dir = Path(args.rawdir).resolve().parent
    assign_dir = base_dir / "5-AssignRaw"
    reports_dir = assign_dir / "Reports"
    os.makedirs(assign_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    save_csv = True  # Enable saving output to CSV

    read_len_min, read_len_max = args.read_len_min, args.read_len_max
    read_len_range = f"{read_len_min}-{read_len_max}"
    mapping_type = args.mapping
    sample_names = args.names.split()

    for sample in sample_names:
        bam_file_path = base_dir / "4-Aligned" / f"{sample}.bam"
        bam_file = pysam.AlignmentFile(bam_file_path, "rb")

        forward_raw_csv = (
            assign_dir / f"{sample}_{mapping_type}-end_{read_len_range}_raw_For.txt"
        )
        reverse_raw_csv = (
            assign_dir / f"{sample}_{mapping_type}-end_{read_len_range}_raw_Rev.txt"
        )
        forward_rpm_csv = (
            assign_dir / f"{sample}_{mapping_type}-end_{read_len_range}_rpm_For.txt"
        )
        reverse_rpm_csv = (
            assign_dir / f"{sample}_{mapping_type}-end_{read_len_range}_rpm_Rev.txt"
        )

        output_hdf5 = assign_dir / f"{sample}_{mapping_type}-end_{read_len_range}.h5"

        log_file_path = (
            reports_dir / f"{sample}_{mapping_type}-end_{read_len_range}.log"
        )
        log_file = open(log_file_path, "wt")

        total_reads, reads_mapped_once, reads_mapped_twice = 0, 0, 0
        forward_summary, reverse_summary = pd.DataFrame(), pd.DataFrame()

        report = "\nBamFile: {}\nrlmin: {}\nrlmax: {}\nName: {}\nMapping: {}".format(
            bam_file_path, read_len_min, read_len_max, sample, mapping_type
        )
        log_file.write(report + "\n")

        for chromosome in yeast_chromosomes():
            chrom_reads, forward_counts, reverse_counts = (
                0,
                defaultdict(list),
                defaultdict(list),
            )

            for read in bam_file.fetch(chromosome):
                total_reads += 1
                chrom_reads += 1
                read_length = read.query_length
                nh_tag = read.get_tag("NH")

                if nh_tag == 1:
                    reads_mapped_once += 1
                elif nh_tag == 2:
                    reads_mapped_twice += 1

                if nh_tag == 1 or (nh_tag == 2 and args.mapped_twice):
                    if not read.is_reverse:
                        start = read.reference_start
                        end = read.reference_end - 1
                    else:
                        start = read.reference_end - 1
                        end = read.reference_start

                    position = start if mapping_type == "5" else end
                    if read.is_reverse:
                        reverse_counts[read_length].append(position)
                    else:
                        forward_counts[read_length].append(position)

            # Initialize a default value for missing keys
            default_value = [0]
            # Populate Forward and Reverse dictionaries for sequence lengths within the specified range
            for length in range(read_len_min, read_len_max + 1):
                # Populate forward dictionary with counts for the current read length, defaulting to [0] if missing
                forward_counts[length] = Counter(
                    forward_counts.get(length, default_value)
                )

                # Populate reverse dictionary with counts for the current read length, defaulting to [0] if missing
                reverse_counts[length] = Counter(
                    reverse_counts.get(length, default_value)
                )

            forward_summary = pd.concat(
                [
                    forward_summary,
                    update_dataframe(pd.DataFrame(forward_counts), chromosome, "+"),
                ],
                ignore_index=True,
            )
            reverse_summary = pd.concat(
                [
                    reverse_summary,
                    update_dataframe(pd.DataFrame(reverse_counts), chromosome, "-"),
                ],
                ignore_index=True,
            )

            log_file.write(f"{chromosome}: {chrom_reads} reads processed.\n")

        forward_summary.rename(columns=str, inplace=True)
        reverse_summary.rename(columns=str, inplace=True)

        # Save Raw data
        if save_csv:
            save_to_csv(forward_summary, forward_raw_csv)
            save_to_csv(reverse_summary, reverse_raw_csv)

        with pd.HDFStore(output_hdf5, complevel=5, complib="zlib", mode="w") as store:
            store.put("Forward_Raw", forward_summary, format="table", data_columns=True)
            store.put("Reverse_Raw", reverse_summary, format="table", data_columns=True)

        report = "\nTotal No of reads {:>11,} mapped to genome\n".format(total_reads)
        report += "Number of reads {:>11,d} mapped once to genome\n".format(
            reads_mapped_once
        )
        report += "Number of reads {:>11,d} mapped twice to genome\n".format(
            reads_mapped_twice
        )
        log_file.write(report + "\n")

        # Normalize and save RPM data
        normalization_factor = (
            sum(
                1
                for read in bam_file.fetch()
                if read.get_tag("NH") <= (2 if args.mapped_twice else 1)
            )
            / 1e6
        )
        for col in list(map(str, range(read_len_min, read_len_max + 1))) + ["sum"]:
            forward_summary[col] /= normalization_factor
            reverse_summary[col] /= normalization_factor

        if save_csv:
            save_to_csv(forward_summary, forward_rpm_csv)
            save_to_csv(reverse_summary, reverse_rpm_csv)

        with pd.HDFStore(output_hdf5, complevel=5, complib="zlib", mode="a") as store:
            store.put("Forward_RPM", forward_summary, format="table", data_columns=True)
            store.put("Reverse_RPM", reverse_summary, format="table", data_columns=True)

        log_file.close()
        bam_file.close()


def main():
    args = parse_arguments()
    process_raw_assignment(args)


if __name__ == "__main__":
    main()
