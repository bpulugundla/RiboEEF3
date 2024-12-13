#!/usr/bin/env python3
"""
Author: Bhargav Pulugundla
Last Updated: 13-Dec-2024
"""

import argparse
import csv
import pysam
from multiprocessing import Pool, cpu_count


def parse_gtf(gtf_file):
    """
    Parse a GTF file to extract start and stop codon positions.
    Handles lines specifying 'start_codon' and 'stop_codon' features.

    Args:
        gtf_file (str): Path to the GTF file.

    Returns:
        tuple: Lists of start codon positions and stop codon positions.
    """
    start_codons = []
    stop_codons = []

    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue  # Skip header lines
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # Skip invalid lines

            chrom, feature, start, end, strand = (
                fields[0],
                fields[2],
                int(fields[3]),
                int(fields[4]),
                fields[6],
            )

            if feature == "start_codon":
                start_pos = start if strand == "+" else end
                start_codons.append((chrom, start_pos, strand))
            elif feature == "stop_codon":
                stop_pos = start if strand == "+" else end
                stop_codons.append((chrom, stop_pos, strand))

    return start_codons, stop_codons


def get_total_mapped_reads(bam_file):
    """
    Get the total number of mapped reads from a BAM file.

    Args:
        bam_file (str): Path to the BAM file.

    Returns:
        int: Total number of mapped reads.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return bam.mapped


def count_reads_for_codon(codon_data):
    """
    Count reads for a single codon position.

    Args:
        codon_data (tuple): Data containing chromosome, position, strand, range, and BAM file path.

    Returns:
        dict: Read counts at relative positions for the codon.
    """
    chrom, position, strand, rel_range, bam_file = codon_data
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        counts = {}
        for rel_pos in range(rel_range[0], rel_range[1] + 1):
            pos = position + rel_pos if strand == "+" else position - rel_pos
            counts[rel_pos] = bam.count(chrom, pos, pos + 1)
    return counts


def count_reads_at_relpos(bam_file, codons, ranges):
    """
    Count reads at relative positions for start and stop codons using parallel processing.

    Args:
        bam_file (str): Path to the BAM file.
        codons (dict): Codon data with 'start_codon' and 'stop_codon' keys.
        ranges (dict): Relative ranges for start and stop codons.

    Returns:
        dict: Read counts at relative positions for all codons.
    """
    tasks = [
        (chrom, position, strand, ranges[codon_type], bam_file)
        for codon_type, codon_list in codons.items()
        for chrom, position, strand in codon_list
    ]

    with Pool(processes=cpu_count()) as pool:
        results = pool.map(count_reads_for_codon, tasks)

    read_counts = {
        key: {rel_pos: 0 for rel_pos in range(ranges[key][0], ranges[key][1] + 1)}
        for key in codons
    }
    idx = 0
    for codon_type, codon_list in codons.items():
        for _ in codon_list:
            counts = results[idx]
            for rel_pos, count in counts.items():
                read_counts[codon_type][rel_pos] += count
            idx += 1

    return read_counts


def normalize_counts_to_rpm(read_counts, total_reads):
    """
    Normalize raw counts to Reads Per Million (RPM).

    Args:
        read_counts (dict): Raw read counts.
        total_reads (int): Total number of mapped reads.

    Returns:
        dict: Normalized read counts.
    """
    for codon_type in read_counts:
        for rel_pos in read_counts[codon_type]:
            read_counts[codon_type][rel_pos] = (
                read_counts[codon_type][rel_pos] / total_reads * 1e6
            )
    return read_counts


def save_counts_to_csv(read_counts, output_file):
    """
    Save normalized read counts to a CSV file.

    Args:
        read_counts (dict): Normalized read counts.
        output_file (str): Path to the output CSV file.
    """
    with open(output_file, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Codon Type", "Relative Position", "Counts RPM"])
        for codon_type, counts in read_counts.items():
            for rel_pos, count in sorted(counts.items()):
                writer.writerow([codon_type, rel_pos, count])


def main():
    """
    Main function to process GTF and BAM files and generate normalized read counts.
    """
    parser = argparse.ArgumentParser(
        description="Process BAM and GTF files to count reads at codon positions."
    )
    parser.add_argument("bam_file", type=str, help="Path to the BAM file.")
    parser.add_argument("gtf_file", type=str, help="Path to the GTF file.")
    parser.add_argument("output_csv", type=str, help="Path to the output CSV file.")
    args = parser.parse_args()

    # Define relative ranges for codons
    relative_ranges = {"start_codon": (-20, 40), "stop_codon": (-40, 20)}

    # Parse GTF file
    start_codons, stop_codons = parse_gtf(args.gtf_file)
    codons = {"start_codon": start_codons, "stop_codon": stop_codons}

    # Get total reads from BAM file
    total_reads = get_total_mapped_reads(args.bam_file)

    # Count reads for codons
    read_counts = count_reads_at_relpos(args.bam_file, codons, relative_ranges)

    # Normalize counts and save
    read_counts = normalize_counts_to_rpm(read_counts, total_reads)
    save_counts_to_csv(read_counts, args.output_csv)

    print(f"Normalized read counts (RPM) have been saved to {args.output_csv}")


if __name__ == "__main__":
    main()
