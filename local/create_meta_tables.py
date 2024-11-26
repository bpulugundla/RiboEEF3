import os
import pandas as pd
import argparse
from pathlib import Path
import pysam

from process_raw_assignment import yeast_chromosomes

import warnings
import tables

warnings.filterwarnings("ignore", category=tables.NaturalNameWarning)


def parse_args():
    """
    Parse command-line arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description=(
            "This script generates meta-gene tables from sequencing data. "
            "It aggregates data around start and stop codons for a set of genes "
            "using various parameters such as read lengths, span, and normalization."
        )
    )

    parser.add_argument("--rawdir", type=str, help="Raw FASTQ file directory")
    parser.add_argument("--names", type=str, help="Names of samples")
    parser.add_argument(
        "--metagene_span",
        type=int,
        default=60,
        help="Number of nucleotides before and after start/stop",
    )
    parser.add_argument(
        "--mapped_twice",
        type=int,
        default=0,
        choices=[0, 1],
        help="Include reads mapped twice (1) or only once (0).",
    )
    parser.add_argument(
        "--mapping",
        type=str,
        default="5",
        choices=["5", "3"],
        help="Map reads according to their 5' or 3' end",
    )
    parser.add_argument(
        "--metagene_threshold",
        type=int,
        default=30,
        help="Number of nucleotides before and after start/stop",
    )
    parser.add_argument(
        "--read_len_min", type=int, default=25, help="Minimum length of reads"
    )
    parser.add_argument(
        "--read_len_max", type=int, default=35, help="Maximum length of reads"
    )
    parser.add_argument(
        "--normalised",
        type=str,
        default="rpm",
        choices=["raw", "rpm"],
        help="Normalization type: 'raw' or 'rpm'",
    )
    return parser.parse_args()


def not_enough_data(df, threshold=12):
    """
    Check if the DataFrame has sufficient data based on the threshold.
    """
    return True if df.sum().sum() < threshold else False


def reads_count_in_bam(bam_name, mapped_twice):
    """
    Count reads in a BAM file based on mapping criteria.
    """
    with pysam.AlignmentFile(bam_name, "rb") as bamfile:
        if mapped_twice:
            count = sum(1 for read in bamfile.fetch() if read.get_tag("NH") <= 2)
            report = f"No of reads mapped once and twice: {count:,}"
        else:
            count = sum(1 for read in bamfile.fetch() if read.get_tag("NH") == 1)
            report = f"No of reads mapped once: {count:,}"
    return count


def normalisation_factor_from_bam(bam_name, mapped_twice):
    """
    Calculate the normalization factor from a BAM file.
    """
    return reads_count_in_bam(bam_name, mapped_twice) / 1_000_000


def raw_metag_threshold_to_rpm(bam_name, threshold, mapped_twice):
    """
    Convert a raw threshold to RPM-normalized threshold.
    """
    normalization_factor = normalisation_factor_from_bam(bam_name, mapped_twice)
    return threshold / normalization_factor


def process_bam_metrics(bam_name, mapped_twice, raw_threshold=None):
    """
    Compute read count, normalization factor, and optionally convert a raw threshold to RPM.

    Parameters:
        bam_name: Path to the BAM file
        mapped_twice: Whether to include reads mapped up to twice
        raw_threshold: Optional raw threshold to convert to RPM (default: None)

    Returns:
        result: A dictionary containing read count, normalization factor, and converted threshold if provided
    """
    bamfile = pysam.AlignmentFile(bam_name, "rb")  # Open BAM file

    # Calculate read count
    read_count = sum(
        1
        for read in bamfile.fetch()
        if read.get_tag("NH") <= (2 if mapped_twice else 1)
    )

    normalization_factor = read_count / (10**6)  # Calculate normalization factor

    result = {"read_count": read_count, "normalization_factor": normalization_factor}

    # Convert raw threshold to RPM if provided
    if raw_threshold is not None:
        rpm_threshold = raw_threshold / normalization_factor
        result["rpm_threshold"] = rpm_threshold

    return result


def df_framing(df, index, columns, strand="+"):
    """
    Frame the DataFrame with all positions in the specified range, filling missing values.
    """
    full_df = pd.DataFrame(0, index=index, columns=columns)
    df = df.add(full_df, fill_value=0, axis=1)
    df.reset_index(inplace=True)
    return df if strand == "+" else df.iloc[::-1]


# Function to process codons for a given strand and feature
def process_codon(
    gtf, df_source, meta_df, count, span, threshold, columns, strand_type
):
    """
    Process codons for given feature, strand, and parameters.

    Parameters:
        gtf: GTF record for the feature.
        df_source: DataFrame to extract data from.
        meta_df: Accumulator DataFrame to store results.
        count: Counter to track processed features.
        span: Number of nucleotides around the feature to include.
        threshold: Minimum data threshold to proceed.
        columns: Relevant columns to extract.
        strand_type: Strand direction ('+' or '-').

    Returns:
        Updated meta_df and count.
    """
    # Determine start and index range based on strand type
    if strand_type == "+":
        start = gtf.start
        index = range(start - span, start + span + 1)
    elif strand_type == "-":
        start = gtf.end - 1  # -1 correction for reverse strand
        index = range(gtf.end - span - 1, gtf.end + span)

    # Extract a sub-dataframe around the feature
    df = df_source[start - span : start + span][columns].copy()

    # Skip if there's not enough data
    if not_enough_data(df, threshold):
        return meta_df, count

    # Reframe and align the dataframe with the given index
    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)

    # Accumulate results and increment the count
    meta_df += df
    count += 1

    return meta_df, count


def create_meta_gene_tables(args):
    """
    Generate meta-gene tables based on the provided arguments.
    """
    base_dir = Path(args.rawdir).resolve().parent
    metagene_dir = base_dir / "6-MetaGeneTables"
    reports_dir = metagene_dir / "Reports"
    os.makedirs(metagene_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    read_len_min, read_len_max = args.read_len_min, args.read_len_max
    span, mapping, normalization = args.metagene_span, args.mapping, args.normalised

    columns = [str(i) for i in range(read_len_min, read_len_max + 1)] + ["sum"]
    read_len_range = f"{read_len_min}-{read_len_max}"
    log_file_path = reports_dir / f"MetaGeneTables_{mapping}-end_{read_len_range}.log"

    with open(log_file_path, "wt") as log_file:
        sample_names = args.names.split()
        for sample in sample_names:
            cf1 = cr1 = cf2 = cr2 = 0  # counters
            log_file.write(f"\nProcessing sample: {sample}\n")

            # File paths
            filename = f"{sample}_{mapping}-end_{read_len_range}"
            output_start = (
                metagene_dir / f"{filename}_{normalization}_Start_Meta_Sum.txt"
            )
            output_stop = metagene_dir / f"{filename}_{normalization}_Stop_Meta_Sum.txt"
            input_h5 = base_dir / f"5-AssignRaw/{filename}.h5"

            # Initialize empty meta-gene DataFrames
            meta_start_f = pd.DataFrame(0, index=range(2 * span + 1), columns=columns)
            meta_start_r = pd.DataFrame(0, index=range(2 * span + 1), columns=columns)
            meta_stop_f = pd.DataFrame(0, index=range(2 * span + 1), columns=columns)
            meta_stop_r = pd.DataFrame(0, index=range(2 * span + 1), columns=columns)

            # Threshold adjustment for normalization
            threshold = args.metagene_threshold
            if normalization == "rpm":
                bam_file = base_dir / f"4-Aligned/{sample}.bam"
                bam_metrics = process_bam_metrics(
                    bam_file, args.mapped_twice, threshold
                )
                threshold = bam_metrics["rpm_threshold"]

            log_file.write(f"Read Count: {bam_metrics['read_count']}")
            log_file.write(
                f"Normalization Factor: {bam_metrics['normalization_factor']}"
            )
            if "rpm_threshold" in bam_metrics:
                log_file.write(f"RPM Threshold: {bam_metrics['rpm_threshold']}")

            # Process chromosomes and GTF annotations
            tabix_file = pysam.TabixFile(
                "References/genome.gtf.gz", parser=pysam.asGTF()
            )
            df_forward = pd.read_hdf(
                input_h5, "Forward_RPM" if normalization == "rpm" else "Forward_Raw"
            )
            df_reverse = pd.read_hdf(
                input_h5, "Forward_RPM" if normalization == "rpm" else "Forward_Raw"
            )

            for chromosome in yeast_chromosomes():
                # Filter by chromosome
                chr_f_df = df_forward[df_forward.Chromosome == chromosome].set_index(
                    "Position"
                )[columns]
                chr_r_df = df_reverse[df_reverse.Chromosome == chromosome].set_index(
                    "Position"
                )[columns]

                for gtf in tabix_file.fetch(reference=chromosome):
                    # Stop codon processing
                    if gtf.feature == "stop_codon":
                        if gtf.strand == "+":  # Forward strand
                            meta_stop_f, cf2 = process_codon(
                                gtf=gtf,
                                df_source=chr_f_df,
                                meta_df=meta_stop_f,
                                count=cf2,
                                span=span,
                                threshold=threshold,
                                columns=columns,
                                strand_type="+",
                            )
                        elif gtf.strand == "-":  # Reverse strand
                            meta_stop_r, cr2 = process_codon(
                                gtf=gtf,
                                df_source=chr_r_df,
                                meta_df=meta_stop_r,
                                count=cr2,
                                span=span,
                                threshold=threshold,
                                columns=columns,
                                strand_type="-",
                            )

                    # Start codon processing
                    elif gtf.feature == "start_codon":
                        if gtf.strand == "+":  # Forward strand
                            meta_start_f, cf1 = process_codon(
                                gtf=gtf,
                                df_source=chr_f_df,
                                meta_df=meta_start_f,
                                count=cf1,
                                span=span,
                                threshold=threshold,
                                columns=columns,
                                strand_type="+",
                            )
                        elif gtf.strand == "-":  # Reverse strand
                            meta_start_r, cr1 = process_codon(
                                gtf=gtf,
                                df_source=chr_r_df,
                                meta_df=meta_start_r,
                                count=cr1,
                                span=span,
                                threshold=threshold,
                                columns=columns,
                                strand_type="-",
                            )

            # Summing up and saving results
            meta_start_sum = meta_start_f + meta_start_r
            meta_stop_sum = meta_stop_f + meta_stop_r

            meta_start_sum["rel_pos"] = list(range(-span, span + 1))
            meta_start_sum.to_csv(output_start, sep="\t", header=True, index=True)

            meta_stop_sum["rel_pos"] = list(range(-span, span + 1))
            meta_stop_sum.to_csv(output_stop, sep="\t", header=True, index=True)

            log_file.write(f"Meta-gene tables saved for {sample}.\n")

        log_file.close()


def main():
    args = parse_args()
    create_meta_gene_tables(args)


if __name__ == "__main__":
    main()
