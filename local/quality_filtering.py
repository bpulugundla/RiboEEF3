import os
import subprocess
import gzip
import argparse
from pathlib import Path


def parse_args():
    """
    Parse command-line arguments for the script.
    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="A script for quality control of FASTQ sequences"
    )

    parser.add_argument(
        "--rawdir", type=str, help="Directory containing raw FASTQ files."
    )
    parser.add_argument("--names", type=str, help="Sample names separated by space.")
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.80,
        help="Quality score threshold.",
    )
    parser.add_argument(
        "--read_len_min",
        type=int,
        default=25,
        help="Minimum length of reads.",
    )
    parser.add_argument(
        "--read_len_max",
        type=int,
        default=35,
        help="Maximum length of reads.",
    )

    return parser.parse_args()


def quality_filtering(args):
    """
    Filters FASTQ sequences based on PHRED score and read length.

    Args:
        args (argparse.Namespace): Parsed arguments containing filtering parameters.
    """
    # Mapping PHRED quality scores to probabilities
    PHREDDict = {
        "!": 9.999999e-01,
        '"': 7.943282e-01,
        "#": 6.309573e-01,
        "$": 5.011872e-01,
        "%": 3.981072e-01,
        "&": 3.162278e-01,
        "'": 2.511886e-01,
        "(": 1.995262e-01,
        ")": 1.584893e-01,
        "*": 1.258925e-01,
        "+": 1.000000e-01,
        ",": 7.943282e-02,
        "-": 6.309573e-02,
        ".": 5.011872e-02,
        "/": 3.981072e-02,
        "0": 3.162278e-02,
        "1": 2.511886e-02,
        "2": 1.995262e-02,
        "3": 1.584893e-02,
        "4": 1.258925e-02,
        "5": 1.000000e-02,
        "6": 7.943282e-03,
        "7": 6.309573e-03,
        "8": 5.011872e-03,
        "9": 3.981072e-03,
        ":": 3.162278e-03,
        ";": 2.511886e-03,
        "<": 1.995262e-03,
        "=": 1.584893e-03,
        ">": 1.258925e-03,
        "?": 1.000000e-03,
        "@": 7.943282e-04,
        "A": 6.309573e-04,
        "B": 5.011872e-04,
        "C": 3.981072e-04,
        "D": 3.162278e-04,
        "E": 2.511886e-04,
        "F": 1.995262e-04,
        "G": 1.584893e-04,
        "H": 1.258925e-04,
        "I": 1.000000e-04,
        "J": 7.943282e-05,
    }

    # Setting up directories for output and logs
    base_dir = str(Path(args.rawdir).resolve().parent)
    filtered_dir = f"{base_dir}/2-Filtered"
    reports_dir = f"{filtered_dir}/Reports"
    os.makedirs(filtered_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    # Log file to record quality filtering details
    log_file_path = f"{reports_dir}/Quality_filtering_iv_log.txt"
    log_file = open(log_file_path, "wt")

    # Process each sample specified in the `--names` argument
    sample_names = args.names.split()
    for name in sample_names:
        # Initialize counters for filtering statistics
        low_qual, short_reads, long_reads, included_reads = 0, 0, 0, 0

        # Input and output file paths for the sample
        trimmed_file = f"{base_dir}/1-Trimmed/{name}_trimmed.fastq.gz"
        filtered_file = f"{filtered_dir}/{name}_filtered.fastq"

        # Calculate the total number of reads in the sample
        cmd = f"echo $(zcat {trimmed_file} | wc -l)/4 | bc"
        total_reads = int(subprocess.check_output(cmd, shell=True, text=True).strip())

        # Log the total reads for the sample
        log_file.write(f"{name:16}: {total_reads:>12,} reads\n")
        print(f"{name:16}: {total_reads:>12,} reads")

        # Open input (compressed) and output (filtered) files
        with gzip.open(trimmed_file, "rt") as infile, open(
            filtered_file, "w"
        ) as outfile:
            for _ in range(total_reads):
                # Read one record (4 lines) from the FASTQ file
                identifier = infile.readline().strip()
                sequence = infile.readline().strip()
                q_identifier = infile.readline().strip()
                phred = infile.readline().strip()

                read_length = len(phred)

                # Classify reads based on length
                if read_length < args.read_len_min:
                    short_reads += 1
                elif read_length > args.read_len_max:
                    long_reads += 1
                else:
                    # Calculate the probability score for the read
                    score = 1.0
                    for phred_character in phred:
                        score *= 1 - PHREDDict[phred_character]

                    # Include or discard the read based on the quality score
                    if score > args.threshold:
                        included_reads += 1
                        outfile.write(
                            f"{identifier}\n{sequence}\n{q_identifier}\n{phred}\n"
                        )
                    else:
                        low_qual += 1

        # Generate a summary report for the sample
        report = (
            f" Reads len < {args.read_len_min:>3}: {short_reads:>12,}\n"
            f" Reads len > {args.read_len_max:>3}: {long_reads:>12,}\n"
            f" Quality < {args.threshold:>3}: {low_qual:>12,}\n"
            f" Reads left : {included_reads:>12,}\n"
        )

        # Log and print the summary
        log_file.write(report + "\n")
        print(report)

    # Close the log file
    log_file.close()


def main():
    args = parse_args()
    quality_filtering(args)


if __name__ == "__main__":
    main()
