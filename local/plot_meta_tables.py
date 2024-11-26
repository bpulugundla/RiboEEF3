import os
import argparse
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import set_matplotlib_formats

import warnings
import tables

warnings.filterwarnings("ignore", category=tables.NaturalNameWarning)


def parse_args():
    """
    Parse command-line arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description=("This script generates meta-gene plots from meta-gene tables")
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
        "--mapping",
        type=str,
        default="5",
        choices=["5", "3"],
        help="Map reads according to their 5' or 3' end",
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


def trim_df(df, span, position, inside_gene=33, outside_gene=18):
    """
    Trim a DataFrame to fit in a 5' prime mapping figure.

    Parameters:
    - df: Input pandas DataFrame indexed by relative positions ('rel_Pos').
    - span: Integer, the range of positions to include.
    - position: String, either 'Start' or 'Stop', indicating gene boundary type.
    - inside_gene: Integer, number of positions to retain inside the gene boundary.
    - outside_gene: Integer, number of positions to retain outside the gene boundary.

    Returns:
    - Trimmed pandas DataFrame.
    """
    if inside_gene > span or outside_gene > span:
        print(
            "Warning: Inside- or outside-gene parameters exceed the span! "
            "Querying out-of-range data."
        )
        return df

    if position == "Start":
        return df.loc[-outside_gene:inside_gene,]
    elif position == "Stop":
        return df.loc[-inside_gene:outside_gene,]
    else:
        print("Table is not modified. Unknown mapping type!")
        return df


def generate_metag_plot_pdf(args):
    """
    Generate metagene plots and save as PDF files.
    """

    # Set output formats for matplotlib
    set_matplotlib_formats("pdf", "svg")

    base_dir = Path(args.rawdir).resolve().parent
    metagene_dir = base_dir / "6-MetaGeneTables"
    plots_dir = base_dir / "7-MetaGenePlots"
    os.makedirs(plots_dir, exist_ok=True)

    # Set seaborn style for aesthetics
    sns.set_style("white")
    sns.set_context("paper")

    # Extract parameters
    read_len_min, read_len_max = args.read_len_min, args.read_len_max
    span, mapping, normalization = args.metagene_span, args.mapping, args.normalised
    read_len_range = f"{read_len_min}-{read_len_max}"
    read_lengths = [str(i) for i in range(read_len_min, read_len_max + 1)] + ["sum"]

    # Define plot colors for different read lengths
    color_map = {
        "25": "fuchsia",
        "26": "blueviolet",
        "27": "darkblue",
        "28": "b",
        "29": "r",
        "30": "salmon",
        "31": "orange",
        "32": "olive",
        "33": "g",
        "34": "tan",
        "35": "y",
        "sum": "brown",
    }

    sample_names = args.names.split()
    for sample in sample_names:
        for position_type in ["Start", "Stop"]:
            # Input and output file paths
            input_file = (
                f"{metagene_dir}/{sample}_{mapping}-end_{read_len_range}_"
                f"{normalization}_{position_type}_Meta_Sum.txt"
            )
            output_file = (
                f"{plots_dir}/{sample}-{mapping}-end-{read_len_range}-"
                f"{normalization}-{position_type}.pdf"
            )

            plot_title = (
                f"{sample.replace('_', '-')} {position_type} {mapping}' mapping"
            )
            legend_location = "upper right" if position_type == "Stop" else "upper left"

            if os.path.isfile(input_file):
                # Configure figure dimensions
                fig_width = 8
                fig_height = 1.2 * len(read_lengths)
                fig, axes = plt.subplots(
                    nrows=len(read_lengths), figsize=(fig_width, fig_height)
                )
                fig.suptitle(plot_title, y=0.9, fontsize=12)

                # Load and preprocess data
                data_frame = pd.read_csv(input_file, index_col=0, sep="\t")
                data_frame.set_index("rel_pos", inplace=True)

                # Trim data based on mapping type and position
                data_frame = trim_dataframe(data_frame, mapping, position_type, span)

                # Generate subplots for each read length
                for i, read_length in enumerate(read_lengths):
                    alpha_value = 0.6
                    color_map = validate_color_map(color_map, read_length)
                    x_values = data_frame.index
                    y_values = data_frame[read_length]

                    # Create bar plot
                    axes[i].bar(
                        x_values,
                        y_values,
                        color=color_map[read_length],
                        alpha=alpha_value,
                    )
                    axes[i].legend([read_length], loc=legend_location)

                    # Add guide lines based on mapping
                    add_guidelines(axes[i], x_values, mapping, position_type)

                    # Set y-axis label
                    axes[i].set_ylabel(normalization)

                # Finalize plot aesthetics
                sns.despine()
                fig.savefig(output_file, format="pdf", dpi=300, bbox_inches="tight")
                print(f"Saved: {output_file}")
            else:
                print(f"Missing input file: {input_file}")


def trim_dataframe(data_frame, mapping, position_type, span):
    """
    Trim the data frame based on the mapping type and position.

    Parameters:
    - data_frame: Input pandas DataFrame
    - mapping: Mapping type ('5' or '3')
    - position_type: Position type ('Start' or 'Stop')
    - span: Integer span value

    Returns:
    - Trimmed pandas DataFrame
    """
    if mapping == "5":
        if position_type == "Start":
            return trim_df(
                data_frame, span, position_type, inside_gene=39, outside_gene=21
            )
        elif position_type == "Stop":
            return trim_df(
                data_frame, span, position_type, inside_gene=60, outside_gene=3
            )
    elif mapping == "3":
        if position_type == "Start":
            return trim_df(
                data_frame, span, position_type, inside_gene=60, outside_gene=3
            )
        elif position_type == "Stop":
            return trim_df(
                data_frame, span, position_type, inside_gene=39, outside_gene=30
            )
    return data_frame


def validate_color_map(color_map, read_length):
    """
    Validate and update color map for the given read length.
    """
    if read_length not in color_map:
        color_map[read_length] = "gray"  # Default color if missing
    return color_map


def add_guidelines(axis, x_values, mapping, position_type):
    """
    Add guide lines to the plot based on mapping type and position.

    Parameters:
    - axis: Matplotlib axis object
    - x_values: Range of x values in the plot
    - mapping: Mapping type ('5' or '3')
    - position_type: Position type ('Start' or 'Stop')
    """
    min_x, max_x = x_values.min(), x_values.max()
    for position in range(min_x, max_x + 1, 3):
        color = "gray"
        alpha = 0.2
        if mapping == "5" and position_type == "Start" and position == -12:
            color, alpha = "g", 0.5
        elif mapping == "5" and position_type == "Start" and position == 0:
            color, alpha = "r", 0.4
        elif mapping == "3" and position_type == "Stop" and position == 12:
            color, alpha = "g", 0.5
        elif mapping == "3" and position_type == "Stop" and position == 0:
            color, alpha = "r", 0.4
        axis.axvline(x=position, linewidth=1, alpha=alpha, color=color)


def main():
    args = parse_args()
    generate_metag_plot_pdf(args)


if __name__ == "__main__":
    main()
