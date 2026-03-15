import typer
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re

app = typer.Typer()

# Define all IUPAC nucleotide codes, including degenerate bases
nucleotide_replacements = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


def generate_reverse_complements(search_seqs):
    reverse_complements = []
    for seq in search_seqs:
        seq_obj = Seq(seq)  # Correctly create a Seq object
        reverse_complements.append(str(seq_obj.reverse_complement()))
    return reverse_complements


def search_sequences(genome_path, search_seqs):
    results = []
    # Create a mapping of reverse complements back to the original sequence
    original_to_reverse = {
        seq: str(Seq(seq).reverse_complement()) for seq in search_seqs
    }
    # Include both original sequences and their reverse complements in the search
    search_seqs_extended = list(search_seqs) + list(original_to_reverse.values())
    regex_to_original = {}

    for seq in search_seqs_extended:
        regex_version = "".join([nucleotide_replacements.get(nuc, nuc) for nuc in seq])
        # Map both original and reverse complement sequences back to the original sequence
        original_seq = next(
            (
                orig
                for orig, rev in original_to_reverse.items()
                if seq == orig or seq == rev
            ),
            seq,
        )
        regex_to_original[regex_version] = original_seq

    for record in SeqIO.parse(genome_path, "fasta"):
        for search_seq_regex, original_seq in regex_to_original.items():
            for match in re.finditer(search_seq_regex, str(record.seq)):
                results.append(
                    {
                        "contig_id": record.id,
                        "Oligo": original_seq,  # Use the original sequence for matching
                        "Start": match.start(),
                        "End": match.end(),
                    }
                )

    return pd.DataFrame(results)


def get_genome_info(genome_path):
    # assuming there is only a single chromosome
    basename = os.path.basename(genome_path).split(".")[0]
    for record in SeqIO.parse(genome_path, "fasta"):
        genome_size = len(record.seq) / 1_000_000
        gc_content = SeqUtils.GC(record.seq)
    return basename, genome_size, gc_content


def draw_circular_genome_with_info(dataframe, genome_path):
    basename, genome_size, gc_content = get_genome_info(genome_path)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, aspect="equal")
    genome_circle = plt.Circle(
        (0.5, 0.5), 0.4, edgecolor="#b3b3b3", facecolor="none", linewidth=15
    )
    ax.add_artist(genome_circle)
    ax.text(
        0.5,
        0.5,
        f"{basename}\n{genome_size:.2f} Mb\nGC%: {gc_content:.2f}",
        ha="center",
        va="center",
        fontsize=20,
    )
    genome_length = genome_size * 1_000_000

    positions_normalized = dataframe["Start"] / genome_length
    positions_angles = positions_normalized * 360

    for angle in positions_angles:
        angle_rad = np.deg2rad(angle)
        start_x, start_y = 0.5 + 0.4 * np.cos(angle_rad), 0.5 + 0.4 * np.sin(angle_rad)
        end_x, end_y = 0.5 + 0.45 * np.cos(angle_rad), 0.5 + 0.45 * np.sin(angle_rad)
        ax.plot([start_x, end_x], [start_y, end_y], color="red", lw=1)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    plt.show()


def calculate_distribution_score(positions, genome_length):
    # Assuming positions are sorted
    distances = [positions[i + 1] - positions[i] for i in range(len(positions) - 1)]
    # Include distances to genome start and end for completeness
    distances.append(positions[0] + genome_length - positions[-1])
    # Score could be the minimum distance, to maximize the minimum spacing
    return min(distances)


def select_distributed_oligos(
    dataframe, genome_length, noligos, oligo_match_counts=None
):
    # Group by Oligo and aggregate Start positions
    grouped = dataframe.groupby("Oligo")["Start"].apply(list).reset_index()
    grouped["Start"] = grouped["Start"].apply(lambda x: sorted(x))

    # Calculate distribution scores for each oligo
    grouped["DistributionScore"] = grouped["Start"].apply(
        lambda x: calculate_distribution_score(x, genome_length)
    )

    # Add non-target match counts to the dataframe, if provided
    if oligo_match_counts:
        grouped["NonTargetMatches"] = grouped["Oligo"].apply(
            lambda x: oligo_match_counts.get(x, 0)
        )
        # Prefer oligos with fewer non-target matches
        grouped = grouped[
            grouped["NonTargetMatches"] <= grouped["NonTargetMatches"].quantile(0.1)
        ]

    # Select n oligos with the best (highest) distribution scores
    selected = grouped.nlargest(noligos, "DistributionScore")

    # Get all rows from the original dataframe that match the selected oligos
    selected_oligos_df = dataframe[dataframe["Oligo"].isin(selected["Oligo"])]

    return selected_oligos_df


def summarize_oligos(df, oligo_match_counts):
    # Group by 'Oligo' and aggregate 'Start' and 'End' into a list
    grouped = (
        df.groupby("Oligo")
        .agg(
            {
                "Start": lambda x: list(sorted(x)),
                "End": lambda x: list(sorted(x)),
                "og_id": "first",
                "gene_id": "first",
                "sequence": "first",
                "tm": "first",
                "gc%": "first",
                "conservation": "first",
                "composite_score": "first",
                "nd_perc": "first",
            }
        )
        .reset_index()
    )

    # Create position ranges by pairing Start and End values
    grouped["Position"] = grouped.apply(
        lambda row: ";".join(
            [f"{start}-{end}" for start, end in zip(row["Start"], row["End"])]
        ),
        axis=1,
    )
    # Drop the intermediate Start and End columns
    grouped = grouped.drop(["Start", "End"], axis=1)

    # Add a new column for the number of non-target matches
    grouped["#nontarget-matches"] = grouped["Oligo"].apply(
        lambda x: oligo_match_counts.get(x, 0)
    )

    # Calculate '#matches' based on the Position column
    grouped["#matches"] = grouped["Position"].apply(lambda x: len(x.split(";")))

    return grouped


def test_oligos_against_nontarget(non_target_path, search_seqs):
    # Read the non-target sequences efficiently
    non_target_sequences = []
    for record in SeqIO.parse(non_target_path, "fasta"):
        non_target_sequences.append(str(record.seq))

    # Join all sequences once to avoid repeated concatenation
    non_target_seq = "".join(non_target_sequences)

    # Count matches for each oligo in the non-target sequence
    oligo_match_counts = {}
    for seq in search_seqs:
        regex_version = "".join([nucleotide_replacements.get(nuc, nuc) for nuc in seq])
        matches = re.findall(regex_version, non_target_seq)
        oligo_match_counts[seq] = len(matches)

    return oligo_match_counts


@app.command()
def main(
    genome_path: str = typer.Argument(
        ...,
        help="Path to the genome file.",
    ),
    df_initial_path: str = typer.Argument(..., help="Path to the initial CSV file."),
    tm_low: float = typer.Option(
        30, "--tm-low", help="Low threshold for melting temperature."
    ),
    tm_high: float = typer.Option(
        35, "--tm-high", help="High threshold for melting temperature."
    ),
    length_low: int = typer.Option(
        9, "--length-low", help="Low threshold for oligo length."
    ),
    length_high: int = typer.Option(
        12, "--length-high", help="High threshold for oligo length."
    ),
    conservation_low: float = typer.Option(
        80, "--conservation-low", help="Low threshold for conservation."
    ),
    conservation_high: float = typer.Option(
        100, "--conservation-high", help="High threshold for conservation."
    ),
    nd_perc_low: float = typer.Option(
        70, "--nd-perc-low", help="Low threshold for nd_perc."
    ),
    nd_perc_high: float = typer.Option(
        100, "--nd-perc-high", help="High threshold for nd_perc."
    ),
    noligos: int = typer.Option(10, "--noligos", help="Number of oligos to select."),
    non_target_path: str = typer.Option(
        "", "--non-target-path", help="Path to the non-target sequence file."
    ),
    output_file: str = typer.Option(
        "selected_oligos.csv",
        "--output-file",
        help="Path to save the summarized oligos CSV.",
    ),
):

    df_initial = pd.read_csv(df_initial_path, sep=",")
    # Filter the DataFrame based on specified criteria
    filtered_df_initial = df_initial[
        (df_initial["tm"] >= tm_low)
        & (df_initial["tm"] <= tm_high)
        & (df_initial["length"] >= length_low)
        & (df_initial["length"] <= length_high)
        & (df_initial["conservation"] >= conservation_low)
        & (df_initial["conservation"] <= conservation_high)
        & (df_initial["nd_perc"] >= nd_perc_low)
        & (df_initial["nd_perc"] <= nd_perc_high)
    ]

    # Extract oligo sequences from the filtered DataFrame
    search_seqs = filtered_df_initial["sequence"].unique().tolist()

    # Retrieve genome information
    basename, genome_size_mb, gc_content = get_genome_info(genome_path)
    genome_length_bp = (
        genome_size_mb * 1_000_000
    )  # Convert Mb to base pairs for detailed calculations

    # Proceed with searching these sequences in the target genome
    matched_sequences_df = search_sequences(genome_path, search_seqs)

    # Filtering and sequence extraction logic remains the same

    # Test oligos against the non-target sequence if provided
    oligo_match_counts = {}
    if non_target_path:
        oligo_match_counts = test_oligos_against_nontarget(non_target_path, search_seqs)

    selected_oligos_df = select_distributed_oligos(
        matched_sequences_df, genome_length_bp, noligos, oligo_match_counts
    )

    # Continue with plotting, summarizing, and saving
    draw_circular_genome_with_info(selected_oligos_df, genome_path)
    summarized_oligos_df = summarize_oligos(selected_oligos_df, oligo_match_counts)
    summarized_oligos_df.to_csv(output_file, sep="\t", index=None)


if __name__ == "__main__":
    app()
