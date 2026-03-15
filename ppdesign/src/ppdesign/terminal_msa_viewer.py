#!/usr/bin/env python3
"""Terminal-based MSA viewer for gRNA alignments with rich formatting."""

from pathlib import Path
from typing import List, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from rich.console import Console
from rich.table import Table
from rich.text import Text
from rich.panel import Panel
from rich import box
import typer

app = typer.Typer(help="Terminal-based MSA viewer for gRNA alignments")
console = Console()


def find_grna_matches(
    grna_seq: str, target_seq: str
) -> List[Tuple[int, int, str, str]]:
    """Find all gRNA matches in target sequence."""
    matches = []
    target_seq = target_seq.upper()
    grna_seq = grna_seq.upper()

    # Forward strand
    start = 0
    while True:
        pos = target_seq.find(grna_seq, start)
        if pos == -1:
            break
        # Check PAM
        if pos + len(grna_seq) + 2 <= len(target_seq):
            pam = target_seq[pos + len(grna_seq) : pos + len(grna_seq) + 3]
            if len(pam) >= 2 and pam[1:3] in ["GG", "AG"]:
                matches.append((pos, pos + len(grna_seq), "+", pam))
        start = pos + 1

    # Reverse strand
    rev_grna = str(Seq(grna_seq).reverse_complement())
    start = 0
    while True:
        pos = target_seq.find(rev_grna, start)
        if pos == -1:
            break
        # Check PAM
        if pos >= 3:
            pam = target_seq[pos - 3 : pos]
            if len(pam) >= 2 and pam[0:2] in ["CC", "CT"]:
                matches.append((pos, pos + len(grna_seq), "-", pam))
        start = pos + 1

    return matches


def create_colored_sequence(
    seq: str,
    grna_start: int,
    grna_end: int,
    pam_start: int,
    pam_end: int,
    window: int = 30,
) -> Text:
    """Create a colored sequence with highlighted gRNA and PAM."""
    colored = Text()

    # Add position ruler
    for i, nt in enumerate(seq):
        if grna_start <= i < grna_end:
            # gRNA region - bold red
            colored.append(nt, style="bold red on dark_red")
        elif pam_start <= i < pam_end:
            # PAM region - bold cyan
            colored.append(nt, style="bold cyan on dark_cyan")
        else:
            # Normal sequence - dim
            if nt in "ATCG":
                colored.append(nt, style="dim white")
            else:
                colored.append(nt, style="dim yellow")

    return colored


def create_position_ruler(length: int, start_pos: int = 0, interval: int = 10) -> Text:
    """Create a position ruler for alignment."""
    ruler = Text()
    ruler_line1 = ""
    ruler_line2 = ""

    for i in range(length):
        actual_pos = start_pos + i
        if actual_pos % interval == 0:
            pos_str = str(actual_pos)
            ruler_line1 += pos_str[0] if len(pos_str) > 0 else " "
            ruler_line2 += "|"
        else:
            ruler_line1 += " "
            ruler_line2 += "."

    ruler.append(ruler_line1 + "\n", style="dim cyan")
    ruler.append(ruler_line2, style="dim cyan")
    return ruler


@app.command()
def view(
    grna_csv: Path = typer.Option(
        "results/results/validation_run/guide_rnas.csv",
        "--csv",
        "-c",
        help="Path to guide RNA CSV file",
    ),
    fasta_file: Path = typer.Option(
        "tests/test.fna", "--fasta", "-f", help="Path to FASTA file with sequences"
    ),
    num_grnas: int = typer.Option(
        3, "--num", "-n", help="Number of top gRNAs to display"
    ),
    num_targets: int = typer.Option(
        5, "--targets", "-t", help="Number of target sequences to show per gRNA"
    ),
    window: int = typer.Option(
        40, "--window", "-w", help="Window size around gRNA match"
    ),
):
    """View gRNA alignments in the terminal with colored highlighting."""

    # Load sequences
    console.print("[bold cyan]Loading sequences...[/bold cyan]")
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    # Load gRNA results
    df = pd.read_csv(grna_csv)

    # Load target mapping for complete target list
    target_mapping_file = grna_csv.parent / "target_mapping.tsv"
    if target_mapping_file.exists():
        target_df = pd.read_csv(target_mapping_file, sep="\t")
    else:
        target_df = None

    # Create main title
    console.print(
        Panel.fit(
            "[bold magenta]gRNA Terminal MSA Viewer[/bold magenta]\n"
            + f"[dim]Showing top {num_grnas} gRNAs with {num_targets} targets each[/dim]",
            border_style="bright_blue",
        )
    )

    # Process each gRNA
    for idx in range(min(num_grnas, len(df))):
        row = df.iloc[idx]
        grna_id = row["Guide_ID"]
        grna_seq = row["Sequence"]
        conservation = row["Conservation"]
        quality_score = row["Quality_Score"]

        # Get complete target list from target_mapping if available
        if target_df is not None:
            grna_targets = target_df[target_df["Guide_ID"] == grna_id][
                "Target_Sequence"
            ].unique()
            target_list = list(grna_targets)
        else:
            target_list = row["Target_Sequences"].split(";")

        # Calculate how many targets to actually show
        actual_targets_to_show = min(num_targets, len(target_list))

        # Create gRNA info panel
        info_table = Table(show_header=False, box=None, padding=0)
        info_table.add_column("Field", style="cyan")
        info_table.add_column("Value", style="white")

        info_table.add_row("Sequence:", f"[bold red]{grna_seq}[/bold red]")
        info_table.add_row("PAM:", f"[cyan]{row['PAM']}[/cyan]")
        info_table.add_row(
            "Conservation:",
            f"[green]{conservation}[/green] ({len(target_list)}/{len(sequences)} sequences)",
        )
        info_table.add_row("Quality Score:", f"{quality_score:.3f}")
        info_table.add_row("Matches Found:", f"{len(target_list)} sequences")
        info_table.add_row("Showing:", f"{actual_targets_to_show} sequences")

        console.print(
            Panel(
                info_table,
                title=f"[bold yellow]{grna_id}[/bold yellow]",
                border_style="yellow",
                box=box.ROUNDED,
            )
        )

        # Create alignment table
        alignment_table = Table(
            title="[bold]Sequence Alignments (all shown in forward orientation)[/bold]",
            show_header=True,
            header_style="bold cyan",
            box=box.SIMPLE_HEAD,
            padding=(0, 1),
            collapse_padding=True,
        )

        alignment_table.add_column("Target", style="dim", width=30, overflow="ellipsis")
        alignment_table.add_column("Pos", style="green", justify="right")
        alignment_table.add_column("Strand", style="yellow", justify="center")
        alignment_table.add_column("Alignment", style="white", no_wrap=False)

        # Show alignments for this gRNA
        shown_targets = 0

        # Add note if fewer targets available than requested
        if num_targets > len(target_list):
            console.print(
                f"[yellow]Note: Showing all {len(target_list)} available targets (requested {num_targets})[/yellow]"
            )

        for target_id in target_list[:actual_targets_to_show]:
            if target_id in sequences:
                seq = sequences[target_id]
                matches = find_grna_matches(grna_seq, seq)

                for match in matches[:1]:  # Show first match per target
                    start, end, strand, pam = match

                    # Extract region with window
                    if strand == "+":
                        # Forward strand - extract normally
                        region_start = max(0, start - window)
                        region_end = min(len(seq), end + window + 3)
                        region_seq = seq[region_start:region_end]

                        # Calculate positions in region
                        grna_start_in_region = start - region_start
                        grna_end_in_region = end - region_start
                        pam_start = grna_end_in_region
                        pam_end = min(pam_start + 3, len(region_seq))
                    else:
                        # Reverse strand - need to reverse complement the region
                        # so gRNA appears in forward orientation
                        region_start = max(0, start - 3 - window)  # Include PAM before
                        region_end = min(len(seq), end + window)
                        region_seq = seq[region_start:region_end]

                        # Reverse complement the region
                        region_seq = str(Seq(region_seq).reverse_complement())

                        # After reverse complement, positions are flipped
                        total_len = len(region_seq)
                        grna_start_in_region = total_len - (end - region_start)
                        grna_end_in_region = total_len - (start - region_start)

                        # PAM is now at the end (after gRNA)
                        pam_start = grna_end_in_region
                        pam_end = min(pam_start + 3, len(region_seq))

                    # Create colored sequence
                    colored_seq = create_colored_sequence(
                        region_seq,
                        grna_start_in_region,
                        grna_end_in_region,
                        pam_start,
                        pam_end,
                    )

                    # Add to table
                    target_short = (
                        target_id if len(target_id) <= 30 else target_id[:27] + "..."
                    )
                    alignment_table.add_row(
                        target_short,
                        str(start + 1),  # 1-based position
                        strand,
                        colored_seq,
                    )

                    shown_targets += 1
                    break

            # Stop if we've shown enough targets
            if shown_targets >= actual_targets_to_show:
                break

        console.print(alignment_table)
        console.print()

    # Add legend
    legend_text = Text()
    legend_text.append("Legend: ", style="bold")
    legend_text.append("█", style="bold red on dark_red")
    legend_text.append(" gRNA  ", style="dim")
    legend_text.append("█", style="bold cyan on dark_cyan")
    legend_text.append(" PAM  ", style="dim")
    legend_text.append("█", style="dim white")
    legend_text.append(" Flanking", style="dim")

    console.print(Panel(legend_text, border_style="dim"))


@app.command()
def compare(
    grna_csv: Path = typer.Option(
        "results/results/validation_run/guide_rnas.csv",
        "--csv",
        "-c",
        help="Path to guide RNA CSV file",
    ),
    grna_ids: str = typer.Option(
        "gRNA_0001,gRNA_0002", "--ids", "-i", help="Comma-separated gRNA IDs to compare"
    ),
):
    """Compare multiple gRNAs side by side."""

    # Load data
    df = pd.read_csv(grna_csv)

    # Parse gRNA IDs
    ids_to_compare = [id.strip() for id in grna_ids.split(",")]

    # Create comparison table
    compare_table = Table(
        title="[bold magenta]gRNA Comparison[/bold magenta]",
        show_header=True,
        header_style="bold cyan",
        box=box.DOUBLE_EDGE,
    )

    compare_table.add_column("Property", style="yellow")
    for grna_id in ids_to_compare:
        compare_table.add_column(grna_id, style="white")

    # Add comparison rows
    properties = ["Sequence", "PAM", "Conservation", "Target Count", "Quality Score"]

    for prop in properties:
        row_data = [prop]
        for grna_id in ids_to_compare:
            grna_row = (
                df[df["Guide_ID"] == grna_id].iloc[0]
                if grna_id in df["Guide_ID"].values
                else None
            )
            if grna_row is not None:
                if prop == "Sequence":
                    row_data.append(f"[red]{grna_row['Sequence']}[/red]")
                elif prop == "PAM":
                    row_data.append(f"[cyan]{grna_row['PAM']}[/cyan]")
                elif prop == "Conservation":
                    row_data.append(f"[green]{grna_row['Conservation']}[/green]")
                elif prop == "Target Count":
                    row_data.append(str(grna_row["Target_Count"]))
                elif prop == "Quality Score":
                    row_data.append(f"{grna_row['Quality_Score']:.3f}")
            else:
                row_data.append("N/A")
        compare_table.add_row(*row_data)

    console.print(compare_table)


@app.command()
def stats(
    grna_csv: Path = typer.Option(
        "results/results/validation_run/guide_rnas.csv",
        "--csv",
        "-c",
        help="Path to guide RNA CSV file",
    ),
):
    """Display statistics about the gRNA results."""

    df = pd.read_csv(grna_csv)

    # Create stats panels
    stats_table = Table(
        title="[bold cyan]gRNA Statistics[/bold cyan]",
        show_header=False,
        box=box.ROUNDED,
    )
    stats_table.add_column("Metric", style="yellow")
    stats_table.add_column("Value", style="green")

    # Calculate statistics
    stats_table.add_row("Total gRNAs", str(len(df)))
    stats_table.add_row(
        "Avg Conservation",
        f"{df['Conservation'].str.rstrip('%').astype(float).mean():.1f}%",
    )
    stats_table.add_row(
        "Max Conservation",
        f"{df['Conservation'].str.rstrip('%').astype(float).max():.1f}%",
    )
    stats_table.add_row(
        "Min Conservation",
        f"{df['Conservation'].str.rstrip('%').astype(float).min():.1f}%",
    )

    # Count PAM types
    pam_counts = df["PAM"].value_counts()
    for pam, count in pam_counts.items():
        stats_table.add_row(f"PAM {pam}", f"{count} ({count/len(df)*100:.1f}%)")

    # Conservation distribution
    cons_values = df["Conservation"].str.rstrip("%").astype(float)
    stats_table.add_row("≥70% Conservation", str(len(df[cons_values >= 70])))
    stats_table.add_row("≥60% Conservation", str(len(df[cons_values >= 60])))
    stats_table.add_row("≥50% Conservation", str(len(df[cons_values >= 50])))

    console.print(stats_table)

    # Create conservation histogram
    console.print("\n[bold]Conservation Distribution:[/bold]")

    # Create bins
    bins = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    hist_data = pd.cut(cons_values, bins=bins).value_counts().sort_index()

    for interval, count in hist_data.items():
        bar_length = int(count / len(df) * 50)
        bar = "█" * bar_length
        console.print(f"{interval}: [green]{bar}[/green] {count}")


if __name__ == "__main__":
    app()
