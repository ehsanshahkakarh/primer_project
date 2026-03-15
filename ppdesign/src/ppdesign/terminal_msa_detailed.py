#!/usr/bin/env python3
"""Detailed MSA viewer showing actual gRNA sequence alignment."""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from rich.console import Console
from rich.text import Text
from rich.panel import Panel
import typer

app = typer.Typer(help="Detailed MSA viewer for gRNA alignments")
console = Console()


def find_and_extract_grna(grna_seq: str, target_seq: str, window: int = 10):
    """Find gRNA in sequence and extract aligned regions."""
    results = []
    target_seq = target_seq.upper()
    grna_seq = grna_seq.upper()

    # Check forward strand
    pos = target_seq.find(grna_seq)
    if pos != -1:
        # Extract with context
        start = max(0, pos - window)
        end = min(len(target_seq), pos + len(grna_seq) + 3 + window)

        before = target_seq[start:pos]
        grna_match = target_seq[pos : pos + len(grna_seq)]
        pam = (
            target_seq[pos + len(grna_seq) : pos + len(grna_seq) + 3]
            if pos + len(grna_seq) + 3 <= len(target_seq)
            else ""
        )
        after = (
            target_seq[pos + len(grna_seq) + 3 : end]
            if pos + len(grna_seq) + 3 <= len(target_seq)
            else ""
        )

        results.append(
            {
                "strand": "+",
                "position": pos + 1,
                "before": before,
                "grna": grna_match,
                "pam": pam,
                "after": after,
                "full": before + grna_match + pam + after,
            }
        )

    # Check reverse strand
    rev_grna = str(Seq(grna_seq).reverse_complement())
    pos = target_seq.find(rev_grna)
    if pos != -1:
        # For reverse strand, we'll show it reverse complemented to align with gRNA
        start = max(0, pos - 3 - window)
        end = min(len(target_seq), pos + len(grna_seq) + window)

        region = target_seq[start:end]
        region_rc = str(Seq(region).reverse_complement())

        # Find gRNA in reverse complemented region
        grna_pos = region_rc.find(grna_seq)
        if grna_pos != -1:
            before = region_rc[:grna_pos]
            grna_match = region_rc[grna_pos : grna_pos + len(grna_seq)]
            pam = region_rc[grna_pos + len(grna_seq) : grna_pos + len(grna_seq) + 3]
            after = region_rc[grna_pos + len(grna_seq) + 3 :]

            results.append(
                {
                    "strand": "-",
                    "position": pos + 1,
                    "before": before[-window:] if len(before) > window else before,
                    "grna": grna_match,
                    "pam": pam,
                    "after": after[:window],
                    "full": before + grna_match + pam + after,
                }
            )

    return results


@app.command()
def align(
    grna_csv: Path = typer.Option(
        "results/results/validation_run/guide_rnas.csv",
        "--csv",
        "-c",
        help="Path to guide RNA CSV file",
    ),
    fasta_file: Path = typer.Option(
        "tests/test.fna", "--fasta", "-f", help="Path to FASTA file with sequences"
    ),
    grna_id: str = typer.Option(
        "gRNA_0001", "--id", "-i", help="Specific gRNA ID to analyze"
    ),
    max_targets: int = typer.Option(
        10, "--max", "-m", help="Maximum number of targets to show"
    ),
):
    """Show detailed alignment of a specific gRNA across all targets."""

    # Load sequences
    console.print("[bold cyan]Loading sequences...[/bold cyan]")
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    # Load gRNA data
    df = pd.read_csv(grna_csv)
    grna_row = df[df["Guide_ID"] == grna_id]

    if grna_row.empty:
        console.print(f"[red]gRNA {grna_id} not found![/red]")
        return

    grna_row = grna_row.iloc[0]
    grna_seq = grna_row["Sequence"]

    # Try to get complete target list from target_mapping.tsv
    target_mapping_file = grna_csv.parent / "target_mapping.tsv"
    if target_mapping_file.exists():
        target_df = pd.read_csv(target_mapping_file, sep="\t")
        grna_targets = target_df[target_df["Guide_ID"] == grna_id][
            "Target_Sequence"
        ].unique()
        targets = list(grna_targets)
    else:
        targets = grna_row["Target_Sequences"].split(";")

    # Create header
    console.print(
        Panel.fit(
            f"[bold magenta]Detailed MSA for {grna_id}[/bold magenta]\n"
            + f"[yellow]gRNA Sequence: {grna_seq}[/yellow]\n"
            + f"[green]Conservation: {grna_row['Conservation']}[/green]",
            border_style="bright_blue",
        )
    )

    # Create consensus line
    console.print("\n[bold]Aligned sequences (all in forward orientation):[/bold]\n")

    # Collect all alignments
    alignments = []
    for target_id in targets[:max_targets]:
        if target_id in sequences:
            matches = find_and_extract_grna(grna_seq, sequences[target_id])
            for match in matches:
                alignments.append(
                    {
                        "target": (
                            target_id[:25] + "..." if len(target_id) > 25 else target_id
                        ),
                        **match,
                    }
                )

    # Display as aligned text
    if alignments:
        # Find max lengths for padding
        max_before = max(len(a["before"]) for a in alignments)
        max_after = max(len(a["after"]) for a in alignments)

        # Print reference gRNA
        console.print("[bold cyan]Reference:[/bold cyan]")
        console.print(
            f"            {'.' * max_before}[bold red]{grna_seq}[/bold red][bold cyan]NGG[/bold cyan]{'.' * max_after}"
        )
        console.print()

        # Print each alignment
        for i, align in enumerate(alignments, 1):
            # Pad sequences for alignment
            before_padded = align["before"].rjust(max_before)
            after_padded = align["after"].ljust(max_after)

            # Create colored text
            line = Text()
            line.append(f"{i:2}. {align['strand']} ", style="dim")
            line.append(before_padded, style="dim white")
            line.append(align["grna"], style="bold red on dark_red")
            line.append(align["pam"], style="bold cyan on dark_cyan")
            line.append(after_padded, style="dim white")
            line.append(f"  {align['target']}", style="dim yellow")

            console.print(line)

        # Show conservation
        console.print(f"\n[green]✓ {len(alignments)} sequences aligned[/green]")

        # Check if all gRNAs match exactly
        all_match = all(a["grna"] == grna_seq for a in alignments)
        if all_match:
            console.print(
                "[bold green]✓ All gRNA sequences match perfectly![/bold green]"
            )
        else:
            mismatches = [a for a in alignments if a["grna"] != grna_seq]
            console.print(
                f"[yellow]⚠ {len(mismatches)} sequences have mismatches[/yellow]"
            )


@app.command()
def consensus(
    grna_csv: Path = typer.Option(
        "results/results/validation_run/guide_rnas.csv",
        "--csv",
        "-c",
        help="Path to guide RNA CSV file",
    ),
    fasta_file: Path = typer.Option(
        "tests/test.fna", "--fasta", "-f", help="Path to FASTA file"
    ),
    num_grnas: int = typer.Option(
        3, "--num", "-n", help="Number of top gRNAs to analyze"
    ),
):
    """Show consensus view across multiple gRNAs."""

    # Load data
    sequences = {
        record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")
    }
    df = pd.read_csv(grna_csv)

    console.print(
        Panel.fit(
            "[bold magenta]Consensus MSA View[/bold magenta]",
            border_style="bright_blue",
        )
    )

    for idx in range(min(num_grnas, len(df))):
        row = df.iloc[idx]
        grna_id = row["Guide_ID"]
        grna_seq = row["Sequence"]
        targets = row["Target_Sequences"].split(";")

        # Count perfect matches
        perfect_matches = 0
        total_targets = 0

        for target_id in targets[:5]:  # Check first 5 targets
            if target_id in sequences:
                matches = find_and_extract_grna(grna_seq, sequences[target_id])
                if matches and matches[0]["grna"] == grna_seq:
                    perfect_matches += 1
                total_targets += 1

        # Display summary
        console.print(f"\n[bold yellow]{grna_id}:[/bold yellow]")
        console.print(f"  Sequence: [red]{grna_seq}[/red]")
        console.print(
            f"  Perfect matches: [green]{perfect_matches}/{total_targets}[/green]"
        )
        console.print(f"  Conservation: {row['Conservation']}")


if __name__ == "__main__":
    app()
