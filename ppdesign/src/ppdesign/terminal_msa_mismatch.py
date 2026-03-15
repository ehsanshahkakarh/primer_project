#!/usr/bin/env python3
"""MSA viewer with mismatch visualization for gRNA alignments."""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from rich.console import Console
from rich.text import Text
from rich.panel import Panel
from rich.table import Table
from rich import box
import typer
from typing import List, Dict, Optional

app = typer.Typer(help="MSA viewer with mismatch visualization")
console = Console()


def find_best_match(
    grna_seq: str, target_seq: str, max_mismatches: int = 3
) -> Optional[Dict]:
    """Find best gRNA match in sequence allowing mismatches."""
    target_seq = target_seq.upper()
    grna_seq = grna_seq.upper()
    best_match = None
    min_mismatches = max_mismatches + 1

    # Check forward strand
    for i in range(len(target_seq) - len(grna_seq) + 1):
        segment = target_seq[i : i + len(grna_seq)]
        mismatches = sum(1 for a, b in zip(grna_seq, segment) if a != b)

        if mismatches <= max_mismatches:
            # Check for PAM
            if i + len(grna_seq) + 2 <= len(target_seq):
                pam = target_seq[i + len(grna_seq) : i + len(grna_seq) + 3]
                if len(pam) >= 2 and pam[1:3] in ["GG", "AG"]:
                    if mismatches < min_mismatches:
                        min_mismatches = mismatches
                        best_match = {
                            "position": i,
                            "strand": "+",
                            "sequence": segment,
                            "pam": pam,
                            "mismatches": mismatches,
                            "mismatch_positions": [
                                j
                                for j, (a, b) in enumerate(zip(grna_seq, segment))
                                if a != b
                            ],
                        }

    # Check reverse strand
    rev_comp = str(Seq(target_seq).reverse_complement())
    for i in range(len(rev_comp) - len(grna_seq) + 1):
        segment = rev_comp[i : i + len(grna_seq)]
        mismatches = sum(1 for a, b in zip(grna_seq, segment) if a != b)

        if mismatches <= max_mismatches:
            # Check for PAM
            if i + len(grna_seq) + 2 <= len(rev_comp):
                pam = rev_comp[i + len(grna_seq) : i + len(grna_seq) + 3]
                if len(pam) >= 2 and pam[1:3] in ["GG", "AG"]:
                    if mismatches < min_mismatches:
                        min_mismatches = mismatches
                        # Convert position back to original sequence coordinates
                        orig_pos = len(target_seq) - i - len(grna_seq)
                        best_match = {
                            "position": orig_pos,
                            "strand": "-",
                            "sequence": segment,
                            "pam": pam,
                            "mismatches": mismatches,
                            "mismatch_positions": [
                                j
                                for j, (a, b) in enumerate(zip(grna_seq, segment))
                                if a != b
                            ],
                        }

    return best_match


def create_colored_alignment(
    grna_seq: str,
    match_seq: str,
    mismatch_positions: List[int],
    before: str = "",
    after: str = "",
    pam: str = "",
) -> Text:
    """Create colored alignment showing matches and mismatches."""
    colored = Text()

    # Add context before
    if before:
        colored.append(before[-10:], style="dim white")
        colored.append(" ", style="dim")

    # Add gRNA region with match/mismatch coloring
    for i, (ref, obs) in enumerate(zip(grna_seq, match_seq)):
        if i in mismatch_positions:
            # Mismatch - red background
            colored.append(obs, style="bold white on red")
        else:
            # Match - green background
            colored.append(obs, style="bold white on green")

    # Add PAM
    if pam:
        colored.append(" ", style="dim")
        colored.append(pam, style="bold cyan on dark_cyan")

    # Add context after
    if after:
        colored.append(" ", style="dim")
        colored.append(after[:10], style="dim white")

    return colored


@app.command()
def show_all(
    grna_csv: Path = typer.Option(
        "results/results/validation_run/guide_rnas.csv",
        "--csv",
        "-c",
        help="Path to guide RNA CSV file",
    ),
    fasta_file: Path = typer.Option(
        "tests/test.fna", "--fasta", "-f", help="Path to FASTA file"
    ),
    grna_id: str = typer.Option("gRNA_0001", "--id", "-i", help="gRNA ID to analyze"),
    max_mismatches: int = typer.Option(
        2, "--mismatches", "-m", help="Maximum mismatches to show"
    ),
):
    """Show all sequences including those with mismatches."""

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

    # Get target list from target_mapping if available
    target_mapping_file = grna_csv.parent / "target_mapping.tsv"
    perfect_match_ids = set()
    if target_mapping_file.exists():
        target_df = pd.read_csv(target_mapping_file, sep="\t")
        perfect_match_ids = set(
            target_df[target_df["Guide_ID"] == grna_id]["Target_Sequence"].unique()
        )

    # Create header
    console.print(
        Panel.fit(
            "[bold magenta]Complete MSA with Mismatch Visualization[/bold magenta]\n"
            + f"[yellow]gRNA: {grna_seq}[/yellow]\n"
            + f"[green]Perfect matches: {len(perfect_match_ids)}/26[/green]\n"
            + f"[red]Showing mismatches up to {max_mismatches} differences[/red]",
            border_style="bright_blue",
        )
    )

    # Create results table
    results_table = Table(
        title="[bold]All Sequences Alignment[/bold]",
        show_header=True,
        header_style="bold cyan",
        box=box.SIMPLE_HEAD,
        padding=(0, 1),
    )

    results_table.add_column("#", style="dim", width=3)
    results_table.add_column(
        "Sequence ID", style="yellow", width=30, overflow="ellipsis"
    )
    results_table.add_column("Pos", style="green", justify="right", width=6)
    results_table.add_column("Str", style="cyan", justify="center", width=3)
    results_table.add_column("MM", style="red", justify="right", width=3)
    results_table.add_column("Alignment", no_wrap=False)

    # Process all sequences
    perfect_matches = []
    near_matches = []
    no_matches = []

    for seq_id, seq in sequences.items():
        match = find_best_match(grna_seq, seq, max_mismatches)

        if match:
            if match["mismatches"] == 0:
                perfect_matches.append((seq_id, match))
            else:
                near_matches.append((seq_id, match))
        else:
            no_matches.append(seq_id)

    # Sort by number of mismatches
    perfect_matches.sort(key=lambda x: x[0])
    near_matches.sort(key=lambda x: (x[1]["mismatches"], x[0]))

    row_num = 1

    # Show perfect matches first
    if perfect_matches:
        results_table.add_row(
            "", "[bold green]PERFECT MATCHES[/bold green]", "", "", "", ""
        )
        for seq_id, match in perfect_matches:
            # Extract context
            start = max(0, match["position"] - 10)
            end = min(len(sequences[seq_id]), match["position"] + len(grna_seq) + 13)
            before = sequences[seq_id][start : match["position"]]
            after = sequences[seq_id][match["position"] + len(grna_seq) + 3 : end]

            colored_seq = create_colored_alignment(
                grna_seq,
                match["sequence"],
                match["mismatch_positions"],
                before,
                after,
                match["pam"],
            )

            seq_id_short = seq_id if len(seq_id) <= 30 else seq_id[:27] + "..."
            results_table.add_row(
                str(row_num),
                seq_id_short,
                str(match["position"] + 1),
                match["strand"],
                "0",
                colored_seq,
            )
            row_num += 1

    # Show near matches
    if near_matches:
        results_table.add_row(
            "", "[bold yellow]NEAR MATCHES[/bold yellow]", "", "", "", ""
        )
        for seq_id, match in near_matches:
            # Extract context
            start = max(0, match["position"] - 10)
            end = min(len(sequences[seq_id]), match["position"] + len(grna_seq) + 13)
            before = (
                sequences[seq_id][start : match["position"]]
                if match["strand"] == "+"
                else ""
            )
            after = (
                sequences[seq_id][match["position"] + len(grna_seq) + 3 : end]
                if match["strand"] == "+"
                else ""
            )

            # For reverse strand, need to extract from reverse complement
            if match["strand"] == "-":
                rev_seq = str(Seq(sequences[seq_id]).reverse_complement())
                rev_pos = len(sequences[seq_id]) - match["position"] - len(grna_seq)
                start = max(0, rev_pos - 10)
                end = min(len(rev_seq), rev_pos + len(grna_seq) + 13)
                before = rev_seq[start:rev_pos]
                after = rev_seq[rev_pos + len(grna_seq) + 3 : end]

            colored_seq = create_colored_alignment(
                grna_seq,
                match["sequence"],
                match["mismatch_positions"],
                before,
                after,
                match["pam"],
            )

            seq_id_short = seq_id if len(seq_id) <= 30 else seq_id[:27] + "..."
            results_table.add_row(
                str(row_num),
                seq_id_short,
                str(match["position"] + 1),
                match["strand"],
                str(match["mismatches"]),
                colored_seq,
            )
            row_num += 1

    # Show sequences with no matches
    if no_matches:
        results_table.add_row("", "[bold red]NO MATCHES[/bold red]", "", "", "", "")
        for seq_id in no_matches:
            seq_id_short = seq_id if len(seq_id) <= 30 else seq_id[:27] + "..."
            results_table.add_row(
                str(row_num),
                seq_id_short,
                "-",
                "-",
                "-",
                Text("No match found", style="dim red"),
            )
            row_num += 1

    console.print(results_table)

    # Summary
    console.print("\n[bold]Summary:[/bold]")
    console.print(f"  [green]Perfect matches: {len(perfect_matches)}[/green]")
    console.print(
        f"  [yellow]Near matches (1-{max_mismatches} mismatches): {len(near_matches)}[/yellow]"
    )
    console.print(f"  [red]No matches: {len(no_matches)}[/red]")
    console.print(f"  [cyan]Total sequences: {len(sequences)}[/cyan]")


if __name__ == "__main__":
    app()
