#!/usr/bin/env python3
"""
Terminal-based alignment viewer for specific gRNA sequences.
Shows how a selected gRNA aligns to all target sequences.
"""

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple
import sys
from Bio import SeqIO
from rich.console import Console
from rich.table import Table
from rich.text import Text
from rich.panel import Panel
from rich import box


def reverse_complement(seq: str) -> str:
    """Return reverse complement of sequence."""
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
    }
    return "".join(complement.get(base, base) for base in seq[::-1])


def find_grna_in_sequence(
    grna_seq: str, target_seq: str, pam_type: str = "NGG"
) -> List[Tuple[int, str, str]]:
    """Find all occurrences of gRNA in target sequence."""
    matches = []
    grna_len = len(grna_seq)
    pam_len = len(pam_type)

    # Search forward strand
    for i in range(len(target_seq) - grna_len - pam_len + 1):
        window = target_seq[i : i + grna_len]
        pam = target_seq[i + grna_len : i + grna_len + pam_len]

        # Check if gRNA matches (allow for some mismatches)
        mismatches = sum(1 for a, b in zip(grna_seq, window) if a != b)
        if mismatches <= 2:  # Allow up to 2 mismatches
            if pam_type == "NGG" and pam[1:] == "GG":
                matches.append((i, window, "+", mismatches))
            elif pam_type == "NAG" and pam[1:] == "AG":
                matches.append((i, window, "+", mismatches))

    # Search reverse strand
    rc_grna = reverse_complement(grna_seq)
    for i in range(len(target_seq) - grna_len - pam_len + 1):
        window = target_seq[i : i + grna_len]
        pam = target_seq[i + grna_len : i + grna_len + pam_len]

        mismatches = sum(1 for a, b in zip(rc_grna, window) if a != b)
        if mismatches <= 2:
            if pam_type == "NGG" and pam[1:] == "GG":
                matches.append((i, window, "-", mismatches))
            elif pam_type == "NAG" and pam[1:] == "AG":
                matches.append((i, window, "-", mismatches))

    return matches


def load_grna_data(results_dir: Path) -> Dict:
    """Load gRNA data from results directory."""
    grna_data = {}

    # Load guide_rnas.csv
    csv_file = results_dir / "guide_rnas.csv"
    if not csv_file.exists():
        raise FileNotFoundError(f"guide_rnas.csv not found in {results_dir}")

    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            grna_id = row["Guide_ID"]
            grna_data[grna_id] = {
                "sequence": row["Sequence"],
                "pam": row["PAM"],
                "strand": row["Strand"],
                "conservation": row["Conservation"],
                "targets": (
                    row.get("Target_Sequences", "").split(";")
                    if row.get("Target_Sequences")
                    else []
                ),
            }

    return grna_data


def load_target_sequences(fasta_path: Path) -> Dict[str, str]:
    """Load target sequences from FASTA file."""
    sequences = {}

    if fasta_path.is_file():
        files = [fasta_path]
    else:
        files = (
            list(fasta_path.glob("*.fna"))
            + list(fasta_path.glob("*.fasta"))
            + list(fasta_path.glob("*.fa"))
        )

    for fasta_file in files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = str(record.seq).upper()

    return sequences


def create_alignment_display(
    grna_seq: str,
    target_seq: str,
    position: int,
    strand: str,
    mismatches: int,
    full_sequence: str,
    window: int = 30,
) -> Text:
    """Create colored alignment display showing gRNA in genomic context."""

    # Get the full sequence context
    grna_len = len(grna_seq)

    if strand == "+":
        # Forward strand
        region_start = max(0, position - window)
        region_end = min(
            len(full_sequence), position + grna_len + window + 3
        )  # Include PAM
        region_seq = full_sequence[region_start:region_end]

        # Calculate positions in the extracted region
        grna_start_in_region = position - region_start
        grna_end_in_region = grna_start_in_region + grna_len
        pam_start = grna_end_in_region
        pam_end = min(pam_start + 3, len(region_seq))

    else:
        # Reverse strand - show reverse complement
        from Bio.Seq import Seq

        region_start = max(0, position - window - 3)  # Include PAM
        region_end = min(len(full_sequence), position + grna_len + window)
        region_seq = full_sequence[region_start:region_end]

        # Reverse complement the region so gRNA appears in forward orientation
        region_seq = str(Seq(region_seq).reverse_complement())

        # After reverse complement, positions are flipped
        total_len = len(region_seq)
        grna_start_in_region = total_len - (position + grna_len - region_start)
        grna_end_in_region = grna_start_in_region + grna_len
        pam_start = grna_end_in_region
        pam_end = min(pam_start + 3, len(region_seq))

    # Create colored sequence with highlighted gRNA and PAM
    colored = Text()

    for i, nt in enumerate(region_seq):
        if grna_start_in_region <= i < grna_end_in_region:
            # gRNA region - check for mismatches
            grna_pos = i - grna_start_in_region
            if grna_pos < len(grna_seq):
                if region_seq[i] == grna_seq[grna_pos]:
                    colored.append(nt, style="bold green on dark_green")
                else:
                    colored.append(nt, style="bold red on dark_red")
            else:
                colored.append(nt, style="bold red on dark_red")
        elif pam_start <= i < pam_end:
            # PAM region
            colored.append(nt, style="bold cyan on dark_cyan")
        else:
            # Context sequence
            colored.append(nt, style="dim white")

    # Add position and strand info
    result = Text()
    result.append(f"Pos {position + 1:5d} [{strand}] ", style="dim cyan")
    result.append(colored)

    if mismatches == 0:
        result.append("  ✓", style="bold green")
    else:
        result.append(f"  {mismatches}⨯", style="yellow")

    return result


def view_grna_alignment(
    grna_id: str, results_dir: Path, fasta_path: Path, max_targets: int = 20
):
    """Display alignment of specific gRNA to target sequences."""
    console = Console()

    # Load data
    try:
        grna_data = load_grna_data(results_dir)
        if grna_id not in grna_data:
            console.print(f"[red]Error: {grna_id} not found in results[/red]")
            console.print(f"Available gRNAs: {', '.join(sorted(grna_data.keys()))}")
            return

        target_sequences = load_target_sequences(fasta_path)

    except Exception as e:
        console.print(f"[red]Error loading data: {e}[/red]")
        return

    # Get gRNA info
    grna_info = grna_data[grna_id]
    grna_seq = grna_info["sequence"]

    # Create header panel
    header = Panel(
        f"[bold cyan]Guide RNA: {grna_id}[/bold cyan]\n"
        f"Sequence: [yellow]{grna_seq}[/yellow]\n"
        f"PAM: {grna_info['pam']} | Strand: {grna_info['strand']} | Conservation: {grna_info['conservation']}",
        title="gRNA Information",
        box=box.ROUNDED,
    )
    console.print(header)
    console.print()

    # Find and display matches
    all_matches = []

    for target_id in grna_info["targets"][:max_targets]:
        if target_id not in target_sequences:
            continue

        target_seq = target_sequences[target_id]
        matches = find_grna_in_sequence(grna_seq, target_seq, grna_info["pam"])

        for pos, match_seq, strand, mismatches in matches:
            all_matches.append((target_id, pos, match_seq, strand, mismatches))

    # Sort by number of mismatches (perfect matches first)
    all_matches.sort(key=lambda x: x[4])

    # Create alignment table
    table = Table(
        title=f"Target Sequence Alignments (showing top {min(len(all_matches), max_targets)})",
        show_lines=True,
    )
    table.add_column("Target ID", style="cyan", width=50)
    table.add_column("Alignment", style="white", no_wrap=False)

    for target_id, pos, match_seq, strand, mismatches in all_matches[:max_targets]:
        # Truncate long target IDs for display
        display_id = target_id if len(target_id) <= 50 else target_id[:47] + "..."
        # Get the full sequence for context display
        full_seq = target_sequences.get(target_id, "")
        alignment = create_alignment_display(
            grna_seq, match_seq, pos, strand, mismatches, full_seq
        )
        table.add_row(display_id, alignment)

    console.print(table)

    # Summary statistics
    perfect_matches = sum(1 for m in all_matches if m[4] == 0)
    one_mismatch = sum(1 for m in all_matches if m[4] == 1)
    two_mismatches = sum(1 for m in all_matches if m[4] == 2)

    summary = Panel(
        f"[green]Perfect matches: {perfect_matches}[/green]\n"
        f"[yellow]1 mismatch: {one_mismatch}[/yellow]\n"
        f"[orange1]2 mismatches: {two_mismatches}[/orange1]\n"
        f"Total targets: {len(grna_info['targets'])}",
        title="Match Statistics",
        box=box.ROUNDED,
    )
    console.print()
    console.print(summary)


def main():
    parser = argparse.ArgumentParser(
        description="View gRNA alignments to target sequences"
    )
    parser.add_argument("grna_id", help="gRNA identifier (e.g., gRNA_0008)")
    parser.add_argument(
        "--results-dir",
        "-r",
        required=True,
        help="Results directory from ppdesign-grna",
    )
    parser.add_argument(
        "--fasta-input", "-f", required=True, help="Original FASTA file or directory"
    )
    parser.add_argument(
        "--max-targets", "-m", type=int, default=20, help="Maximum targets to display"
    )

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    fasta_path = Path(args.fasta_input)

    if not results_dir.exists():
        print(f"Error: Results directory {results_dir} not found")
        sys.exit(1)

    if not fasta_path.exists():
        print(f"Error: FASTA path {fasta_path} not found")
        sys.exit(1)

    view_grna_alignment(args.grna_id, results_dir, fasta_path, args.max_targets)


if __name__ == "__main__":
    main()
