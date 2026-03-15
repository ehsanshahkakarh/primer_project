#!/usr/bin/env python3
"""
PCR primer pair design pipeline.
Finds conserved regions and designs forward/reverse primer pairs for amplification.
"""

import logging
from pathlib import Path
from typing import Optional, List
import typer
from Bio import SeqIO

from .conserved_finder import KmerBasedFinder, MSAConservedRegionScanner
from .primer_candidates import PrimerCandidateGenerator
from .primer_pairing import PrimerPairMatcher
from .primer_validation import PrimerValidator
from .primer_scoring import PrimerPairScorer
from .primer_output import write_all_outputs
from .output_utils import resolve_output_dir

app = typer.Typer()

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


@app.command()
def main(
    # Input
    fasta_input: Path = typer.Option(
        ..., "--fasta-input", help="Input FASTA file or directory"
    ),
    output_dir: str = typer.Option(
        ..., "-o", "--output-dir", help="Output directory for results"
    ),
    # Amplicon constraints
    amplicon_min: int = typer.Option(
        100, "--amplicon-min", help="Minimum amplicon size (bp)"
    ),
    amplicon_max: int = typer.Option(
        2000, "--amplicon-max", help="Maximum amplicon size (bp), supports long-read"
    ),
    # Primer constraints
    primer_min_length: int = typer.Option(
        18, "--min-length", help="Minimum primer length (bp)"
    ),
    primer_max_length: int = typer.Option(
        25, "--max-length", help="Maximum primer length (bp)"
    ),
    # Thermodynamic constraints
    tm_min: float = typer.Option(55.0, "--tm-min", help="Minimum melting temp (°C)"),
    tm_max: float = typer.Option(65.0, "--tm-max", help="Maximum melting temp (°C)"),
    tm_diff_max: float = typer.Option(
        5.0, "--tm-diff-max", help="Maximum Tm difference between pairs (°C)"
    ),
    # GC constraints
    gc_min: float = typer.Option(40.0, "--gc-min", help="Minimum GC content (%)"),
    gc_max: float = typer.Option(60.0, "--gc-max", help="Maximum GC content (%)"),
    # Conservation
    conservation: float = typer.Option(
        0.8, "--conservation", help="Minimum conservation (0-1)"
    ),
    # K-mer settings
    kmer_size: int = typer.Option(
        20, "--kmer-size", help="K-mer size for conserved region finding"
    ),
    # Method selection
    method: str = typer.Option(
        "msa",
        "--method",
        help="Conserved region finding method: 'msa' (MAFFT L-INS-i, default) or 'kmer' (fast, no alignment)",
    ),
    structure_aware: bool = typer.Option(
        False,
        "--structure-aware",
        help="Use Q-INS-i structure-aware alignment for rRNA (16S/18S). Requires MAFFT extensions.",
    ),
    mafft_auto: bool = typer.Option(
        False,
        "--mafft-auto",
        help="Use MAFFT --auto mode (faster but less accurate than L-INS-i). Good for >100 sequences.",
    ),
    # Deprecated, kept for backward compatibility
    align_regions: bool = typer.Option(
        False,
        "--align-regions",
        hidden=True,
        help="Deprecated: use --method msa instead",
    ),
    max_degenerate_positions: int = typer.Option(
        2,
        "--max-degenerate-positions",
        help="Maximum IUPAC degenerate positions allowed per primer (0-5, default: 2)",
    ),
    # Performance
    threads: int = typer.Option(4, "--threads", help="Number of threads"),
):
    """
    Design PCR primer pairs for amplifying conserved regions.

    This pipeline:
    1. Finds conserved regions using k-mer analysis
    2. Generates forward and reverse primer candidates
    3. Pairs primers within amplicon size constraints
    4. Validates primers against PCR best practices
    5. Scores and ranks primer pairs by quality

    Examples:
        # MSA mode (default, MAFFT L-INS-i alignment)
        ppdesign primer main --fasta-input sequences.fna \\
                            --output-dir primer_results \\
                            --conservation 0.8 \\
                            --amplicon-min 200 \\
                            --amplicon-max 800

        # Structure-aware mode for rRNA (Q-INS-i)
        ppdesign primer main --fasta-input 16s_sequences.fna \\
                            --output-dir primer_results \\
                            --structure-aware

        # Fast k-mer mode (no alignment)
        ppdesign primer main --fasta-input sequences.fna \\
                            --output-dir primer_results \\
                            --method kmer
    """
    logger = logging.getLogger(__name__)

    # Setup output directory
    output_path = resolve_output_dir(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    typer.echo("=" * 60)
    typer.echo("PCR Primer Pair Design Pipeline")
    typer.echo("=" * 60)

    # Validate parameter compatibility
    min_required_region_length = amplicon_min + primer_max_length
    typer.echo(f"\nParameter Validation:")
    typer.echo(f"  Amplicon range: {amplicon_min}-{amplicon_max} bp")
    typer.echo(f"  Primer length: {primer_min_length}-{primer_max_length} bp")
    typer.echo(f"  Min conserved region needed: {min_required_region_length} bp")

    # Load sequences
    typer.echo(f"\n[1/7] Loading sequences from {fasta_input}...")
    sequences = _load_sequences(fasta_input)
    typer.echo(f"  ✓ Loaded {len(sequences)} sequences")

    # Convert to format expected by conserved_finder
    seq_dict = {seq_id: str(record.seq) for seq_id, record in sequences.items()}

    # Find conserved regions
    use_msa = method == "msa" or align_regions

    if use_msa:
        if structure_aware:
            strategy = "qinsi"
        elif mafft_auto:
            strategy = "auto"
        else:
            strategy = "linsi"
        typer.echo(
            f"\n[2/7] Finding conserved regions (MAFFT {strategy}"
            f"{', structure-aware' if structure_aware else ''})..."
        )
        scanner = MSAConservedRegionScanner(
            min_conservation=conservation, mafft_strategy=strategy
        )
        conserved_regions = scanner.find_conserved_regions(
            seq_dict, min_length=primer_min_length, max_length=primer_max_length
        )
    else:
        typer.echo(f"\n[2/7] Finding conserved regions (k-mer, k={kmer_size})...")
        finder = KmerBasedFinder(kmer_size=kmer_size, min_conservation=conservation)
        min_region_length = amplicon_min + primer_max_length
        max_region_length = amplicon_max + primer_max_length * 2
        conserved_regions = finder.find_conserved_regions(
            seq_dict,
            min_length=min_region_length,
            max_length=max_region_length,
        )

    typer.echo(f"  ✓ Found {len(conserved_regions)} conserved regions")

    if not conserved_regions:
        typer.echo(
            "\n⚠ No conserved regions found. Try lowering --conservation threshold.",
            err=True,
        )
        raise typer.Exit(1)

    # Generate primer candidates
    typer.echo("\n[3/7] Generating primer candidates...")
    generator = PrimerCandidateGenerator(
        primer_min_length=primer_min_length,
        primer_max_length=primer_max_length,
        tm_min=tm_min,
        tm_max=tm_max,
        gc_min=gc_min,
        gc_max=gc_max,
        conservation_threshold=conservation,
        max_degenerate_positions=max_degenerate_positions,
    )

    forward_primers = generator.generate_forward_primers(conserved_regions)
    reverse_primers = generator.generate_reverse_primers(conserved_regions)
    typer.echo(
        f"  ✓ Generated {len(forward_primers)} forward and {len(reverse_primers)} reverse primers"
    )

    if not forward_primers or not reverse_primers:
        typer.echo(
            "\n⚠ No valid primers found. Try relaxing Tm/GC constraints.", err=True
        )
        raise typer.Exit(1)

    # Pair primers
    typer.echo("\n[4/7] Pairing primers...")
    matcher = PrimerPairMatcher(
        amplicon_min_size=amplicon_min,
        amplicon_max_size=amplicon_max,
        tm_difference_max=tm_diff_max,
    )

    pairs = matcher.pair_primers(forward_primers, reverse_primers)
    typer.echo(f"  ✓ Created {len(pairs)} primer pairs")

    if not pairs:
        typer.echo(
            "\n⚠ No valid primer pairs found. Try adjusting amplicon size range.",
            err=True,
        )
        raise typer.Exit(1)

    # Validate pairs
    typer.echo("\n[5/7] Validating primer pairs...")
    validator = PrimerValidator()
    valid_pairs = []

    for pair in pairs:
        is_valid, issues = validator.validate_pair(pair)
        if is_valid:
            valid_pairs.append(pair)

    typer.echo(f"  ✓ {len(valid_pairs)} pairs passed validation")

    if not valid_pairs:
        typer.echo(
            "\n⚠ No pairs passed validation. Check primer design parameters.",
            err=True,
        )
        raise typer.Exit(1)

    # Score and rank pairs
    typer.echo("\n[6/7] Scoring and ranking primer pairs...")
    scorer = PrimerPairScorer(tm_difference_max=tm_diff_max)
    ranked_pairs = scorer.score_pairs(valid_pairs)
    typer.echo(f"  ✓ Ranked {len(ranked_pairs)} primer pairs by quality")

    # Write outputs
    typer.echo("\n[7/7] Writing output files...")
    write_all_outputs(ranked_pairs, output_path)
    typer.echo(f"  ✓ primer_pairs.csv")
    typer.echo(f"  ✓ primer_pairs.fasta")
    typer.echo(f"  ✓ amplicon_predictions.tsv")
    typer.echo(f"  ✓ summary.txt")

    # Summary
    typer.echo("\n" + "=" * 60)
    typer.echo("Pipeline Complete!")
    typer.echo("=" * 60)
    typer.echo(f"\nResults saved to: {output_path}")

    if ranked_pairs:
        top_pair = ranked_pairs[0]
        typer.echo(f"\nTop primer pair:")
        typer.echo(f"  Quality score: {top_pair.quality_score:.1f}/100")
        typer.echo(f"  Amplicon size: {top_pair.amplicon_size} bp")
        typer.echo(f"  Conservation: {top_pair.avg_conservation:.1%}")
        typer.echo(f"  Tm difference: {top_pair.tm_difference:.1f}°C")
        typer.echo(f"\n  Forward: {top_pair.forward.sequence}")
        typer.echo(f"  Reverse: {top_pair.reverse.sequence}")


def _load_sequences(fasta_input: Path) -> dict:
    """
    Load sequences from FASTA file or directory.

    Args:
        fasta_input: Path to FASTA file or directory

    Returns:
        Dictionary mapping sequence ID to SeqRecord

    Raises:
        ValueError: If path not found or no sequences loaded
    """
    if not fasta_input.exists():
        raise ValueError(f"FASTA path not found: {fasta_input}")

    sequences = {}

    if fasta_input.is_file():
        # Single file
        for record in SeqIO.parse(fasta_input, "fasta"):
            sequences[record.id] = record
    elif fasta_input.is_dir():
        # Directory of FASTA files
        fasta_files = list(fasta_input.glob("*.fna")) + list(fasta_input.glob("*.fa"))
        fasta_files += list(fasta_input.glob("*.fasta"))

        if not fasta_files:
            raise ValueError(f"No FASTA files found in: {fasta_input}")

        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = f"{fasta_file.stem}_{record.id}"
                sequences[seq_id] = record

    if not sequences:
        raise ValueError(f"No sequences loaded from: {fasta_input}")

    return sequences


if __name__ == "__main__":
    app()
