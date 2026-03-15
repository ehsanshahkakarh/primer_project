#!/usr/bin/env python3
"""
Command-line interface for CRISPR guide RNA design.
Identifies conserved guide RNAs with SpCas9 PAM sites across multiple sequences.
"""

import logging
from pathlib import Path
from typing import Optional, List
import typer
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from .guide_rna_finder import GuideRNAFinder, GuideRNA
from .output_utils import resolve_output_dir

app = typer.Typer()

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

logger = logging.getLogger(__name__)


class GuideRNAPipeline:
    """Pipeline for guide RNA design with CLI output handling."""

    def __init__(self, output_dir: str):
        """Initialize the pipeline."""
        self.output_dir = resolve_output_dir(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)

    def load_sequences(self, fasta_input: str) -> dict:
        """Load sequences from FASTA file or directory."""
        fasta_path = Path(fasta_input)

        if not fasta_path.exists():
            raise ValueError(f"FASTA path not found: {fasta_input}")

        # Check if input is a file or directory
        if fasta_path.is_file():
            # Single file provided
            fasta_files = [fasta_path]
        else:
            # Directory provided - find all FASTA files
            fasta_files = (
                list(fasta_path.glob("*.fna"))
                + list(fasta_path.glob("*.fasta"))
                + list(fasta_path.glob("*.fa"))
            )

        if not fasta_files:
            raise ValueError(f"No FASTA files found in: {fasta_input}")

        self.logger.info(f"Found {len(fasta_files)} FASTA file(s)")

        # Load sequences
        sequences = {}
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                # Use full record ID to avoid truncation issues
                seq_id = record.id
                if seq_id in sequences:
                    self.logger.warning(f"Duplicate sequence ID found: {seq_id}")
                sequences[seq_id] = str(record.seq).upper()

        self.logger.info(f"Loaded {len(sequences)} sequences")
        return sequences

    def save_results(
        self,
        guides: List[GuideRNA],
        sequences: dict,
        pam_type: str,
        conservation_threshold: float,
        max_degenerate: int,
        perfect_coverage: bool = False,
        coverage_stats: Optional[dict] = None,
    ):
        """Save guide RNA results in multiple formats."""

        # Save as CSV
        csv_file = self.output_dir / "guide_rnas.csv"
        df_data = []

        for i, guide in enumerate(guides, 1):
            row_data = {
                "Guide_ID": f"gRNA_{i:04d}",
                "Sequence": guide.consensus or guide.sequence,
                "PAM": guide.pam,
                "Strand": guide.strand,
                "Conservation": f"{guide.conservation:.2%}",
                "Target_Count": len(guide.positions),
            }

            # Add required target information if applicable
            if guide.required_targets_total > 0:
                row_data["Required_Targets"] = (
                    f"{guide.required_targets_hit}/{guide.required_targets_total}"
                )
                row_data["Required_Coverage"] = (
                    f"{(guide.required_targets_hit/guide.required_targets_total):.1%}"
                )

            row_data.update(
                {
                    "Degenerate_Bases": guide.degenerate_count,
                    "Specificity_Score": f"{guide.specificity_score:.3f}",
                    "Quality_Score": f"{guide.quality_score:.3f}",
                    "Target_Sequences": ";".join(list(guide.positions.keys())[:10]),
                }
            )

            df_data.append(row_data)

        df = pd.DataFrame(df_data)
        df.to_csv(csv_file, index=False)
        self.logger.info(f"Saved CSV results to {csv_file}")

        # Save as FASTA
        fasta_file = self.output_dir / "guide_rnas.fasta"
        fasta_records = []

        for i, guide in enumerate(guides, 1):
            record = SeqRecord(
                Seq(guide.consensus or guide.sequence),
                id=f"gRNA_{i:04d}",
                description=(
                    f"PAM={guide.pam} strand={guide.strand} "
                    f"conservation={guide.conservation:.2%} "
                    f"targets={len(guide.positions)} "
                    f"quality={guide.quality_score:.3f}"
                ),
            )
            fasta_records.append(record)

        SeqIO.write(fasta_records, fasta_file, "fasta")
        self.logger.info(f"Saved FASTA results to {fasta_file}")

        # Save summary statistics
        summary_file = self.output_dir / "summary.txt"
        with open(summary_file, "w") as f:
            f.write("CRISPR Guide RNA Design Summary\n")
            f.write("=" * 60 + "\n\n")

            f.write("Input Parameters:\n")
            f.write(f"  PAM Type: {pam_type}\n")
            f.write(f"  Conservation Threshold: {conservation_threshold:.1%}\n")
            f.write(f"  Max Degenerate Bases: {max_degenerate}\n")
            f.write(f"  Input Sequences: {len(sequences)}\n")
            if perfect_coverage:
                f.write("  Perfect Coverage Mode: Enabled\n")
                if coverage_stats:
                    min_cov = coverage_stats.get('min_coverage_per_target')
                    max_guides = coverage_stats.get('max_total_grnas')
                    if min_cov is not None:
                        f.write(f"  Min Coverage Per Target: {min_cov}\n")
                    if max_guides is not None:
                        f.write(f"  Max Total gRNAs: {max_guides}\n")
            f.write("\n")

            f.write("Results:\n")
            f.write(f"  Total Guide RNAs Found: {len(guides)}\n")

            if guides:
                avg_conservation = sum(g.conservation for g in guides) / len(guides)
                f.write(f"  Average Conservation: {avg_conservation:.1%}\n")

                no_degen = sum(1 for g in guides if g.degenerate_count == 0)
                f.write(
                    f"  Guides with No Degeneracy: {no_degen} ({no_degen/len(guides):.1%})\n"
                )

                ngg_count = sum(1 for g in guides if g.pam == "NGG")
                f.write(f"  Guides with NGG PAM: {ngg_count}\n")

                nag_count = sum(1 for g in guides if g.pam == "NAG")
                f.write(f"  Guides with NAG PAM: {nag_count}\n")

                # Coverage statistics
                all_targets = set()
                for guide in guides:
                    all_targets.update(guide.positions.keys())
                coverage = len(all_targets) / len(sequences)
                f.write(
                    f"  Sequence Coverage: {len(all_targets)}/{len(sequences)} ({coverage:.1%})\n"
                )

                # Perfect coverage specific statistics
                if perfect_coverage and coverage_stats:
                    f.write("\nPerfect Coverage Statistics:\n")
                    f.write(
                        f"  Targets with Sufficient Coverage: {coverage_stats.get('covered_targets', 0)}/{coverage_stats.get('total_targets', 0)} ({coverage_stats.get('coverage_percentage', 0):.1%})\n"
                    )
                    f.write(
                        f"  Average Guides per Target: {coverage_stats.get('average_guides_per_target', 0):.1f}\n"
                    )
                    f.write(
                        f"  Min Guides per Target: {coverage_stats.get('min_guides_per_target', 0)}\n"
                    )
                    f.write(
                        f"  Max Guides per Target: {coverage_stats.get('max_guides_per_target', 0)}\n"
                    )

                f.write("\n")

                # Top 5 guides
                f.write("Top 5 Guide RNAs by Quality Score:\n")
                f.write("-" * 40 + "\n")
                for i, guide in enumerate(guides[:5], 1):
                    f.write(f"\n{i}. {guide.consensus or guide.sequence}\n")
                    f.write(f"   PAM: {guide.pam} ({guide.strand} strand)\n")
                    f.write(
                        f"   Conservation: {guide.conservation:.1%} ({len(guide.positions)} sequences)\n"
                    )
                    f.write(f"   Quality Score: {guide.quality_score:.3f}\n")

                # Most conserved guide
                if guides:
                    best = max(guides, key=lambda x: x.conservation)
                    f.write("\n\nMost Conserved Guide RNA:\n")
                    f.write(f"  Sequence: {best.consensus or best.sequence}\n")
                    f.write(f"  Conservation: {best.conservation:.1%}\n")
                    f.write(
                        f"  Targets: {len(best.positions)}/{len(sequences)} sequences\n"
                    )

        self.logger.info(f"Saved summary to {summary_file}")

        # Save detailed target mapping
        targets_file = self.output_dir / "target_mapping.tsv"
        with open(targets_file, "w") as f:
            f.write("Guide_ID\tGuide_Sequence\tTarget_Sequence\tPosition\tStrand\n")

            for i, guide in enumerate(guides, 1):
                guide_id = f"gRNA_{i:04d}"
                guide_seq = guide.consensus or guide.sequence

                for seq_id, positions in guide.positions.items():
                    for pos, strand in positions:
                        f.write(f"{guide_id}\t{guide_seq}\t{seq_id}\t{pos}\t{strand}\n")

        self.logger.info(f"Saved target mapping to {targets_file}")

        # Save coverage matrix for perfect coverage mode
        if (
            perfect_coverage
            and coverage_stats
            and coverage_stats.get("coverage_matrix")
        ):
            coverage_file = self.output_dir / "coverage_matrix.tsv"
            with open(coverage_file, "w") as f:
                # Header: Guide_ID followed by all sequence IDs
                seq_ids = list(sequences.keys())
                f.write("Guide_ID\t" + "\t".join(seq_ids) + "\n")

                # Data rows: each guide with 1/0 coverage for each sequence
                matrix = coverage_stats["coverage_matrix"]
                for i, row in enumerate(matrix):
                    guide_id = f"gRNA_{i+1:04d}"
                    f.write(guide_id + "\t" + "\t".join(map(str, row)) + "\n")

            self.logger.info(f"Saved coverage matrix to {coverage_file}")

            # Save per-target coverage details
            target_coverage_file = self.output_dir / "target_coverage_details.tsv"
            with open(target_coverage_file, "w") as f:
                f.write(
                    "Target_ID\tGuide_Count\tSufficient_Coverage\tCovering_Guides\n"
                )

                target_stats = coverage_stats.get("target_stats", {})
                for seq_id, stats in target_stats.items():
                    guide_indices = [
                        f"gRNA_{i+1:04d}" for i in stats.get("guide_indices", [])
                    ]
                    f.write(
                        f"{seq_id}\t{stats.get('guide_count', 0)}\t{stats.get('sufficient_coverage', False)}\t{';'.join(guide_indices)}\n"
                    )

            self.logger.info(f"Saved target coverage details to {target_coverage_file}")


@app.command()
def main(
    fasta_input: str = typer.Option(
        None,
        "--fasta-input",
        "-f",
        help="Path to FASTA file or directory containing target sequences",
    ),
    fasta_dir: str = typer.Option(
        None,
        "--fasta-dir",
        help="[Deprecated: use --fasta-input] Directory containing FASTA files",
    ),
    output_dir: str = typer.Option(
        ..., "--output-dir", "-o", help="Output directory for results"
    ),
    pam_type: str = typer.Option(
        "spCas9",
        "--pam-type",
        help="PAM type: spCas9 (NGG/NAG), custom patterns comma-separated",
    ),
    conservation: float = typer.Option(
        0.8,
        "--conservation",
        "-c",
        min=0.0,
        max=1.0,
        help="Minimum conservation threshold (0.0-1.0). In perfect coverage mode, guides must still meet this threshold to be considered.",
    ),
    max_degenerate: int = typer.Option(
        3,
        "--max-degenerate",
        "-d",
        min=0,
        help="Maximum number of degenerate bases allowed",
    ),
    no_degenerate: bool = typer.Option(
        False,
        "--no-degenerate",
        help="Disable degenerate bases - use exact sequences only (prioritizes required targets)",
    ),
    allow_mismatches: int = typer.Option(
        2,
        "--mismatches",
        "-m",
        min=0,
        help="Number of mismatches allowed when clustering guides",
    ),
    threads: int = typer.Option(
        1, "--threads", "-t", help="Number of threads (currently single-threaded)"
    ),
    top_n: Optional[int] = typer.Option(
        None, "--top-n", help="Output only top N guides by quality score"
    ),
    min_targets: Optional[int] = typer.Option(
        None, "--min-targets", help="Minimum number of target sequences required"
    ),
    required_targets: Optional[List[str]] = typer.Option(
        None,
        "--required-targets",
        help="Sequence IDs that must be targeted (space-separated)",
    ),
    required_targets_file: Optional[str] = typer.Option(
        None,
        "--required-targets-file",
        help="File with sequence IDs (one per line) that must be targeted",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose output"
    ),
    min_guides: int = typer.Option(
        50,
        "--min-guides",
        help="Target number of guides to return (returns all available if fewer exist)",
    ),
    max_guides: Optional[int] = typer.Option(
        None,
        "--max-guides",
        help="Maximum number of guides to return (None = no limit)",
    ),
    min_additional_coverage: float = typer.Option(
        0.05,
        "--min-additional-coverage",
        min=0.0,
        max=1.0,
        help="Minimum additional coverage threshold for selecting extra guides (0.0-1.0)",
    ),
    perfect_coverage: bool = typer.Option(
        False,
        "--perfect-coverage",
        help="Enable perfect coverage mode: ensure each target has N perfect gRNA matches",
    ),
    min_coverage_per_target: int = typer.Option(
        3,
        "--min-coverage-per-target",
        min=1,
        help="Minimum number of perfect gRNA matches required per target sequence (perfect coverage mode)",
    ),
    max_total_grnas: int = typer.Option(
        20,
        "--max-total-grnas",
        min=1,
        help="Maximum total number of gRNAs to select (perfect coverage mode)",
    ),
    min_gc: Optional[float] = typer.Option(
        None,
        "--min-gc",
        min=0.0,
        max=100.0,
        help="Minimum GC content percentage for guide RNAs (0-100)",
    ),
    max_gc: Optional[float] = typer.Option(
        None,
        "--max-gc",
        min=0.0,
        max=100.0,
        help="Maximum GC content percentage for guide RNAs (0-100)",
    ),
    cluster_threshold: Optional[float] = typer.Option(
        None,
        "--cluster-threshold",
        min=0.0,
        max=1.0,
        help="Similarity threshold for clustering sequences (0.0-1.0)",
    ),
):
    """
    Design CRISPR guide RNAs targeting conserved regions across multiple sequences.

    This tool identifies guide RNAs with SpCas9-compatible PAM sites (NGG/NAG)
    that are conserved across input sequences, allowing for targeted editing
    of multiple related sequences with a single guide.

    Example:
        ppdesign-grna -f sequences/ -o results/guides/ -c 0.8 --max-degenerate 2
    """

    try:
        # Set logging level
        if verbose:
            logging.getLogger().setLevel(logging.DEBUG)

        # Handle backward compatibility for fasta_dir parameter
        if fasta_input is None and fasta_dir is None:
            typer.echo(
                "Error: Must provide either --fasta-input or --fasta-dir", err=True
            )
            raise typer.Exit(1)

        # Use fasta_dir if fasta_input not provided (backward compatibility)
        fasta_path = fasta_input if fasta_input else fasta_dir
        # Emit a one-time deprecation notice if using --fasta-dir
        if fasta_input is None and fasta_dir is not None:
            typer.echo(
                "[DEPRECATION] --fasta-dir is deprecated; use --fasta-input for files or directories.",
                err=True,
            )

        # Initialize pipeline
        pipeline = GuideRNAPipeline(output_dir)

        # Load sequences
        typer.echo(f"Loading sequences from {fasta_path}...")
        sequences = pipeline.load_sequences(fasta_path)

        if len(sequences) < 2:
            typer.echo(
                "Error: Need at least 2 sequences for conservation analysis", err=True
            )
            raise typer.Exit(1)

        typer.echo(f"Loaded {len(sequences)} sequences")

        # Collect required targets
        required_target_set = set()
        if required_targets:
            required_target_set.update(required_targets)
            typer.echo(f"Required targets from command line: {required_targets}")

        if required_targets_file:
            req_file = Path(required_targets_file)
            if not req_file.exists():
                typer.echo(
                    f"Error: Required targets file not found: {required_targets_file}",
                    err=True,
                )
                raise typer.Exit(1)
            with open(req_file) as f:
                file_targets = [line.strip() for line in f if line.strip()]
                required_target_set.update(file_targets)
                typer.echo(f"Loaded {len(file_targets)} required targets from file")

        # Validate required targets exist in sequences
        if required_target_set:
            missing_targets = required_target_set - set(sequences.keys())
            if missing_targets:
                typer.echo(
                    f"Warning: Required targets not found in input sequences: {missing_targets}",
                    err=True,
                )
                # Remove missing targets from required set
                required_target_set = required_target_set & set(sequences.keys())

            if required_target_set:
                typer.echo(
                    f"✓ Will prioritize {len(required_target_set)} required target sequences"
                )

        # Parse PAM patterns
        if pam_type.lower() == "spcas9":
            pam_patterns = ["NGG", "NAG"]
        else:
            # Custom PAM patterns
            pam_patterns = [p.strip().upper() for p in pam_type.split(",")]

        typer.echo(f"Using PAM patterns: {', '.join(pam_patterns)}")

        # Handle no_degenerate mode
        if no_degenerate:
            max_degenerate = 0
            typer.echo("✓ Non-degenerate mode enabled - exact sequences only")
            if required_target_set:
                typer.echo(
                    f"  Prioritizing 100% coverage of {len(required_target_set)} required targets"
                )

        # Handle perfect coverage mode
        if perfect_coverage:
            typer.echo("✓ Perfect coverage mode enabled")
            typer.echo(f"  Target: {min_coverage_per_target} guides per sequence")
            typer.echo(f"  Maximum total guides: {max_total_grnas}")
            typer.echo(
                f"  Conservation filter: {conservation:.1%} (guides below this won't be selected)"
            )
            if no_degenerate:
                typer.echo(
                    "  Note: Perfect coverage mode automatically uses exact sequences"
                )

        # Initialize guide RNA finder with all parameters
        finder = GuideRNAFinder(
            min_conservation=conservation,
            max_degenerate=max_degenerate,
            pam_patterns=pam_patterns,
            required_targets=list(required_target_set) if required_target_set else None,
            no_degenerate=no_degenerate,
            min_guides=min_guides,
            max_guides=max_guides,
            min_additional_coverage=min_additional_coverage,
            perfect_coverage=perfect_coverage,
            min_coverage_per_target=min_coverage_per_target,
            max_total_grnas=max_total_grnas,
            min_gc=min_gc,
            max_gc=max_gc,
            cluster_threshold=cluster_threshold,
        )

        # Find guide RNAs
        typer.echo("\nSearching for conserved guide RNAs...")
        guides = finder.find_guide_rnas(sequences, allow_mismatches=allow_mismatches)

        if not guides:
            typer.echo(
                "No conserved guide RNAs found with specified parameters", err=True
            )
            raise typer.Exit(1)

        typer.echo(f"Found {len(guides)} guide RNA candidates")

        # Provide feedback if fewer guides than requested
        if no_degenerate and len(guides) < min_guides:
            typer.echo(
                f"\n⚠️  Warning: Only {len(guides)} guides available (requested {min_guides})",
                err=True,
            )
            typer.echo("💡 Tip: Try one of these options to get more guides:", err=True)
            typer.echo(
                "   1. Remove --no-degenerate flag to allow degenerate bases", err=True
            )
            typer.echo("   2. Lower the --conservation threshold", err=True)
            typer.echo("   3. Increase --mismatches for more clustering", err=True)
            typer.echo("")  # Blank line for readability

        # Filter by minimum targets if specified
        if min_targets:
            guides = [g for g in guides if len(g.positions) >= min_targets]
            typer.echo(f"Filtered to {len(guides)} guides with ≥{min_targets} targets")

        # Limit to top N if specified
        if top_n and len(guides) > top_n:
            guides = guides[:top_n]
            typer.echo(f"Limiting output to top {top_n} guides")

        # Get coverage statistics for perfect coverage mode
        coverage_stats = None
        if perfect_coverage and guides:
            coverage_stats = finder._calculate_perfect_coverage_stats(guides, sequences)
            # Add configuration parameters to stats for output
            coverage_stats["min_coverage_per_target"] = min_coverage_per_target
            coverage_stats["max_total_grnas"] = max_total_grnas

        # Save results
        typer.echo("\nSaving results...")
        pipeline.save_results(
            guides=guides,
            sequences=sequences,
            pam_type=pam_type,
            conservation_threshold=conservation,
            max_degenerate=max_degenerate,
            perfect_coverage=perfect_coverage,
            coverage_stats=coverage_stats,
        )

        # Display summary
        typer.echo("\n" + "=" * 60)
        typer.echo("SUMMARY")
        typer.echo("=" * 60)
        typer.echo(f"Total guide RNAs found: {len(guides)}")

        if guides:
            avg_cons = sum(g.conservation for g in guides) / len(guides)
            typer.echo(f"Average conservation: {avg_cons:.1%}")

            no_degen = sum(1 for g in guides if g.degenerate_count == 0)
            typer.echo(f"Guides with no degeneracy: {no_degen}")

            # Show perfect coverage statistics
            if perfect_coverage and coverage_stats:
                typer.echo("\nPerfect Coverage Results:")
                typer.echo(
                    f"  Targets with sufficient coverage: {coverage_stats.get('covered_targets', 0)}/{coverage_stats.get('total_targets', 0)} ({coverage_stats.get('coverage_percentage', 0):.1%})"
                )
                typer.echo(
                    f"  Average guides per target: {coverage_stats.get('average_guides_per_target', 0):.1f}"
                )

            # Show top guide
            top_guide = guides[0]
            typer.echo("\nTop guide RNA:")
            typer.echo(f"  Sequence: {top_guide.consensus or top_guide.sequence}")
            typer.echo(f"  Conservation: {top_guide.conservation:.1%}")
            typer.echo(f"  Quality score: {top_guide.quality_score:.3f}")

        typer.echo(f"\nResults saved to: {pipeline.output_dir}")
        typer.echo("  - guide_rnas.csv: Detailed guide information")
        typer.echo("  - guide_rnas.fasta: Guide sequences in FASTA format")
        typer.echo("  - summary.txt: Analysis summary and statistics")
        typer.echo("  - target_mapping.tsv: Detailed target positions")

        if perfect_coverage and coverage_stats:
            typer.echo(
                "  - coverage_matrix.tsv: Binary coverage matrix (perfect coverage mode)"
            )
            typer.echo(
                "  - target_coverage_details.tsv: Per-target coverage statistics"
            )

    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback

            traceback.print_exc()
        raise typer.Exit(1)


if __name__ == "__main__":
    app()


# Allow calling `ppdesign grna` without explicitly specifying `main`
@app.callback(invoke_without_command=True)
def _default(
    ctx: typer.Context,
    fasta_input: str = typer.Option(
        None,
        "--fasta-input",
        "-f",
        help="Path to FASTA file or directory containing target sequences",
    ),
    fasta_dir: str = typer.Option(
        None,
        "--fasta-dir",
        help="[Deprecated: use --fasta-input] Directory containing FASTA files",
    ),
    output_dir: str = typer.Option(
        None, "--output-dir", "-o", help="Output directory for results"
    ),
    pam_type: str = typer.Option(
        "spCas9",
        "--pam-type",
        help="PAM type: spCas9 (NGG/NAG), custom patterns comma-separated",
    ),
    conservation: float = typer.Option(
        0.8,
        "--conservation",
        "-c",
        min=0.0,
        max=1.0,
        help="Minimum conservation threshold (0.0-1.0). In perfect coverage mode, guides must still meet this threshold to be considered.",
    ),
    max_degenerate: int = typer.Option(
        3,
        "--max-degenerate",
        "-d",
        min=0,
        help="Maximum number of degenerate bases allowed",
    ),
    no_degenerate: bool = typer.Option(
        False,
        "--no-degenerate",
        help="Disable degenerate bases - use exact sequences only (prioritizes required targets)",
    ),
    allow_mismatches: int = typer.Option(
        2,
        "--mismatches",
        "-m",
        min=0,
        help="Number of mismatches allowed when clustering guides",
    ),
    threads: int = typer.Option(
        1, "--threads", "-t", help="Number of threads (currently single-threaded)"
    ),
    top_n: Optional[int] = typer.Option(
        None, "--top-n", help="Output only top N guides by quality score"
    ),
    min_targets: Optional[int] = typer.Option(
        None, "--min-targets", help="Minimum number of target sequences required"
    ),
    required_targets: Optional[List[str]] = typer.Option(
        None,
        "--required-targets",
        help="Sequence IDs that must be targeted (space-separated)",
    ),
    required_targets_file: Optional[str] = typer.Option(
        None,
        "--required-targets-file",
        help="File with sequence IDs (one per line) that must be targeted",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose output"
    ),
    min_guides: int = typer.Option(
        50,
        "--min-guides",
        help="Target number of guides to return (returns all available if fewer exist)",
    ),
    max_guides: Optional[int] = typer.Option(
        None,
        "--max-guides",
        help="Maximum limit on guide count (default: no limit)",
    ),
    min_additional_coverage: float = typer.Option(
        0.05,
        "--min-additional-coverage",
        help="Coverage threshold for adding guides (0.0-1.0)",
    ),
    perfect_coverage: bool = typer.Option(
        False,
        "--perfect-coverage",
        help="Ensure each target has at least N perfect gRNA matches",
    ),
    min_coverage_per_target: int = typer.Option(
        3,
        "--min-coverage-per-target",
        help="Minimum perfect gRNAs per target (perfect coverage mode)",
    ),
    max_total_grnas: int = typer.Option(
        20, "--max-total-grnas", help="Max total gRNAs (perfect coverage mode)"
    ),
    min_gc: Optional[float] = typer.Option(
        None, "--min-gc", min=0.0, max=100.0, help="Minimum GC% for guides"
    ),
    max_gc: Optional[float] = typer.Option(
        None, "--max-gc", min=0.0, max=100.0, help="Maximum GC% for guides"
    ),
    cluster_threshold: Optional[float] = typer.Option(
        None,
        "--cluster-threshold",
        min=0.0,
        max=1.0,
        help="Similarity threshold for clustering sequences (0.0-1.0)",
    ),
):
    """Run gRNA design without specifying the 'main' subcommand."""
    if ctx.invoked_subcommand is None:
        if output_dir is None:
            raise typer.BadParameter("--output-dir is required")
        return main(
            fasta_input=fasta_input,
            fasta_dir=fasta_dir,
            output_dir=output_dir,
            pam_type=pam_type,
            conservation=conservation,
            max_degenerate=max_degenerate,
            no_degenerate=no_degenerate,
            allow_mismatches=allow_mismatches,
            threads=threads,
            top_n=top_n,
            min_targets=min_targets,
            required_targets=required_targets,
            required_targets_file=required_targets_file,
            verbose=verbose,
            min_guides=min_guides,
            max_guides=max_guides,
            min_additional_coverage=min_additional_coverage,
            perfect_coverage=perfect_coverage,
            min_coverage_per_target=min_coverage_per_target,
            max_total_grnas=max_total_grnas,
            min_gc=min_gc,
            max_gc=max_gc,
            cluster_threshold=cluster_threshold,
        )
    return None
