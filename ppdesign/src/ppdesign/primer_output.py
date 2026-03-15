"""Output formatter for primer pair results."""

import csv
from pathlib import Path
from typing import List
from .primer_types import PrimerPair


def write_primer_pairs_csv(pairs: List[PrimerPair], output_path: Path):
    """
    Write primer pairs to CSV format.

    Output columns:
    - pair_id: Unique identifier
    - forward_seq, forward_tm, forward_gc, forward_pos: Forward primer details
    - reverse_seq, reverse_tm, reverse_gc, reverse_pos: Reverse primer details
    - amplicon_size: Predicted amplicon size
    - tm_difference: Tm difference between primers
    - cross_dimer_score: Cross-dimer score
    - quality_score: Composite quality score
    - conservation: Average conservation

    Args:
        pairs: List of primer pairs to write
        output_path: Output CSV file path
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)

        # Write header
        writer.writerow(
            [
                "pair_id",
                "forward_seq",
                "forward_length",
                "forward_tm",
                "forward_gc",
                "forward_hairpin_dg",
                "forward_pos",
                "forward_strand",
                "reverse_seq",
                "reverse_length",
                "reverse_tm",
                "reverse_gc",
                "reverse_hairpin_dg",
                "reverse_pos",
                "reverse_strand",
                "amplicon_size",
                "tm_difference",
                "cross_dimer_score",
                "quality_score",
                "conservation",
            ]
        )

        def _fmt_dg(dg):
            return f"{dg:.2f}" if dg is not None else ""

        # Write data
        for i, pair in enumerate(pairs, 1):
            writer.writerow(
                [
                    f"pair_{i:03d}",
                    pair.forward.sequence,
                    len(pair.forward.sequence),
                    f"{pair.forward.tm:.2f}",
                    f"{pair.forward.gc_content:.2f}",
                    _fmt_dg(pair.forward.hairpin_dg),
                    pair.forward.position,
                    pair.forward.strand,
                    pair.reverse.sequence,
                    len(pair.reverse.sequence),
                    f"{pair.reverse.tm:.2f}",
                    f"{pair.reverse.gc_content:.2f}",
                    _fmt_dg(pair.reverse.hairpin_dg),
                    pair.reverse.position,
                    pair.reverse.strand,
                    pair.amplicon_size,
                    f"{pair.tm_difference:.2f}",
                    pair.cross_dimer_score,
                    f"{pair.quality_score:.2f}",
                    f"{pair.avg_conservation:.4f}",
                ]
            )


def write_primer_fasta(pairs: List[PrimerPair], output_path: Path):
    """
    Write primers to FASTA format.

    Each pair produces two entries:
    - >pair_XXX_forward
    - >pair_XXX_reverse

    Args:
        pairs: List of primer pairs to write
        output_path: Output FASTA file path
    """
    with open(output_path, "w") as f:
        for i, pair in enumerate(pairs, 1):
            pair_id = f"pair_{i:03d}"

            # Forward primer
            f.write(f">{pair_id}_forward\n")
            f.write(f"{pair.forward.sequence}\n")

            # Reverse primer
            f.write(f">{pair_id}_reverse\n")
            f.write(f"{pair.reverse.sequence}\n")


def write_amplicon_predictions(pairs: List[PrimerPair], output_path: Path):
    """
    Write amplicon predictions to TSV format.

    Output columns:
    - pair_id: Unique identifier
    - amplicon_start: Start position on forward strand
    - amplicon_end: End position on forward strand
    - amplicon_size: Size in bp
    - target_sequences: Comma-separated list of target sequence IDs

    Args:
        pairs: List of primer pairs to write
        output_path: Output TSV file path
    """
    with open(output_path, "w") as f:
        # Write header
        f.write("pair_id\tamplicon_start\tamplicon_end\tamplicon_size\ttarget_sequences\n")

        # Write data
        for i, pair in enumerate(pairs, 1):
            pair_id = f"pair_{i:03d}"
            start, end = pair.amplicon_range

            # Get common targets (sequences where both primers are found)
            fwd_targets = set(pair.forward.target_sequences)
            rev_targets = set(pair.reverse.target_sequences)
            common_targets = fwd_targets & rev_targets
            target_str = ",".join(sorted(common_targets))

            f.write(
                f"{pair_id}\t{start}\t{end}\t{pair.amplicon_size}\t{target_str}\n"
            )


def write_summary(pairs: List[PrimerPair], output_path: Path):
    """
    Write summary statistics to text file.

    Includes:
    - Total primer pairs
    - Amplicon size range
    - Average Tm difference
    - Top 10 primer pairs

    Args:
        pairs: List of primer pairs
        output_path: Output text file path
    """
    with open(output_path, "w") as f:
        f.write("Primer Pair Design Summary\n")
        f.write("=" * 50 + "\n\n")

        # Basic statistics
        f.write(f"Total primer pairs: {len(pairs)}\n\n")

        if pairs:
            # Amplicon size statistics
            amplicon_sizes = [p.amplicon_size for p in pairs]
            mean_amplicon = sum(amplicon_sizes) / len(amplicon_sizes)
            f.write("Amplicon Size Statistics:\n")
            f.write(f"  Range: {min(amplicon_sizes)}-{max(amplicon_sizes)} bp\n")
            f.write(f"  Mean: {mean_amplicon:.1f} bp\n\n")

            # Tm statistics
            tm_diffs = [p.tm_difference for p in pairs]
            mean_tm = sum(tm_diffs) / len(tm_diffs)
            f.write("Tm Difference Statistics:\n")
            f.write(f"  Range: {min(tm_diffs):.2f}-{max(tm_diffs):.2f}°C\n")
            f.write(f"  Mean: {mean_tm:.2f}°C\n\n")

            # Conservation statistics
            conservations = [p.avg_conservation for p in pairs]
            mean_conservation = sum(conservations) / len(conservations)
            f.write("Conservation Statistics:\n")
            f.write(
                f"  Range: {min(conservations):.2%}-{max(conservations):.2%}\n"
            )
            f.write(f"  Mean: {mean_conservation:.2%}\n\n")

            # Top 10 primer pairs
            f.write("Top 10 Primer Pairs (by quality score):\n")
            f.write("-" * 70 + "\n")
            for i, pair in enumerate(pairs[:10], 1):
                fwd = pair.forward
                rev = pair.reverse
                f.write(f"\n{i}. Score: {pair.quality_score:.1f}/100  "
                        f"Amplicon: {pair.amplicon_size}bp  "
                        f"Conservation: {pair.avg_conservation:.0%}\n")
                f.write(f"   Fwd: {fwd.sequence}  "
                        f"({len(fwd.sequence)}bp, Tm={fwd.tm:.1f}°C, "
                        f"GC={fwd.gc_content:.0f}%")
                if fwd.hairpin_dg is not None:
                    f.write(f", hairpin={fwd.hairpin_dg:.1f} kcal/mol")
                f.write(f", pos={fwd.position})\n")
                f.write(f"   Rev: {rev.sequence}  "
                        f"({len(rev.sequence)}bp, Tm={rev.tm:.1f}°C, "
                        f"GC={rev.gc_content:.0f}%")
                if rev.hairpin_dg is not None:
                    f.write(f", hairpin={rev.hairpin_dg:.1f} kcal/mol")
                f.write(f", pos={rev.position})\n")
                f.write(f"   ΔTm={pair.tm_difference:.1f}°C  "
                        f"cross-dimer={pair.cross_dimer_score}\n")
        else:
            f.write("No valid primer pairs found.\n")


def write_all_outputs(pairs: List[PrimerPair], output_dir: Path):
    """
    Write all output files to the output directory.

    Creates:
    - primer_pairs.csv: Detailed primer pair information
    - primer_pairs.fasta: FASTA format for synthesis
    - amplicon_predictions.tsv: Amplicon coordinates
    - summary.txt: Summary statistics

    Args:
        pairs: List of primer pairs
        output_dir: Output directory path
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    write_primer_pairs_csv(pairs, output_dir / "primer_pairs.csv")
    write_primer_fasta(pairs, output_dir / "primer_pairs.fasta")
    write_amplicon_predictions(pairs, output_dir / "amplicon_predictions.tsv")
    write_summary(pairs, output_dir / "summary.txt")
