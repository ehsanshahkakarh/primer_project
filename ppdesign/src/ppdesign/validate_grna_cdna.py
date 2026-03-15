#!/usr/bin/env python3
"""
cDNA-aware validation for CRISPR guide RNA alignments.

This module provides specialized validation for guide RNAs targeting single-stranded
cDNA sequences, with strand bias analysis critical for cDNA applications.

Key Features:
    - Separate tracking of forward (mRNA-sense) and reverse (antisense) strand hits
    - Strand bias detection to identify preferential targeting patterns
    - Conservation calculation across target sequences
    - Parasail-based semi-global alignment for accurate guide matching

Usage:
    pixi run validate-grna-cdna \\
        --guide-csv results/guide_rnas.csv \\
        --guide-fasta results/guide_rnas.fasta \\
        --target-fasta targets.fna \\
        --output cdna_validation.json

Output:
    JSON report containing:
    - Per-guide strand bias metrics
    - Conservation percentages
    - Forward vs reverse hit counts
    - Summary statistics for cDNA targeting

Author: PPDesign Development Team
License: MIT
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import parasail
from rich.console import Console
from rich.table import Table
from rich.progress import track


class CDNAGuideValidator:
    """Validate guide RNAs for cDNA applications."""

    def __init__(
        self,
        guide_csv: Path,
        guide_fasta: Path,
        target_fasta: Path,
        required_targets: List[str] = None,
    ):
        """Initialize validator.

        Args:
            guide_csv: Path to guide RNA CSV file
            guide_fasta: Path to guide RNA FASTA file
            target_fasta: Path to target sequences FASTA file
            required_targets: List of sequence IDs that must be perfectly matched
        """
        self.console = Console()
        self.guides_df = pd.read_csv(guide_csv)
        self.guides = {r.id: str(r.seq) for r in SeqIO.parse(guide_fasta, "fasta")}
        self.targets = {r.id: str(r.seq) for r in SeqIO.parse(target_fasta, "fasta")}
        self.required_targets = set(required_targets) if required_targets else set()

        # Validate required targets exist
        if self.required_targets:
            missing = self.required_targets - set(self.targets.keys())
            if missing:
                self.console.print(
                    f"[red]Warning: Required targets not found: {missing}[/red]"
                )
                self.required_targets = self.required_targets & set(self.targets.keys())

        # Alignment parameters
        self.matrix = parasail.matrix_create("ACGT", 3, -2)
        self.gap_open = 5
        self.gap_extend = 2

    def validate_guide(self, guide_id: str, guide_seq: str) -> Dict:
        """Validate a guide RNA with strand awareness."""
        results = {
            "guide_id": guide_id,
            "sequence": guide_seq,
            "forward_hits": 0,
            "reverse_hits": 0,
            "total_hits": 0,
            "perfect_matches": 0,
            "conservation": 0,
            "strand_bias": None,
            "sequences_with_hits": set(),
            "perfect_match_targets": set(),
            "required_targets_hit": set(),
            "required_targets_perfect": set(),
        }

        # Check against all targets
        for target_id, target_seq in self.targets.items():
            # Forward strand (mRNA-sense)
            result_fwd = parasail.sg_trace_scan_sat(
                guide_seq, target_seq, self.gap_open, self.gap_extend, self.matrix
            )

            if result_fwd.score >= len(guide_seq) * 2:  # Good match
                results["forward_hits"] += 1
                results["sequences_with_hits"].add(target_id)

                # Check for perfect match (score = length * match_score)
                if (
                    result_fwd.score == len(guide_seq) * 3
                ):  # Perfect match (3 points per match)
                    results["perfect_matches"] += 1
                    results["perfect_match_targets"].add(target_id)
                    if target_id in self.required_targets:
                        results["required_targets_perfect"].add(target_id)

                # Track any hit for required targets
                if target_id in self.required_targets:
                    results["required_targets_hit"].add(target_id)

            # Reverse strand (antisense)
            rev_seq = str(Seq(target_seq).reverse_complement())
            result_rev = parasail.sg_trace_scan_sat(
                guide_seq, rev_seq, self.gap_open, self.gap_extend, self.matrix
            )

            if result_rev.score >= len(guide_seq) * 2:  # Good match
                results["reverse_hits"] += 1
                results["sequences_with_hits"].add(target_id)

                # Check for perfect match
                if result_rev.score == len(guide_seq) * 3:
                    if (
                        target_id not in results["perfect_match_targets"]
                    ):  # Avoid double counting
                        results["perfect_matches"] += 1
                        results["perfect_match_targets"].add(target_id)
                        if target_id in self.required_targets:
                            results["required_targets_perfect"].add(target_id)

                # Track any hit for required targets
                if target_id in self.required_targets:
                    results["required_targets_hit"].add(target_id)

        # Calculate metrics
        results["total_hits"] = results["forward_hits"] + results["reverse_hits"]
        results["conservation"] = (
            len(results["sequences_with_hits"]) / len(self.targets)
        ) * 100

        # Determine strand bias
        if results["forward_hits"] > results["reverse_hits"] * 1.5:
            results["strand_bias"] = "forward (mRNA-sense)"
        elif results["reverse_hits"] > results["forward_hits"] * 1.5:
            results["strand_bias"] = "reverse (antisense)"
        else:
            results["strand_bias"] = "balanced"

        # Convert sets to counts for JSON serialization
        results["sequences_with_hits"] = len(results["sequences_with_hits"])
        results["perfect_match_count"] = len(results["perfect_match_targets"])
        results["required_targets_hit_count"] = len(results["required_targets_hit"])
        results["required_targets_perfect_count"] = len(
            results["required_targets_perfect"]
        )

        # Store the actual IDs for reporting
        results["perfect_match_targets"] = list(results["perfect_match_targets"])
        results["required_targets_hit"] = list(results["required_targets_hit"])
        results["required_targets_perfect"] = list(results["required_targets_perfect"])

        return results

    def generate_report(self, results: List[Dict], output_file: Path):
        """Generate cDNA-specific report."""
        self.console.print("\n[bold green]cDNA Validation Report[/bold green]\n")

        # Summary table
        table = Table(title="Guide RNA Validation for cDNA Targets")
        table.add_column("Guide ID", style="cyan")
        table.add_column("Conservation", justify="right", style="yellow")
        table.add_column("Perfect\nMatches", justify="right", style="green")
        table.add_column("Forward\nHits", justify="right")
        table.add_column("Reverse\nHits", justify="right")
        table.add_column("Strand Bias", style="blue")

        if self.required_targets:
            table.add_column("Required\nTargets", justify="right", style="magenta")

        for result in results:
            row_data = [
                result["guide_id"],
                f"{result['conservation']:.1f}%",
                str(result.get("perfect_matches", 0)),
                str(result["forward_hits"]),
                str(result["reverse_hits"]),
                result["strand_bias"],
            ]

            if self.required_targets:
                req_info = f"{result.get('required_targets_perfect_count', 0)}/{len(self.required_targets)}"
                row_data.append(req_info)

            table.add_row(*row_data)

        self.console.print(table)

        # Check required targets coverage
        if self.required_targets:
            self.console.print("\n[bold]Required Targets Analysis:[/bold]")
            all_required_perfect = set()
            all_required_hit = set()

            for result in results:
                all_required_perfect.update(result.get("required_targets_perfect", []))
                all_required_hit.update(result.get("required_targets_hit", []))

            missing_perfect = self.required_targets - all_required_perfect
            missing_any = self.required_targets - all_required_hit

            if missing_perfect:
                self.console.print(
                    f"[yellow]⚠ Required targets without perfect match: {missing_perfect}[/yellow]"
                )
            else:
                self.console.print(
                    "[green]✓ All required targets have perfect matches[/green]"
                )

            if missing_any:
                self.console.print(
                    f"[red]✗ Required targets with no hits: {missing_any}[/red]"
                )

        # Save JSON report
        report = {
            "target_type": "cDNA",
            "summary": {
                "total_guides": len(results),
                "total_targets": len(self.targets),
                "required_targets": (
                    list(self.required_targets) if self.required_targets else []
                ),
                "guides_with_forward_bias": sum(
                    1 for r in results if "forward" in r["strand_bias"]
                ),
                "guides_with_reverse_bias": sum(
                    1 for r in results if "reverse" in r["strand_bias"]
                ),
                "avg_conservation": np.mean([r["conservation"] for r in results]),
                "total_perfect_matches": sum(
                    r.get("perfect_matches", 0) for r in results
                ),
            },
            "guides": results,
        }

        if self.required_targets:
            all_req_perfect = set()
            all_req_hit = set()
            for r in results:
                all_req_perfect.update(r.get("required_targets_perfect", []))
                all_req_hit.update(r.get("required_targets_hit", []))

            report["summary"]["required_targets_with_perfect_match"] = list(
                all_req_perfect
            )
            report["summary"]["required_targets_with_any_hit"] = list(all_req_hit)
            report["summary"]["required_targets_missing"] = list(
                self.required_targets - all_req_hit
            )

        with open(output_file, "w") as f:
            json.dump(report, f, indent=2)

        self.console.print(f"\n[green]Report saved to:[/green] {output_file}")

        # Print insights
        self.console.print("\n[bold]cDNA-Specific Insights:[/bold]")
        forward_biased = sum(1 for r in results if "forward" in r["strand_bias"])
        reverse_biased = sum(1 for r in results if "reverse" in r["strand_bias"])

        self.console.print(
            f"• Guides targeting mRNA-sense strand: [cyan]{forward_biased}[/cyan]"
        )
        self.console.print(
            f"• Guides targeting antisense strand: [cyan]{reverse_biased}[/cyan]"
        )

        if forward_biased > reverse_biased:
            self.console.print(
                "• [yellow]Most guides target the mRNA-sense strand (expected for cDNA)[/yellow]"
            )
        elif reverse_biased > forward_biased:
            self.console.print(
                "• [yellow]Most guides target the antisense strand[/yellow]"
            )


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description="cDNA-aware guide RNA validation with required target support",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic validation
  %(prog)s --guide-csv guides.csv --guide-fasta guides.fasta --target-fasta targets.fna
  
  # Require perfect matches for specific sequences
  %(prog)s --guide-csv guides.csv --guide-fasta guides.fasta \\
           --target-fasta targets.fna --required-targets seq1 seq2 seq3
  
  # Read required targets from file
  %(prog)s --guide-csv guides.csv --guide-fasta guides.fasta \\
           --target-fasta targets.fna --required-targets-file critical_seqs.txt
        """,
    )
    parser.add_argument(
        "--guide-csv", type=Path, required=True, help="Path to guide RNA CSV file"
    )
    parser.add_argument(
        "--guide-fasta", type=Path, required=True, help="Path to guide RNA FASTA file"
    )
    parser.add_argument(
        "--target-fasta",
        type=Path,
        required=True,
        help="Path to target sequences FASTA file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cdna_validation.json"),
        help="Output JSON report file (default: cdna_validation.json)",
    )
    parser.add_argument(
        "--required-targets",
        nargs="+",
        help="Sequence IDs that must be perfectly targeted",
    )
    parser.add_argument(
        "--required-targets-file",
        type=Path,
        help="File with sequence IDs (one per line) that must be perfectly targeted",
    )

    args = parser.parse_args()

    # Validate files
    for f in [args.guide_csv, args.guide_fasta, args.target_fasta]:
        if not f.exists():
            print(f"Error: File not found: {f}")
            return 1

    # Collect required targets
    required_targets = []
    if args.required_targets:
        required_targets.extend(args.required_targets)

    if args.required_targets_file:
        if not args.required_targets_file.exists():
            print(
                f"Error: Required targets file not found: {args.required_targets_file}"
            )
            return 1
        with open(args.required_targets_file) as f:
            required_targets.extend([line.strip() for line in f if line.strip()])

    # Remove duplicates
    required_targets = list(set(required_targets)) if required_targets else None

    # Run validation
    validator = CDNAGuideValidator(
        args.guide_csv, args.guide_fasta, args.target_fasta, required_targets
    )

    results = []
    for guide_id, guide_seq in track(
        validator.guides.items(), description="Validating..."
    ):
        results.append(validator.validate_guide(guide_id, guide_seq))

    validator.generate_report(results, args.output)
    return 0


if __name__ == "__main__":
    exit(main())
