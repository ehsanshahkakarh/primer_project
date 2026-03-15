#!/usr/bin/env python3
"""
Validate guide RNA alignments against target sequences.

This script performs local alignments of guide RNAs against target sequences
to validate conservation calculations and analyze hit patterns.
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
from dataclasses import dataclass, field


@dataclass
class AlignmentResult:
    """Store alignment results for a guide RNA against a sequence."""

    seq_id: str
    score: int
    identity: float
    coverage: float
    start_pos: int
    end_pos: int
    mismatches: List[int] = field(default_factory=list)
    cigar: str = ""
    ref_aligned: str = ""
    query_aligned: str = ""

    @property
    def has_central_mismatches(self) -> bool:
        """Check if mismatches are in central region (positions 6-14 of 20bp guide)."""
        if not self.mismatches:
            return False
        central_positions = set(range(6, 15))  # 0-indexed positions 6-14
        return bool(set(self.mismatches) & central_positions)

    @property
    def has_peripheral_mismatches(self) -> bool:
        """Check if mismatches are in peripheral regions (positions 0-5 or 15-19)."""
        if not self.mismatches:
            return False
        peripheral_positions = set(range(0, 6)) | set(range(15, 20))
        return bool(set(self.mismatches) & peripheral_positions)


class GuideRNAValidator:
    """Validate guide RNA alignments and conservation."""

    def __init__(self, guide_csv: Path, guide_fasta: Path, target_fasta: Path):
        """Initialize validator with input files."""
        self.console = Console()
        self.guide_csv = guide_csv
        self.guide_fasta = guide_fasta
        self.target_fasta = target_fasta

        # Load data
        self.guides_df = pd.read_csv(guide_csv)
        self.guides = self._load_guides()
        self.targets = self._load_targets()

        # Alignment parameters
        self.matrix = parasail.matrix_create("ACGT", 2, -1)  # Match: 2, Mismatch: -1
        self.gap_open = 3
        self.gap_extend = 1

    def _load_guides(self) -> Dict[str, str]:
        """Load guide RNA sequences from FASTA."""
        guides = {}
        for record in SeqIO.parse(self.guide_fasta, "fasta"):
            guides[record.id] = str(record.seq).upper()
        return guides

    def _load_targets(self) -> Dict[str, str]:
        """Load target sequences from FASTA."""
        targets = {}
        for record in SeqIO.parse(self.target_fasta, "fasta"):
            targets[record.id] = str(record.seq).upper()
        return targets

    def align_guide_to_target(self, guide_seq: str, target_seq: str) -> AlignmentResult:
        """Perform semi-global alignment of guide RNA to target sequence."""
        # Use semi-global alignment (guide must align completely, target can have overhangs)
        result = parasail.sg_trace_scan_sat(
            guide_seq, target_seq, self.gap_open, self.gap_extend, self.matrix
        )

        # Extract alignment details
        cigar = result.cigar
        trace = result.traceback

        # Parse alignment to find mismatches
        ref_aligned = trace.ref
        query_aligned = trace.query

        mismatches = []
        query_pos = 0
        for i, (r, q) in enumerate(zip(ref_aligned, query_aligned)):
            if q != "-":  # Not a gap in query
                if r != q and r != "-":  # Mismatch (not gap)
                    mismatches.append(query_pos)
                query_pos += 1

        # Calculate identity and coverage
        matches = sum(
            1 for r, q in zip(ref_aligned, query_aligned) if r == q and r != "-"
        )
        identity = (matches / len(guide_seq)) * 100 if len(guide_seq) > 0 else 0
        coverage = (len([c for c in query_aligned if c != "-"]) / len(guide_seq)) * 100

        # Get alignment positions - parasail returns end positions
        # For start position, we need to calculate from alignment length
        end_pos = result.end_ref
        start_pos = end_pos - len([c for c in ref_aligned if c != "-"]) + 1

        return AlignmentResult(
            seq_id="",  # Will be set by caller
            score=result.score,
            identity=identity,
            coverage=coverage,
            start_pos=start_pos,
            end_pos=end_pos,
            mismatches=mismatches,
            cigar=cigar.decode("utf-8") if isinstance(cigar, bytes) else str(cigar),
            ref_aligned=ref_aligned,
            query_aligned=query_aligned,
        )

    def validate_guide(self, guide_id: str, guide_seq: str) -> Dict:
        """Validate a single guide RNA against all target sequences."""
        results = {
            "guide_id": guide_id,
            "sequence": guide_seq,
            "length": len(guide_seq),
            "alignments": [],
            "hits": 0,
            "perfect_matches": 0,
            "central_mismatches": 0,
            "peripheral_mismatches": 0,
            "avg_identity": 0,
            "conservation": 0,
            "conservation_by_hits": 0,  # New: total hits method
            "sequences_with_hits": set(),  # New: track unique sequences
        }

        # Align against all targets
        for target_id, target_seq in self.targets.items():
            found_in_sequence = False
            # Check both strands
            for strand, seq in [
                ("+", target_seq),
                ("-", str(Seq(target_seq).reverse_complement())),
            ]:
                alignment = self.align_guide_to_target(guide_seq, seq)

                # Only consider significant alignments (>80% identity)
                if alignment.identity >= 80:
                    alignment.seq_id = f"{target_id}_{strand}"
                    results["alignments"].append(alignment)
                    results["hits"] += 1
                    found_in_sequence = True

                    if alignment.identity == 100:
                        results["perfect_matches"] += 1

                    if alignment.has_central_mismatches:
                        results["central_mismatches"] += 1
                    elif alignment.has_peripheral_mismatches:
                        results["peripheral_mismatches"] += 1

            if found_in_sequence:
                results["sequences_with_hits"].add(target_id)

        # Calculate statistics
        if results["alignments"]:
            results["avg_identity"] = np.mean(
                [a.identity for a in results["alignments"]]
            )
            # Original method: sequences with at least one hit / total sequences
            results["conservation"] = (
                len(results["sequences_with_hits"]) / len(self.targets)
            ) * 100
            # Alternative method: total hits / (sequences * 2 strands)
            results["conservation_by_hits"] = (
                results["hits"] / (len(self.targets) * 2)
            ) * 100

        return results

    def validate_all(self) -> List[Dict]:
        """Validate all guide RNAs."""
        self.console.print("\n[bold cyan]Validating Guide RNAs[/bold cyan]")

        all_results = []
        for guide_id, guide_seq in track(
            self.guides.items(), description="Processing guides..."
        ):
            result = self.validate_guide(guide_id, guide_seq)
            all_results.append(result)

        return all_results

    def generate_report(self, results: List[Dict], output_file: Path):
        """Generate comprehensive validation report."""
        self.console.print("\n[bold green]Validation Report[/bold green]\n")

        # Create summary table
        table = Table(title="Guide RNA Validation Summary")
        table.add_column("Guide ID", style="cyan")
        table.add_column("Sequence", style="white")
        table.add_column("Reported\nConservation", justify="right", style="yellow")
        table.add_column("Validated\nConservation", justify="right", style="green")
        table.add_column("Sequences\nwith Hits", justify="right")
        table.add_column("Total\nHits", justify="right")
        table.add_column("Perfect\nMatches", justify="right")
        table.add_column("Avg\nIdentity", justify="right")
        table.add_column("Central\nMismatches", justify="right", style="red")

        # Add rows
        for i, result in enumerate(results):
            guide_id = result["guide_id"]

            # Get reported conservation from CSV
            reported_cons = self.guides_df[self.guides_df["Guide_ID"] == guide_id][
                "Conservation"
            ].iloc[0]
            if isinstance(reported_cons, str) and "%" in reported_cons:
                reported_cons = float(reported_cons.strip("%"))

            table.add_row(
                guide_id,
                (
                    result["sequence"][:10] + "..."
                    if len(result["sequence"]) > 10
                    else result["sequence"]
                ),
                f"{reported_cons:.1f}%",
                f"{result['conservation']:.1f}%",
                str(len(result.get("sequences_with_hits", []))),
                str(result["hits"]),
                str(result["perfect_matches"]),
                (
                    f"{result['avg_identity']:.1f}%"
                    if result["avg_identity"] > 0
                    else "N/A"
                ),
                str(result["central_mismatches"]),
            )

        self.console.print(table)

        # Generate detailed JSON report
        detailed_report = {
            "summary": {
                "total_guides": len(results),
                "total_targets": len(self.targets),
                "avg_hits_per_guide": np.mean([r["hits"] for r in results]),
                "guides_with_perfect_matches": sum(
                    1 for r in results if r["perfect_matches"] > 0
                ),
                "guides_with_central_mismatches": sum(
                    1 for r in results if r["central_mismatches"] > 0
                ),
                "guides_with_peripheral_mismatches": sum(
                    1 for r in results if r["peripheral_mismatches"] > 0
                ),
            },
            "guides": [],
        }

        for result in results:
            guide_data = {
                "guide_id": result["guide_id"],
                "sequence": result["sequence"],
                "validation": {
                    "hits": result["hits"],
                    "conservation": result["conservation"],
                    "avg_identity": result["avg_identity"],
                    "perfect_matches": result["perfect_matches"],
                    "central_mismatches": result["central_mismatches"],
                    "peripheral_mismatches": result["peripheral_mismatches"],
                },
                "alignments": [],
            }

            # Add top alignments
            for alignment in sorted(
                result["alignments"], key=lambda x: x.score, reverse=True
            )[:5]:
                guide_data["alignments"].append(
                    {
                        "target": alignment.seq_id,
                        "score": alignment.score,
                        "identity": alignment.identity,
                        "position": f"{alignment.start_pos}-{alignment.end_pos}",
                        "mismatches": alignment.mismatches,
                        "mismatch_type": (
                            "central"
                            if alignment.has_central_mismatches
                            else (
                                "peripheral"
                                if alignment.has_peripheral_mismatches
                                else "none"
                            )
                        ),
                    }
                )

            detailed_report["guides"].append(guide_data)

        # Save JSON report
        with open(output_file, "w") as f:
            json.dump(detailed_report, f, indent=2)

        self.console.print(f"\n[green]Detailed report saved to:[/green] {output_file}")

        # Print key findings
        self.console.print("\n[bold]Key Findings:[/bold]")

        # Compare reported vs validated conservation
        conservation_diffs = []
        for i, result in enumerate(results):
            guide_id = result["guide_id"]
            reported = self.guides_df[self.guides_df["Guide_ID"] == guide_id][
                "Conservation"
            ].iloc[0]
            if isinstance(reported, str) and "%" in reported:
                reported = float(reported.strip("%"))
            validated = result["conservation"]
            diff = abs(reported - validated)
            conservation_diffs.append(diff)

        avg_diff = np.mean(conservation_diffs)
        max_diff = np.max(conservation_diffs)

        self.console.print(
            f"• Average conservation difference: [yellow]{avg_diff:.2f}%[/yellow]"
        )
        self.console.print(
            f"• Maximum conservation difference: [yellow]{max_diff:.2f}%[/yellow]"
        )

        # Mismatch analysis
        total_with_mismatches = sum(
            1
            for r in results
            if r["central_mismatches"] > 0 or r["peripheral_mismatches"] > 0
        )
        if total_with_mismatches > 0:
            central_ratio = (
                sum(r["central_mismatches"] for r in results) / total_with_mismatches
            )
            peripheral_ratio = (
                sum(r["peripheral_mismatches"] for r in results) / total_with_mismatches
            )

            self.console.print(
                f"• Guides with mismatches: [yellow]{total_with_mismatches}/{len(results)}[/yellow]"
            )
            self.console.print(
                f"• Central mismatches (critical): [red]{central_ratio:.1f} per guide[/red]"
            )
            self.console.print(
                f"• Peripheral mismatches: [yellow]{peripheral_ratio:.1f} per guide[/yellow]"
            )

        return detailed_report


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Validate guide RNA alignments")
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
        default=Path("grna_validation_report.json"),
        help="Output file for validation report (default: grna_validation_report.json)",
    )

    args = parser.parse_args()

    # Validate input files
    for file_path in [args.guide_csv, args.guide_fasta, args.target_fasta]:
        if not file_path.exists():
            print(f"Error: File not found: {file_path}")
            return 1

    # Run validation
    validator = GuideRNAValidator(args.guide_csv, args.guide_fasta, args.target_fasta)
    results = validator.validate_all()
    validator.generate_report(results, args.output)

    return 0


if __name__ == "__main__":
    exit(main())
