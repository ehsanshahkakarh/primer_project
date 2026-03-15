#!/usr/bin/env python3
"""
Thermodynamic-enhanced validation for CRISPR guide RNAs.

This module extends validation with thermodynamic analysis to predict:
- Binding stability (ΔG, ΔH, ΔS)
- Melting temperature (Tm) 
- Off-target binding probabilities
- Specificity scores

Usage:
    pixi run validate-grna-thermo \
        --guide-csv results/guide_rnas.csv \
        --guide-fasta results/guide_rnas.fasta \
        --target-fasta targets.fna \
        --output thermo_validation.json
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import parasail
from rich.console import Console
from rich.table import Table
from rich.progress import track
import matplotlib.pyplot as plt
import seaborn as sns

from .thermodynamics import (
    ThermodynamicCalculator,
)


class ThermoGuideValidator:
    """Validate guide RNAs with thermodynamic analysis."""

    def __init__(
        self,
        guide_csv: Path,
        guide_fasta: Path,
        target_fasta: Path,
        required_targets: List[str] = None,
        na_conc: float = 0.150,  # 150 mM - CRISPR standard
        mg_conc: float = 0.010,  # 10 mM - CRISPR standard
        temperature: float = 37.0,
    ):
        """Initialize validator with thermodynamic parameters.

        Args:
            guide_csv: Path to guide RNA CSV file
            guide_fasta: Path to guide RNA FASTA file
            target_fasta: Path to target sequences FASTA file
            required_targets: List of sequence IDs that must be perfectly matched
            na_conc: Sodium concentration (M) - Default 150 mM for CRISPR
            mg_conc: Magnesium concentration (M) - Default 10 mM for CRISPR
            temperature: Temperature (°C) - Default 37°C

        Note: Default buffer conditions match standard CRISPR-Cas9 in vitro
        cleavage assays (150 mM NaCl, 10 mM MgCl2, 37°C, pH 7.5)
        """
        self.console = Console()
        self.guides_df = pd.read_csv(guide_csv)
        self.guides = {r.id: str(r.seq) for r in SeqIO.parse(guide_fasta, "fasta")}
        self.targets = {r.id: str(r.seq) for r in SeqIO.parse(target_fasta, "fasta")}
        self.required_targets = set(required_targets) if required_targets else set()

        # Initialize thermodynamic calculator
        self.thermo_calc = ThermodynamicCalculator(
            na_conc=na_conc, mg_conc=mg_conc, temperature=temperature
        )

        # Alignment parameters for initial screening
        self.matrix = parasail.matrix_create("ACGT", 3, -2)
        self.gap_open = 5
        self.gap_extend = 2

    def extract_target_region(
        self, target_seq: str, hit_pos: int, guide_len: int = 20
    ) -> str:
        """Extract the 20bp target region around a hit position."""
        # Ensure we get a 20bp region
        start = max(0, hit_pos)
        end = min(len(target_seq), start + guide_len)

        # Adjust if near the end
        if end - start < guide_len:
            start = max(0, end - guide_len)

        return target_seq[start:end]

    def find_guide_targets(
        self, guide_seq: str
    ) -> Dict[str, List[Tuple[str, int, str]]]:
        """Find all target sites for a guide RNA.

        Returns:
            Dict mapping target_id to list of (target_region, position, strand)
        """
        targets_found = {}

        for target_id, target_seq in self.targets.items():
            hits = []

            # Check forward strand
            result_fwd = parasail.sg_trace_scan_sat(
                guide_seq, target_seq, self.gap_open, self.gap_extend, self.matrix
            )

            if result_fwd.score >= len(guide_seq) * 2:  # Good match
                # Extract target region
                target_region = self.extract_target_region(
                    target_seq, result_fwd.end_ref - len(guide_seq) + 1
                )
                hits.append(
                    (target_region, result_fwd.end_ref - len(guide_seq) + 1, "+")
                )

            # Check reverse strand
            rev_seq = str(Seq(target_seq).reverse_complement())
            result_rev = parasail.sg_trace_scan_sat(
                guide_seq, rev_seq, self.gap_open, self.gap_extend, self.matrix
            )

            if result_rev.score >= len(guide_seq) * 2:  # Good match
                target_region = self.extract_target_region(
                    rev_seq, result_rev.end_ref - len(guide_seq) + 1
                )
                hits.append(
                    (target_region, result_rev.end_ref - len(guide_seq) + 1, "-")
                )

            if hits:
                targets_found[target_id] = hits

        return targets_found

    def validate_guide_thermo(self, guide_id: str, guide_seq: str) -> Dict:
        """Validate a guide RNA with thermodynamic analysis."""

        # Find all target sites
        target_sites = self.find_guide_targets(guide_seq)

        # Thermodynamic analysis for all targets
        thermo_profiles = []
        perfect_matches = []
        required_hits = []

        for target_id, sites in target_sites.items():
            for target_region, pos, strand in sites:
                # Calculate thermodynamic profile
                profile = self.thermo_calc.calculate_thermodynamics(
                    guide_seq, target_region
                )

                # Check if perfect match
                is_perfect = profile.mismatch_penalty == 0
                if is_perfect:
                    perfect_matches.append(target_id)
                    if target_id in self.required_targets:
                        required_hits.append(target_id)

                thermo_profiles.append(
                    {
                        "target_id": target_id,
                        "position": pos,
                        "strand": strand,
                        "delta_g": profile.delta_g,
                        "delta_h": profile.delta_h,
                        "delta_s": profile.delta_s,
                        "tm": profile.tm,
                        "binding_efficiency": profile.binding_efficiency,
                        "stability_score": profile.stability_score,
                        "is_perfect": is_perfect,
                    }
                )

        # Calculate summary statistics
        if thermo_profiles:
            avg_tm = np.mean([p["tm"] for p in thermo_profiles])
            avg_delta_g = np.mean([p["delta_g"] for p in thermo_profiles])
            avg_binding_eff = np.mean(
                [p["binding_efficiency"] for p in thermo_profiles]
            )
            avg_stability = np.mean([p["stability_score"] for p in thermo_profiles])

            # Find most stable target
            best_profile = min(thermo_profiles, key=lambda x: x["delta_g"])

            # Calculate specificity (ratio of best to average binding)
            if len(thermo_profiles) > 1:
                other_binding = [
                    p["binding_efficiency"]
                    for p in thermo_profiles
                    if p != best_profile
                ]
                specificity = best_profile["binding_efficiency"] / (
                    np.mean(other_binding) + best_profile["binding_efficiency"]
                )
            else:
                specificity = 1.0
        else:
            avg_tm = 0
            avg_delta_g = 0
            avg_binding_eff = 0
            avg_stability = 0
            best_profile = None
            specificity = 0

        result = {
            "guide_id": guide_id,
            "sequence": guide_seq,
            "total_targets": len(target_sites),
            "total_binding_sites": sum(len(sites) for sites in target_sites.values()),
            "perfect_matches": len(set(perfect_matches)),
            "conservation": (len(target_sites) / len(self.targets)) * 100,
            "thermodynamics": {
                "average_tm": avg_tm,
                "average_delta_g": avg_delta_g,
                "average_binding_efficiency": avg_binding_eff,
                "average_stability_score": avg_stability,
                "specificity_score": specificity,
                "best_target": best_profile if best_profile else None,
                "profiles": thermo_profiles[:10],  # Limit to top 10 for space
            },
        }

        # Add required targets info
        if self.required_targets:
            result["required_targets"] = {
                "total": len(self.required_targets),
                "hit": len(set(required_hits)),
                "perfect": len(set(required_hits)),
                "missing": list(self.required_targets - set(required_hits)),
            }

        return result

    def generate_thermo_plots(self, results: List[Dict], output_dir: Path):
        """Generate thermodynamic analysis plots."""
        output_dir.mkdir(exist_ok=True)

        # Extract data for plotting
        guide_ids = [r["guide_id"] for r in results]
        avg_tms = [r["thermodynamics"]["average_tm"] for r in results]
        avg_dgs = [r["thermodynamics"]["average_delta_g"] for r in results]
        specificities = [r["thermodynamics"]["specificity_score"] for r in results]
        stabilities = [r["thermodynamics"]["average_stability_score"] for r in results]

        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # 1. Melting Temperature Distribution
        axes[0, 0].hist(avg_tms, bins=20, edgecolor="black", alpha=0.7)
        axes[0, 0].set_xlabel("Average Tm (°C)")
        axes[0, 0].set_ylabel("Number of Guides")
        axes[0, 0].set_title("Melting Temperature Distribution")
        axes[0, 0].axvline(
            np.mean(avg_tms),
            color="red",
            linestyle="--",
            label=f"Mean: {np.mean(avg_tms):.1f}°C",
        )
        axes[0, 0].legend()

        # 2. Free Energy Distribution
        axes[0, 1].hist(avg_dgs, bins=20, edgecolor="black", alpha=0.7, color="green")
        axes[0, 1].set_xlabel("Average ΔG (kcal/mol)")
        axes[0, 1].set_ylabel("Number of Guides")
        axes[0, 1].set_title("Free Energy Distribution")
        axes[0, 1].axvline(
            np.mean(avg_dgs),
            color="red",
            linestyle="--",
            label=f"Mean: {np.mean(avg_dgs):.1f} kcal/mol",
        )
        axes[0, 1].legend()

        # 3. Specificity vs Stability
        scatter = axes[1, 0].scatter(
            stabilities, specificities, c=avg_tms, cmap="coolwarm", alpha=0.6
        )
        axes[1, 0].set_xlabel("Stability Score")
        axes[1, 0].set_ylabel("Specificity Score")
        axes[1, 0].set_title("Specificity vs Stability")
        plt.colorbar(scatter, ax=axes[1, 0], label="Tm (°C)")

        # 4. Top 10 Guides by Combined Score
        combined_scores = [
            r["thermodynamics"]["average_stability_score"]
            * r["thermodynamics"]["specificity_score"]
            for r in results
        ]
        top_indices = np.argsort(combined_scores)[-10:]
        top_guides = [guide_ids[i] for i in top_indices]
        top_scores = [combined_scores[i] for i in top_indices]

        axes[1, 1].barh(range(len(top_guides)), top_scores)
        axes[1, 1].set_yticks(range(len(top_guides)))
        axes[1, 1].set_yticklabels(top_guides)
        axes[1, 1].set_xlabel("Combined Score (Stability × Specificity)")
        axes[1, 1].set_title("Top 10 Guides by Thermodynamic Score")

        plt.tight_layout()
        plt.savefig(
            output_dir / "thermodynamic_analysis.png", dpi=300, bbox_inches="tight"
        )
        plt.close()

        # Create heatmap of ΔG values for top guides
        if len(results) > 0:
            fig, ax = plt.subplots(figsize=(10, 8))

            # Get top 20 guides
            top_20_indices = np.argsort(combined_scores)[-20:]

            # Create matrix of ΔG values
            dg_matrix = []
            guide_labels = []

            for idx in top_20_indices:
                result = results[idx]
                guide_labels.append(result["guide_id"])

                # Get first 10 target ΔG values
                dg_values = []
                for profile in result["thermodynamics"]["profiles"][:10]:
                    dg_values.append(profile["delta_g"])

                # Pad with NaN if fewer than 10 targets
                while len(dg_values) < 10:
                    dg_values.append(np.nan)

                dg_matrix.append(dg_values)

            # Create heatmap
            sns.heatmap(
                dg_matrix,
                xticklabels=[f"Target {i+1}" for i in range(10)],
                yticklabels=guide_labels,
                cmap="RdYlBu_r",
                center=0,
                cbar_kws={"label": "ΔG (kcal/mol)"},
                fmt=".1f",
                linewidths=0.5,
            )

            plt.title("Free Energy Heatmap (Top 20 Guides)")
            plt.tight_layout()
            plt.savefig(
                output_dir / "delta_g_heatmap.png", dpi=300, bbox_inches="tight"
            )
            plt.close()

    def generate_report(self, results: List[Dict], output_file: Path):
        """Generate thermodynamic validation report."""
        self.console.print(
            "\n[bold green]Thermodynamic Validation Report[/bold green]\n"
        )

        # Summary table
        table = Table(title="Guide RNA Thermodynamic Analysis")
        table.add_column("Guide ID", style="cyan")
        table.add_column("Targets", justify="right")
        table.add_column("Tm (°C)", justify="right", style="yellow")
        table.add_column("ΔG\n(kcal/mol)", justify="right", style="green")
        table.add_column("Binding\nEfficiency", justify="right")
        table.add_column("Specificity", justify="right", style="blue")
        table.add_column("Stability", justify="right", style="magenta")

        for result in results[:20]:  # Show top 20
            thermo = result["thermodynamics"]
            table.add_row(
                result["guide_id"],
                str(result["total_targets"]),
                f"{thermo['average_tm']:.1f}",
                f"{thermo['average_delta_g']:.1f}",
                f"{thermo['average_binding_efficiency']:.2f}",
                f"{thermo['specificity_score']:.2f}",
                f"{thermo['average_stability_score']:.2f}",
            )

        self.console.print(table)

        # Save JSON report
        report = {
            "analysis_type": "thermodynamic",
            "conditions": {
                "temperature": self.thermo_calc.temperature,
                "na_concentration": self.thermo_calc.na_conc,
                "mg_concentration": self.thermo_calc.mg_conc,
            },
            "summary": {
                "total_guides": len(results),
                "total_targets": len(self.targets),
                "average_tm": np.mean(
                    [r["thermodynamics"]["average_tm"] for r in results]
                ),
                "average_delta_g": np.mean(
                    [r["thermodynamics"]["average_delta_g"] for r in results]
                ),
                "average_specificity": np.mean(
                    [r["thermodynamics"]["specificity_score"] for r in results]
                ),
            },
            "guides": results,
        }

        with open(output_file, "w") as f:
            json.dump(report, f, indent=2)

        self.console.print(f"\n[green]Report saved to:[/green] {output_file}")

        # Generate plots
        plot_dir = output_file.parent / "thermo_plots"
        self.generate_thermo_plots(results, plot_dir)
        self.console.print(f"[green]Plots saved to:[/green] {plot_dir}")

        # Print insights
        self.console.print("\n[bold]Thermodynamic Insights:[/bold]")

        # Find best guides by different criteria
        best_tm = max(results, key=lambda x: x["thermodynamics"]["average_tm"])
        best_dg = min(results, key=lambda x: x["thermodynamics"]["average_delta_g"])
        best_spec = max(results, key=lambda x: x["thermodynamics"]["specificity_score"])

        self.console.print(
            f"• Highest Tm: {best_tm['guide_id']} ({best_tm['thermodynamics']['average_tm']:.1f}°C)"
        )
        self.console.print(
            f"• Most stable (lowest ΔG): {best_dg['guide_id']} ({best_dg['thermodynamics']['average_delta_g']:.1f} kcal/mol)"
        )
        self.console.print(
            f"• Most specific: {best_spec['guide_id']} (score: {best_spec['thermodynamics']['specificity_score']:.2f})"
        )


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description="Thermodynamic validation for CRISPR guide RNAs"
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
        default=Path("thermo_validation.json"),
        help="Output JSON report file",
    )
    parser.add_argument(
        "--na-conc",
        type=float,
        default=0.150,
        help="Sodium concentration (M) - Default 150 mM for CRISPR-Cas9",
    )
    parser.add_argument(
        "--mg-conc",
        type=float,
        default=0.010,
        help="Magnesium concentration (M) - Default 10 mM for CRISPR-Cas9",
    )
    parser.add_argument(
        "--temperature", type=float, default=37.0, help="Temperature (°C)"
    )
    parser.add_argument(
        "--required-targets", nargs="+", help="Sequence IDs that must be targeted"
    )

    args = parser.parse_args()

    # Validate files
    for f in [args.guide_csv, args.guide_fasta, args.target_fasta]:
        if not f.exists():
            print(f"Error: File not found: {f}")
            return 1

    # Run validation
    validator = ThermoGuideValidator(
        args.guide_csv,
        args.guide_fasta,
        args.target_fasta,
        args.required_targets,
        args.na_conc,
        args.mg_conc,
        args.temperature,
    )

    results = []
    for guide_id, guide_seq in track(
        validator.guides.items(), description="Analyzing thermodynamics..."
    ):
        results.append(validator.validate_guide_thermo(guide_id, guide_seq))

    validator.generate_report(results, args.output)
    return 0


if __name__ == "__main__":
    exit(main())
