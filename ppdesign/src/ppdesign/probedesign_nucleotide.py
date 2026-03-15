#!/usr/bin/env python3
"""
Direct nucleotide probe design pipeline without gene calling.
Finds conserved regions directly from nucleotide sequences.
"""

import logging
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import typer
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from .conserved_finder import (
    KmerBasedFinder,
    MSABasedFinder,
    Minimap2BasedFinder,
    ConservedRegion,
)
from .output_utils import resolve_output_dir

app = typer.Typer()

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class NucleotideProbeDesign:
    """Pipeline for designing probes directly from nucleotide sequences."""

    def _check_dependencies(self, method: str):
        """Check external tool availability for specified method."""
        required_tools = {
            "minimap2": ["minimap2"],
            "msa": ["mafft"]
        }

        if method in required_tools:
            missing = [tool for tool in required_tools[method] if shutil.which(tool) is None]
            if missing:
                raise RuntimeError(
                    f"Missing required tools for {method} method: {', '.join(missing)}"
                )

    def __init__(self, output_dir: str, threads: int = 4):
        """Initialize the pipeline."""
        self.output_dir = resolve_output_dir(output_dir)
        self.threads = threads
        self.logger = logging.getLogger(__name__)

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_sequences(self, fasta_input: str) -> Dict[str, str]:
        """Load sequences from a FASTA file or directory."""
        fasta_path = Path(fasta_input)

        if not fasta_path.exists():
            raise ValueError(f"FASTA path not found: {fasta_input}")

        # Gather FASTA files
        if fasta_path.is_file():
            fasta_files = [fasta_path]
        else:
            fasta_files = (
                list(fasta_path.glob("*.fna"))
                + list(fasta_path.glob("*.fasta"))
                + list(fasta_path.glob("*.fa"))
            )

            if not fasta_files:
                raise ValueError(f"No FASTA files found in: {fasta_input}")

        self.logger.info(f"Found {len(fasta_files)} FASTA files")

        # Load sequences
        sequences = {}
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = f"{fasta_file.stem}_{record.id}"
                sequences[seq_id] = str(record.seq).upper()

        self.logger.info(f"Loaded {len(sequences)} sequences")
        return sequences

    def find_conserved_regions(
        self,
        sequences: Dict[str, str],
        method: str = "kmer",
        min_length: int = 15,
        max_length: int = 30,
        min_conservation: float = 0.8,
        kmer_size: int = 20,
        window_size: Optional[int] = None,
        step_size: Optional[int] = None,
    ) -> List[ConservedRegion]:
        """Find conserved regions using specified method."""

        # Check dependencies before starting work
        self._check_dependencies(method)

        if method == "kmer":
            self.logger.info(f"Using k-mer based method (k={kmer_size})")
            finder = KmerBasedFinder(
                kmer_size=kmer_size, min_conservation=min_conservation
            )
            regions = finder.find_conserved_regions(
                sequences=sequences,
                min_length=min_length,
                max_length=max_length,
                allow_mismatches=2,
            )

        elif method == "msa":
            self.logger.info("Using MSA based method")
            finder = MSABasedFinder(min_conservation=min_conservation)
            regions = finder.find_conserved_regions(
                sequences=sequences,
                min_length=min_length,
                max_length=max_length,
                window_size=window_size,
                step_size=step_size,
            )

        elif method == "minimap2":
            self.logger.info("Using minimap2 based method")
            finder = Minimap2BasedFinder(
                min_conservation=min_conservation,
                identity_threshold=0.9,  # Default 90% identity for probes
            )
            regions = finder.find_conserved_regions(
                sequences=sequences,
                min_length=min_length,
                max_length=max_length,
                preset="asm20",  # Default preset for <20% divergence
            )

        else:
            raise ValueError(
                f"Unknown method: {method}. Choose from: kmer, msa, minimap2"
            )

        self.logger.info(f"Found {len(regions)} conserved regions")
        return regions

    def filter_oligos(
        self,
        regions: List[ConservedRegion],
        gc_range: Tuple[float, float] = (30, 70),
        tm_range: Tuple[float, float] = (40, 60),
        check_hairpins: bool = True,
        check_dimers: bool = True,
    ) -> pd.DataFrame:
        """Filter conserved regions to create final oligo list."""

        oligos = []

        for i, region in enumerate(regions):
            # Calculate properties for consensus sequence
            gc_content = self._calculate_gc(region.consensus)
            tm = self._calculate_tm(region.consensus)
            nd_perc = self._calculate_nd_perc(region.consensus)

            # Check if within ranges
            if not (gc_range[0] <= gc_content <= gc_range[1]):
                continue

            if not (tm_range[0] <= tm <= tm_range[1]):
                continue

            # Check secondary structures
            has_hairpin = False
            has_dimer = False

            if check_hairpins:
                has_hairpin = self._check_hairpin(region.consensus)

            if check_dimers:
                has_dimer = self._check_dimer(region.consensus)

            # Calculate quality score
            quality_score = self._calculate_quality_score(
                conservation=region.conservation,
                gc_content=gc_content,
                tm=tm,
                nd_perc=nd_perc,
                gc_range=gc_range,
                tm_range=tm_range,
            )

            # Create oligo entry
            oligo = {
                "Oligo_ID": f"oligo_{i+1}",
                "Sequence": region.consensus,
                "Length": len(region.consensus),
                "GC_Content": round(gc_content, 2),
                "Tm": round(tm, 2),
                "Conservation": round(region.conservation, 3),
                "ND_Percentage": round(nd_perc, 2),
                "Has_Hairpin": has_hairpin,
                "Has_Dimer": has_dimer,
                "Quality_Score": round(quality_score, 3),
                "Source_Sequences": len(region.sequences),
                "Region_Start": region.start,
                "Region_End": region.end,
            }

            oligos.append(oligo)

        # Convert to DataFrame and sort by quality score
        oligos_df = pd.DataFrame(oligos)

        if not oligos_df.empty:
            oligos_df = oligos_df.sort_values("Quality_Score", ascending=False)
            oligos_df = oligos_df.reset_index(drop=True)

        self.logger.info(f"Filtered to {len(oligos_df)} oligos")

        return oligos_df

    def _calculate_gc(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        gc_count = sequence.count("G") + sequence.count("C")
        # Count degenerate bases that contribute to GC
        gc_count += sequence.count("S")  # G or C
        gc_count += sequence.count("V") * 2 / 3  # A, C, or G (2/3 are G or C)
        gc_count += sequence.count("B") * 2 / 3  # C, G, or T
        gc_count += sequence.count("D") * 1 / 3  # A, G, or T
        gc_count += sequence.count("H") * 1 / 3  # A, C, or T
        gc_count += sequence.count("N") * 0.5  # Any base

        return (gc_count / len(sequence)) * 100 if sequence else 0

    def _calculate_tm(self, sequence: str) -> float:
        """Calculate melting temperature."""
        # Convert degenerate bases to most common base for Tm calculation
        clean_seq = sequence
        for degen, replacement in [
            ("R", "A"),
            ("Y", "C"),
            ("S", "G"),
            ("W", "A"),
            ("K", "G"),
            ("M", "A"),
            ("B", "C"),
            ("D", "A"),
            ("H", "A"),
            ("V", "A"),
            ("N", "A"),
        ]:
            clean_seq = clean_seq.replace(degen, replacement)

        try:
            tm = mt.Tm_Wallace(Seq(clean_seq), strict=False)
        except Exception:
            # Fallback to simple calculation
            tm = 4 * (clean_seq.count("G") + clean_seq.count("C")) + 2 * (
                clean_seq.count("A") + clean_seq.count("T")
            )

        return tm

    def _calculate_nd_perc(self, sequence: str) -> float:
        """Calculate percentage of non-degenerate bases."""
        non_degenerate = "ATCG"
        nd_count = sum(1 for base in sequence if base in non_degenerate)
        return (nd_count / len(sequence)) * 100 if sequence else 0

    def _check_hairpin(
        self, sequence: str, min_stem: int = 4, min_loop: int = 3
    ) -> bool:
        """Check if sequence can form hairpin."""
        seq_len = len(sequence)

        # Check for potential hairpin formations
        for i in range(seq_len - min_stem * 2 - min_loop):
            for j in range(i + min_stem + min_loop, seq_len - min_stem + 1):
                # Check if stem regions are complementary
                stem1 = sequence[i : i + min_stem]
                stem2 = sequence[j : j + min_stem]

                if self._is_complementary(stem1, stem2[::-1]):
                    return True

        return False

    def _check_dimer(self, sequence: str, min_overlap: int = 6) -> bool:
        """Check if sequence can form dimers."""
        rev_comp = self._reverse_complement(sequence)

        # Check for overlaps between sequence and its reverse complement
        for i in range(len(sequence) - min_overlap + 1):
            for j in range(len(rev_comp) - min_overlap + 1):
                if sequence[i : i + min_overlap] == rev_comp[j : j + min_overlap]:
                    return True

        return False

    def _is_complementary(self, seq1: str, seq2: str) -> bool:
        """Check if two sequences are complementary."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

        if len(seq1) != len(seq2):
            return False

        matches = 0
        for b1, b2 in zip(seq1, seq2):
            if b1 in complement and complement[b1] == b2:
                matches += 1

        # Allow some mismatches
        return matches >= len(seq1) * 0.8

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of sequence."""
        complement = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
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
            "N": "N",
        }

        rev_comp = []
        for base in reversed(sequence):
            rev_comp.append(complement.get(base, "N"))

        return "".join(rev_comp)

    def _calculate_quality_score(
        self,
        conservation: float,
        gc_content: float,
        tm: float,
        nd_perc: float,
        gc_range: Tuple[float, float],
        tm_range: Tuple[float, float],
    ) -> float:
        """Calculate overall quality score for oligo."""
        # Normalize values to 0-1 range
        conservation_score = conservation

        # GC score (closer to middle of range is better)
        gc_mid = (gc_range[0] + gc_range[1]) / 2
        gc_width = (gc_range[1] - gc_range[0]) / 2
        gc_score = 1 - abs(gc_content - gc_mid) / gc_width
        gc_score = max(0, min(1, gc_score))

        # Tm score (closer to middle of range is better)
        tm_mid = (tm_range[0] + tm_range[1]) / 2
        tm_width = (tm_range[1] - tm_range[0]) / 2
        tm_score = 1 - abs(tm - tm_mid) / tm_width
        tm_score = max(0, min(1, tm_score))

        # ND percentage score
        nd_score = nd_perc / 100

        # Weighted average
        weights = {"conservation": 0.3, "gc": 0.2, "tm": 0.2, "nd": 0.3}

        quality_score = (
            weights["conservation"] * conservation_score
            + weights["gc"] * gc_score
            + weights["tm"] * tm_score
            + weights["nd"] * nd_score
        )

        return quality_score

    def save_results(
        self, oligos_df: pd.DataFrame, conserved_regions: List[ConservedRegion]
    ):
        """Save results to files."""
        # Save oligos
        oligos_file = self.output_dir / "selected_oligos.csv"
        oligos_df.to_csv(oligos_file, index=False)
        self.logger.info(f"Saved {len(oligos_df)} oligos to {oligos_file}")

        # Save conserved regions details
        regions_file = self.output_dir / "conserved_regions.txt"
        with open(regions_file, "w") as f:
            f.write("# Conserved Regions Found\n\n")

            for i, region in enumerate(conserved_regions):
                f.write(f"Region {i+1}:\n")
                f.write(f"  Consensus: {region.consensus}\n")
                f.write(f"  Length: {len(region.consensus)}\n")
                f.write(f"  Conservation: {region.conservation:.3f}\n")
                f.write(f"  Found in {len(region.sequences)} sequences\n")
                f.write("  Positions:\n")

                for seq_id, pos in list(region.positions.items())[:5]:  # Show first 5
                    f.write(f"    {seq_id}: {pos}\n")

                if len(region.positions) > 5:
                    f.write(f"    ... and {len(region.positions) - 5} more\n")

                f.write("\n")

        self.logger.info(f"Saved conserved regions details to {regions_file}")

        # Save summary statistics
        stats_file = self.output_dir / "pipeline_stats.txt"
        with open(stats_file, "w") as f:
            f.write("# Nucleotide Probe Design Pipeline Statistics\n\n")
            f.write(f"Total conserved regions found: {len(conserved_regions)}\n")
            f.write(f"Regions passing filters: {len(oligos_df)}\n")

            if not oligos_df.empty:
                f.write("\nOligo Statistics:\n")
                f.write(
                    f"  Length range: {oligos_df['Length'].min()}-{oligos_df['Length'].max()}\n"
                )
                f.write(
                    f"  GC% range: {oligos_df['GC_Content'].min():.1f}-{oligos_df['GC_Content'].max():.1f}\n"
                )
                f.write(
                    f"  Tm range: {oligos_df['Tm'].min():.1f}-{oligos_df['Tm'].max():.1f}\n"
                )
                f.write(
                    f"  Conservation range: {oligos_df['Conservation'].min():.3f}-{oligos_df['Conservation'].max():.3f}\n"
                )
                f.write(f"  With hairpins: {oligos_df['Has_Hairpin'].sum()}\n")
                f.write(f"  With dimers: {oligos_df['Has_Dimer'].sum()}\n")


@app.command()
def main(
    fasta_input: Optional[str] = typer.Option(
        None,
        "--fasta-input",
        help="Path to FASTA file or directory containing nucleotide sequences",
    ),
    fasta_dir: Optional[str] = typer.Option(
        None, "--fasta-dir", help="Directory containing nucleotide FASTA files"
    ),
    output_dir: str = typer.Option(
        ..., "-o", "--output-dir", help="Output directory for results"
    ),
    method: str = typer.Option(
        "kmer",
        "--method",
        help="Method for finding conserved regions: kmer, msa, or minimap2",
    ),
    min_length: int = typer.Option(15, "--min-length", help="Minimum oligo length"),
    max_length: int = typer.Option(30, "--max-length", help="Maximum oligo length"),
    gc_min: float = typer.Option(30, "--gc-min", help="Minimum GC percentage"),
    gc_max: float = typer.Option(70, "--gc-max", help="Maximum GC percentage"),
    tm_min: float = typer.Option(40, "--tm-min", help="Minimum melting temperature"),
    tm_max: float = typer.Option(60, "--tm-max", help="Maximum melting temperature"),
    conservation: float = typer.Option(
        0.8, "-c", "--conservation", help="Minimum conservation level"
    ),
    kmer_size: int = typer.Option(20, "--kmer-size", help="K-mer size for kmer method"),
    window_size: Optional[int] = typer.Option(
        None, "--window-size", help="Window size for MSA method"
    ),
    step_size: Optional[int] = typer.Option(
        None, "--step-size", help="Step size for MSA method"
    ),
    threads: int = typer.Option(4, "-j", "--threads", help="Number of threads"),
    no_hairpin_check: bool = typer.Option(
        False, "--no-hairpin-check", help="Skip hairpin checking"
    ),
    no_dimer_check: bool = typer.Option(
        False, "--no-dimer-check", help="Skip dimer checking"
    ),
):
    """
    Design probes directly from nucleotide sequences without gene calling.

    Two methods available:
    - kmer: Fast k-mer based approach (default)
    - msa: Multiple sequence alignment based approach
    """

    # Initialize pipeline
    pipeline = NucleotideProbeDesign(output_dir, threads)

    try:
        # Validate inputs
        if fasta_input is None and fasta_dir is None:
            typer.echo("Error: Must provide either --fasta-input or --fasta-dir", err=True)
            raise typer.Exit(1)

        fasta_path = fasta_input if fasta_input else fasta_dir
        if fasta_input is None and fasta_dir is not None:
            typer.echo(
                "[INFO] --fasta-dir is supported; consider --fasta-input (accepts files or directories).",
                err=True,
            )

        # Load sequences
        sequences = pipeline.load_sequences(fasta_path)

        if len(sequences) < 2:
            typer.echo(
                "Error: Need at least 2 sequences for conservation analysis", err=True
            )
            raise typer.Exit(1)

        # Find conserved regions
        conserved_regions = pipeline.find_conserved_regions(
            sequences=sequences,
            method=method,
            min_length=min_length,
            max_length=max_length,
            min_conservation=conservation,
            kmer_size=kmer_size,
            window_size=window_size,
            step_size=step_size,
        )

        if not conserved_regions:
            typer.echo("No conserved regions found with specified parameters", err=True)
            raise typer.Exit(1)

        # Filter oligos
        oligos_df = pipeline.filter_oligos(
            regions=conserved_regions,
            gc_range=(gc_min, gc_max),
            tm_range=(tm_min, tm_max),
            check_hairpins=not no_hairpin_check,
            check_dimers=not no_dimer_check,
        )

        if oligos_df.empty:
            typer.echo("No oligos passed the filtering criteria", err=True)
            raise typer.Exit(1)

        # Save results
        pipeline.save_results(oligos_df, conserved_regions)

        typer.echo(
            f"Pipeline completed successfully. Found {len(oligos_df)} suitable oligos."
        )
        typer.echo(f"Results saved to: {output_dir}")

    except Exception as e:
        typer.echo(f"Pipeline failed: {e}", err=True)
        raise typer.Exit(1)


if __name__ == "__main__":
    app()


@app.callback(invoke_without_command=True)
def _default(
    ctx: typer.Context,
    fasta_input: Optional[str] = typer.Option(
        None,
        "--fasta-input",
        help="Path to FASTA file or directory containing nucleotide sequences",
    ),
    fasta_dir: Optional[str] = typer.Option(
        None, "--fasta-dir", help="Directory containing nucleotide FASTA files"
    ),
    output_dir: str = typer.Option(
        None, "-o", "--output-dir", help="Output directory for results"
    ),
    method: str = typer.Option(
        "kmer",
        "--method",
        help="Method for finding conserved regions: kmer, msa, or minimap2",
    ),
    min_length: int = typer.Option(15, "--min-length", help="Minimum oligo length"),
    max_length: int = typer.Option(30, "--max-length", help="Maximum oligo length"),
    gc_min: float = typer.Option(30, "--gc-min", help="Minimum GC percentage"),
    gc_max: float = typer.Option(70, "--gc-max", help="Maximum GC percentage"),
    tm_min: float = typer.Option(40, "--tm-min", help="Minimum melting temperature"),
    tm_max: float = typer.Option(60, "--tm-max", help="Maximum melting temperature"),
    conservation: float = typer.Option(
        0.8, "-c", "--conservation", help="Minimum conservation level"
    ),
    kmer_size: int = typer.Option(20, "--kmer-size", help="K-mer size for kmer method"),
    window_size: Optional[int] = typer.Option(
        None, "--window-size", help="Window size for MSA method"
    ),
    step_size: Optional[int] = typer.Option(
        None, "--step-size", help="Step size for MSA method"
    ),
    threads: int = typer.Option(4, "-j", "--threads", help="Number of threads"),
    no_hairpin_check: bool = typer.Option(
        False, "--no-hairpin-check", help="Skip hairpin checking"
    ),
    no_dimer_check: bool = typer.Option(
        False, "--no-dimer-check", help="Skip dimer checking"
    ),
):
    if ctx.invoked_subcommand is None:
        if output_dir is None or (fasta_input is None and fasta_dir is None):
            raise typer.BadParameter("--fasta-input/--fasta-dir and --output-dir are required")
        return main(
            fasta_input=fasta_input,
            fasta_dir=fasta_dir,
            output_dir=output_dir,
            method=method,
            min_length=min_length,
            max_length=max_length,
            gc_min=gc_min,
            gc_max=gc_max,
            tm_min=tm_min,
            tm_max=tm_max,
            conservation=conservation,
            kmer_size=kmer_size,
            window_size=window_size,
            step_size=step_size,
            threads=threads,
            no_hairpin_check=no_hairpin_check,
            no_dimer_check=no_dimer_check,
        )
    return None
