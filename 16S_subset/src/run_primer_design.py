#!/usr/bin/env python3
"""
Batch Primer Design Pipeline for rRNA Targets

Runs ppdesign primer on each extracted FNA file to design PCR primers
for amplifying conserved regions within each taxonomic group.

Supports two input types:
  - output: Individual taxon FNA files
  - umbrella: Combined umbrella FNA files

Usage:
    python src/run_primer_design.py --input-type output --limit 10
    python src/run_primer_design.py --input-type umbrella --limit 10
"""

import subprocess
import argparse
from pathlib import Path

# Directory paths
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_DIR = SCRIPT_DIR.parent

# ppdesign virtual environment path
PPDESIGN_VENV = PROJECT_DIR.parent / "ppdesign" / ".venv" / "bin" / "ppdesign"


def run_primer_design(
    fna_file: Path,
    output_dir: Path,
    output_name: str,
    conservation: float = 0.6,
    method: str = "kmer",
    amplicon_min: int = 200,
    amplicon_max: int = 800,
    threads: int = 4,
    mafft_auto: bool = False,
    structure_aware: bool = False,
) -> dict:
    """Run ppdesign primer on a single FNA file."""

    result_dir = (output_dir / output_name).absolute()

    cmd = [
        str(PPDESIGN_VENV), "primer", "main",
        "--fasta-input", str(fna_file.absolute()),
        "--output-dir", str(result_dir),
        "--conservation", str(conservation),
        "--amplicon-min", str(amplicon_min),
        "--amplicon-max", str(amplicon_max),
        "--threads", str(threads),
        "--max-degenerate-positions", "2",
        "--method", method,
    ]

    if structure_aware and method == "msa":
        cmd.append("--structure-aware")
    elif mafft_auto and method == "msa":
        cmd.append("--mafft-auto")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    # Check for primer_pairs.csv to determine success
    primer_csv = result_dir / "primer_pairs.csv"
    success = primer_csv.exists()

    num_pairs = 0
    if success:
        with open(primer_csv, 'r') as f:
            num_pairs = sum(1 for _ in f) - 1  # Subtract header

    return {
        "file": fna_file.name,
        "output": output_name,
        "success": success,
        "num_pairs": num_pairs,
        "returncode": result.returncode,
        "stderr": result.stderr[:500] if result.stderr else ""
    }


def main():
    parser = argparse.ArgumentParser(description="Batch primer design for rRNA targets")
    parser.add_argument("--input-type", type=str, default="output", choices=["output", "umbrella"],
                        help="Input type: 'output' for individual taxon files, 'umbrella' for combined files")
    parser.add_argument("--conservation", type=float, default=0.3, help="Min conservation (0-1)")
    parser.add_argument("--method", type=str, default="kmer", choices=["kmer", "msa"],
                        help="Conserved region finding method")
    parser.add_argument("--mafft-auto", action="store_true",
                        help="Use MAFFT --auto (faster than L-INS-i, good for >100 seqs)")
    parser.add_argument("--structure-aware", action="store_true",
                        help="Use Q-INS-i structure-aware alignment for rRNA")
    parser.add_argument("--amplicon-min", type=int, default=200, help="Min amplicon size")
    parser.add_argument("--amplicon-max", type=int, default=800, help="Max amplicon size")
    parser.add_argument("--threads", type=int, default=4, help="Threads per job")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of files to process")
    args = parser.parse_args()

    # Set input/output directories based on input type
    if args.input_type == "umbrella":
        input_dir = PROJECT_DIR / "umbrella_output"
        results_dir = PROJECT_DIR / "umbrella_primer_results"
    else:
        input_dir = PROJECT_DIR / "output"
        results_dir = PROJECT_DIR / "primer_results"

    mafft_auto = getattr(args, 'mafft_auto', False)
    structure_aware = getattr(args, 'structure_aware', False)

    # Determine alignment mode description
    if structure_aware:
        mode_desc = "msa (Q-INS-i structure-aware)"
    elif mafft_auto:
        mode_desc = "msa (mafft-auto)"
    else:
        mode_desc = args.method

    # Detect project type from directory name
    project_name = PROJECT_DIR.name.upper().replace("_SUBSET", "")

    print("=" * 70)
    print(f"🧬 BATCH PRIMER DESIGN FOR {project_name} ({args.input_type.upper()})")
    print(f"   Input: {input_dir}")
    print(f"   Conservation: {args.conservation}")
    print(f"   Amplicon size: {args.amplicon_min}-{args.amplicon_max} bp")
    print(f"   Method: {mode_desc}")
    print("=" * 70)

    results_dir.mkdir(parents=True, exist_ok=True)

    # Get all FNA files sorted by priority (numeric prefix)
    fna_files = sorted(input_dir.glob("*.fna"))
    if args.limit:
        fna_files = fna_files[:args.limit]

    print(f"\n📁 Found {len(fna_files)} FNA files to process\n")

    results = []
    for i, fna_file in enumerate(fna_files, 1):
        output_name = fna_file.stem

        print(f"[{i:3d}/{len(fna_files)}] {fna_file.name}")

        result = run_primer_design(
            fna_file,
            results_dir,
            output_name,
            conservation=args.conservation,
            method=args.method,
            amplicon_min=args.amplicon_min,
            amplicon_max=args.amplicon_max,
            threads=args.threads,
            mafft_auto=mafft_auto,
            structure_aware=structure_aware,
        )
        results.append(result)

        status = "✓" if result["success"] else "✗"
        print(f"         {status} {result['num_pairs']} primer pairs found")

    # Summary
    print("\n" + "=" * 70)
    print("📊 SUMMARY")
    print("=" * 70)
    successful = [r for r in results if r["success"] and r["num_pairs"] > 0]
    print(f"   Total processed: {len(results)}")
    print(f"   With primers: {len(successful)}")
    print(f"   No primers: {len(results) - len(successful)}")

    if successful:
        total_pairs = sum(r["num_pairs"] for r in successful)
        print(f"   Total primer pairs: {total_pairs}")
        print(f"\n   Top targets by primer count:")
        for r in sorted(successful, key=lambda x: x["num_pairs"], reverse=True)[:10]:
            print(f"      {r['output']}: {r['num_pairs']} pairs")

    print(f"\n📁 Results saved to: {results_dir}")


if __name__ == "__main__":
    main()

