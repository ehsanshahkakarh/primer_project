#!/usr/bin/env python3
"""
Batch Primer Design Pipeline for 18S Division Targets

Runs ppdesign primer on each extracted FNA file to design PCR primers
for amplifying conserved regions within each taxonomic division.

Usage:
    python src/run_primer_design.py [--structure-aware] [--conservation 0.6]
"""

import subprocess
import sys
import argparse
from pathlib import Path
from datetime import datetime

# Directory paths
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_DIR = SCRIPT_DIR.parent
INPUT_DIR = PROJECT_DIR / "output"
RESULTS_DIR = PROJECT_DIR / "primer_results"


def run_primer_design(
    fna_file: Path,
    output_name: str,
    conservation: float = 0.6,
    method: str = "kmer",  # Use kmer by default for speed
    amplicon_min: int = 200,
    amplicon_max: int = 800,
    threads: int = 4,
) -> dict:
    """Run ppdesign primer on a single FNA file."""

    # Use absolute path so ppdesign doesn't prepend 'results/'
    output_dir = (RESULTS_DIR / output_name).absolute()

    cmd = [
        "ppdesign", "primer", "main",
        "--fasta-input", str(fna_file.absolute()),
        "--output-dir", str(output_dir),
        "--conservation", str(conservation),
        "--amplicon-min", str(amplicon_min),
        "--amplicon-max", str(amplicon_max),
        "--threads", str(threads),
        "--max-degenerate-positions", "2",
        "--method", method,
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    # Check for primer_pairs.csv to determine success
    primer_csv = output_dir / "primer_pairs.csv"
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
    parser = argparse.ArgumentParser(description="Batch primer design for 18S division targets")
    parser.add_argument("--conservation", type=float, default=0.3, help="Min conservation (0-1)")
    parser.add_argument("--method", type=str, default="msa", choices=["kmer", "msa"], help="Conserved region finding method (msa recommended for rRNA)")
    parser.add_argument("--amplicon-min", type=int, default=200, help="Min amplicon size")
    parser.add_argument("--amplicon-max", type=int, default=800, help="Max amplicon size")
    parser.add_argument("--threads", type=int, default=4, help="Threads per job")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of files to process")
    args = parser.parse_args()

    print("=" * 70)
    print("🧬 BATCH PRIMER DESIGN FOR 18S DIVISION TARGETS")
    print(f"   Conservation: {args.conservation}")
    print(f"   Amplicon size: {args.amplicon_min}-{args.amplicon_max} bp")
    print(f"   Method: {args.method}")
    print("=" * 70)
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Get all FNA files sorted by priority (numeric prefix)
    fna_files = sorted(INPUT_DIR.glob("*.fna"))
    if args.limit:
        fna_files = fna_files[:args.limit]
    
    print(f"\n📁 Found {len(fna_files)} FNA files to process\n")
    
    results = []
    for i, fna_file in enumerate(fna_files, 1):
        # Extract name without extension and prefix number
        output_name = fna_file.stem  # e.g., "01_Arthropoda_division"
        
        print(f"[{i:3d}/{len(fna_files)}] {fna_file.name}")
        
        result = run_primer_design(
            fna_file,
            output_name,
            conservation=args.conservation,
            method=args.method,
            amplicon_min=args.amplicon_min,
            amplicon_max=args.amplicon_max,
            threads=args.threads,
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
    
    print(f"\n📁 Results saved to: {RESULTS_DIR}")


if __name__ == "__main__":
    main()

