#!/usr/bin/env python3
"""
Filter FASTA files by minimum sequence length.
Removes sequences shorter than the specified threshold.

Usage:
    python filter_by_length.py --input-dir <dir> --min-length 1200 [--in-place]
"""

import argparse
from pathlib import Path


def filter_fasta_by_length(input_file, output_file, min_length):
    """Filter sequences by minimum length."""
    kept = 0
    removed = 0
    current_header = None
    current_seq = []
    output_data = []
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Process previous sequence
                if current_header and current_seq:
                    seq = ''.join(current_seq)
                    if len(seq) >= min_length:
                        output_data.append((current_header, seq))
                        kept += 1
                    else:
                        removed += 1
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # Don't forget last sequence
        if current_header and current_seq:
            seq = ''.join(current_seq)
            if len(seq) >= min_length:
                output_data.append((current_header, seq))
                kept += 1
            else:
                removed += 1
    
    # Write output
    with open(output_file, 'w') as f:
        for header, seq in output_data:
            f.write(f"{header}\n{seq}\n")
    
    return kept, removed


def main():
    parser = argparse.ArgumentParser(description='Filter FASTA files by sequence length')
    parser.add_argument('--input-dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('--output-dir', help='Output directory (default: input-dir_filtered)')
    parser.add_argument('--min-length', type=int, default=1200, help='Minimum sequence length (default: 1200)')
    parser.add_argument('--in-place', action='store_true', help='Modify files in place')
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        return
    
    if args.in_place:
        output_dir = input_dir
    else:
        output_dir = Path(args.output_dir) if args.output_dir else input_dir.parent / f"{input_dir.name}_filtered"
        output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print(f"🔬 FILTERING FASTA FILES BY LENGTH (≥{args.min_length} bp)")
    print("=" * 80)
    print(f"   Input: {input_dir}")
    print(f"   Output: {output_dir}")
    print(f"   Min length: {args.min_length} bp")
    print()
    
    fasta_files = list(input_dir.glob("*.fna")) + list(input_dir.glob("*.fasta"))
    total_kept = 0
    total_removed = 0
    
    for fna_file in sorted(fasta_files):
        output_file = output_dir / fna_file.name
        kept, removed = filter_fasta_by_length(fna_file, output_file, args.min_length)
        total_kept += kept
        total_removed += removed
        print(f"  {fna_file.name}: {kept} kept, {removed} removed")
    
    print()
    print("=" * 80)
    print(f"✅ Filtering complete!")
    print(f"   Total kept: {total_kept:,} sequences")
    print(f"   Total removed: {total_removed:,} sequences")
    print(f"   Retention rate: {100*total_kept/(total_kept+total_removed):.1f}%")
    print("=" * 80)


if __name__ == "__main__":
    main()

