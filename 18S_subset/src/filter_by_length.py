#!/usr/bin/env python3
"""
Filter FASTA sequences by minimum length.
Keeps only sequences >= specified length (default 1200bp for full 18S).
"""

import argparse
from pathlib import Path


def parse_fasta(filepath):
    """Parse FASTA file and yield (header, sequence) tuples."""
    with open(filepath, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, ''.join(seq_lines)


def filter_fasta_by_length(input_path, output_path, min_length=1200):
    """Filter FASTA file keeping only sequences >= min_length."""
    kept = 0
    removed = 0
    
    with open(output_path, 'w') as out:
        for header, seq in parse_fasta(input_path):
            if len(seq) >= min_length:
                out.write(f"{header}\n{seq}\n")
                kept += 1
            else:
                removed += 1
    
    return kept, removed


def main():
    parser = argparse.ArgumentParser(description="Filter FASTA by sequence length")
    parser.add_argument("--input-dir", type=str, default="output", help="Input directory with FNA files")
    parser.add_argument("--output-dir", type=str, default="output_filtered", help="Output directory")
    parser.add_argument("--min-length", type=int, default=1200, help="Minimum sequence length (default: 1200)")
    parser.add_argument("--in-place", action="store_true", help="Overwrite original files")
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    
    if args.in_place:
        output_dir = input_dir
    else:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    fna_files = sorted(input_dir.glob("*.fna"))
    
    print(f"=" * 70)
    print(f"🧬 FILTERING SEQUENCES BY LENGTH (>= {args.min_length} bp)")
    print(f"=" * 70)
    print(f"\n📁 Found {len(fna_files)} FNA files\n")
    
    total_kept = 0
    total_removed = 0
    
    for fna_file in fna_files:
        if args.in_place:
            # Write to temp file, then replace
            temp_path = fna_file.with_suffix('.fna.tmp')
            kept, removed = filter_fasta_by_length(fna_file, temp_path, args.min_length)
            temp_path.replace(fna_file)
        else:
            output_path = output_dir / fna_file.name
            kept, removed = filter_fasta_by_length(fna_file, output_path, args.min_length)
        
        total_kept += kept
        total_removed += removed
        
        pct = (kept / (kept + removed) * 100) if (kept + removed) > 0 else 0
        print(f"  {fna_file.name}: {kept} kept, {removed} removed ({pct:.1f}% retained)")
    
    print(f"\n" + "=" * 70)
    print(f"📊 SUMMARY")
    print(f"=" * 70)
    print(f"   Total sequences kept: {total_kept}")
    print(f"   Total sequences removed: {total_removed}")
    pct = (total_kept / (total_kept + total_removed) * 100) if (total_kept + total_removed) > 0 else 0
    print(f"   Retention rate: {pct:.1f}%")
    if not args.in_place:
        print(f"\n📁 Filtered files saved to: {output_dir.absolute()}")


if __name__ == "__main__":
    main()

