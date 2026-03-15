#!/usr/bin/env python3
"""
Generate detailed stats for individual taxon FNA files in the output directory.
Shows sequence counts, lengths, and genus breakdown for each target.
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PROJECT_DIR = Path(__file__).parent.parent.absolute()
PRIORITY_CSV = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "18s" / "nf_abundance_priority_scores.csv"
OUTPUT_DIR = PROJECT_DIR / "output"
STATS_DIR = PROJECT_DIR / "output_stats"


def get_taxa_info():
    """Load priority data with full lineage info."""
    taxa_info = {}
    with open(PRIORITY_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            taxa_info[row['node']] = {
                'rank': row['rank'],
                'path': row['path'],
                'priority': float(row['total_priority_score']),
                'otu': int(row['total_otu']),
                'genera': [g.strip() for g in row.get('all_genera', '').split(';') if g.strip()],
            }
    return taxa_info


def analyze_fna_file(fna_file):
    """Analyze sequences in a FASTA file. Returns list of sequence lengths."""
    lengths = []
    current_seq = []
    
    with open(fna_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    lengths.append(len(''.join(current_seq)))
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_seq:
            lengths.append(len(''.join(current_seq)))
    
    return lengths


def main():
    print("=" * 80)
    print("📊 OUTPUT FNA STATISTICS GENERATOR")
    print("=" * 80)
    
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load taxa info
    taxa_info = get_taxa_info()
    
    # Process each FNA file
    fna_files = sorted(OUTPUT_DIR.glob("*.fna"))
    print(f"\n🔍 Processing {len(fna_files)} FNA files...\n")
    
    results = []
    total_seqs = 0
    
    for fna_file in fna_files:
        # Parse filename like "001_Euamoebida_order.fna"
        parts = fna_file.stem.split('_')
        if len(parts) >= 3:
            taxon_name = '_'.join(parts[1:-1])
            taxon_rank = parts[-1]
        else:
            taxon_name = fna_file.stem
            taxon_rank = 'unknown'
        
        lengths = analyze_fna_file(fna_file)
        
        if not lengths:
            continue
        
        # Get taxa info
        info = taxa_info.get(taxon_name, {})
        
        result = {
            'file': fna_file.name,
            'taxon': taxon_name,
            'rank': taxon_rank,
            'seq_count': len(lengths),
            'avg_length': round(sum(lengths) / len(lengths), 1),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'priority': info.get('priority', 0),
            'expected_otu': info.get('otu', 0),
            'lineage': info.get('path', 'N/A'),
            'genera_count': len(info.get('genera', [])),
        }
        results.append(result)
        total_seqs += len(lengths)
        
        print(f"  {fna_file.name}: {len(lengths)} seqs, avg {result['avg_length']:.0f} bp")
    
    # Write summary CSV
    summary_file = STATS_DIR / "output_summary.csv"
    with open(summary_file, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
    
    # Print summary
    print(f"\n{'='*80}")
    print("📊 SUMMARY")
    print(f"{'='*80}")
    print(f"   Total FNA files: {len(results)}")
    print(f"   Total sequences: {total_seqs:,}")
    if results:
        avg_len = sum(r['avg_length'] * r['seq_count'] for r in results) / total_seqs
        print(f"   Average length: {avg_len:.0f} bp")
    print(f"\n📁 Stats written to: {summary_file}")


if __name__ == "__main__":
    main()

