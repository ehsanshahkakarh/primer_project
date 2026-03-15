#!/usr/bin/env python3
"""
Generate detailed stats for each umbrella FASTA file (16S version).
Shows sequences grouped by source taxon with counts and avg lengths.
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PROJECT_DIR = Path(__file__).parent.parent.absolute()
PRIORITY_DIR = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "16s"
INPUT_DIR = PROJECT_DIR / "output"
UMBRELLA_DIR = PROJECT_DIR / "umbrella_output"
STATS_DIR = PROJECT_DIR / "umbrella_stats"


def load_all_priority_data():
    """Load priority data from both Archaea and Bacteria files."""
    taxa_info = {}
    for domain in ['archaea', 'bacteria']:
        priority_file = PRIORITY_DIR / domain / "nf_abundance_priority_scores.csv"
        if not priority_file.exists():
            continue
        with open(priority_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                taxa_info[row['node']] = {
                    'domain': domain,
                    'rank': row['rank'],
                    'path': row['path'],
                    'priority': float(row['total_priority_score']),
                    'otu': int(row['total_otu']),
                    'genera': [g.strip() for g in row.get('all_genera', '').split(';') if g.strip()],
                }
    return taxa_info


def analyze_umbrella_file(fna_file):
    """Analyze sequences in an umbrella FASTA file."""
    umbrella_name = fna_file.stem.split('_', 1)[1].replace('_umbrella', '')
    all_seqs = []
    current_header = None
    current_seq = []
    
    with open(fna_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header and current_seq:
                    seq_len = len(''.join(current_seq))
                    all_seqs.append({'header': current_header, 'length': seq_len})
                current_header = line.strip()[1:]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_header and current_seq:
            seq_len = len(''.join(current_seq))
            all_seqs.append({'header': current_header, 'length': seq_len})
    
    return umbrella_name, all_seqs


def match_seqs_to_source_files(input_dir):
    """Match umbrella sequences back to their source FNA files."""
    source_index = {}
    for fna_file in input_dir.glob("*.fna"):
        parts = fna_file.stem.split('_')
        if len(parts) >= 3:
            taxon_name = '_'.join(parts[1:-1])
            taxon_rank = parts[-1]
        else:
            taxon_name = fna_file.stem
            taxon_rank = 'unknown'
        
        with open(fna_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]
                    source_index[header] = (taxon_name, taxon_rank, fna_file.name)
    return source_index


def main():
    print("=" * 80)
    print("📊 16S UMBRELLA FASTA STATISTICS GENERATOR")
    print("=" * 80)
    
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    
    print("\n📁 Building sequence source index...")
    source_index = match_seqs_to_source_files(INPUT_DIR)
    print(f"   Indexed {len(source_index):,} sequences")
    
    taxa_info = load_all_priority_data()
    umbrella_files = sorted(UMBRELLA_DIR.glob("*_umbrella.fna"))
    print(f"\n🔍 Processing {len(umbrella_files)} umbrella files...\n")
    
    for fna_file in umbrella_files:
        umbrella_name, all_seqs = analyze_umbrella_file(fna_file)
        by_taxon = defaultdict(list)
        
        for seq in all_seqs:
            clean_header = seq['header'].split('|umbrella=')[0]
            if clean_header in source_index:
                taxon_name, taxon_rank, _ = source_index[clean_header]
                by_taxon[(taxon_name, taxon_rank)].append(seq['length'])
        
        stats_file = STATS_DIR / f"{fna_file.stem}_stats.csv"
        sorted_taxa = sorted(by_taxon.items(), key=lambda x: -len(x[1]))
        
        with open(stats_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['taxon', 'rank', 'domain', 'lineage', 'seq_count', 'avg_length', 'min_length', 'max_length'])
            for (taxon_name, taxon_rank), lengths in sorted_taxa:
                info = taxa_info.get(taxon_name, {})
                writer.writerow([
                    taxon_name, taxon_rank, info.get('domain', 'unknown'),
                    info.get('path', 'N/A'), len(lengths),
                    round(sum(lengths)/len(lengths), 1), min(lengths), max(lengths)
                ])
        
        total_seqs = len(all_seqs)
        avg_len = sum(s['length'] for s in all_seqs) / total_seqs if total_seqs else 0
        print(f"📄 {fna_file.name}: {total_seqs:,} seqs | Avg: {avg_len:.0f} bp | Taxa: {len(by_taxon)}")

    print(f"\n{'='*80}")
    print(f"✅ Stats written to: {STATS_DIR}")
    print("=" * 80)


if __name__ == "__main__":
    main()

