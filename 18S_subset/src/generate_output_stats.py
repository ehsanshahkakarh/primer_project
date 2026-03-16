#!/usr/bin/env python3
"""
Generate detailed stats for individual taxon FNA files in the output directory.
Creates per-taxon stats.csv and genera.csv files, plus an overall summary.
"""

import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)

PROJECT_DIR = Path(__file__).parent.parent.absolute()
PRIORITY_CSV = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "18s" / "nf_abundance_priority_scores.csv"
OUTPUT_DIR = PROJECT_DIR / "output"
STATS_DIR = PROJECT_DIR / "output_stats"


def get_taxa_info():
    """Load priority data with full lineage info and genera with scores."""
    taxa_info = {}
    with open(PRIORITY_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse top_genera to get genus with scores like "Genus1(12.5); Genus2(8.3)"
            top_genera_raw = row.get('top_genera', '')
            genera_with_scores = []
            for match in re.findall(r'(\w+)\(([0-9.]+)\)', top_genera_raw):
                genera_with_scores.append({'name': match[0], 'score': float(match[1])})

            taxa_info[row['node']] = {
                'rank': row['rank'],
                'path': row['path'],
                'priority': float(row['total_priority_score']),
                'otu': int(row['total_otu']),
                'genera': [g.strip() for g in row.get('all_genera', '').split(';') if g.strip()],
                'genera_with_scores': genera_with_scores,
            }
    return taxa_info


def analyze_fna_file(fna_file):
    """Analyze sequences in a FASTA file. Returns list of (header, length) tuples."""
    sequences = []
    current_header = None
    current_seq = []

    with open(fna_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header and current_seq:
                    sequences.append({
                        'header': current_header,
                        'length': len(''.join(current_seq))
                    })
                current_header = line.strip()[1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_header and current_seq:
            sequences.append({
                'header': current_header,
                'length': len(''.join(current_seq))
            })

    return sequences


def main():
    print("=" * 80)
    print("📊 OUTPUT FNA STATISTICS GENERATOR")
    print("   Creates per-taxon stats.csv and genera.csv files")
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
            prefix = parts[0]
            taxon_name = '_'.join(parts[1:-1])
            taxon_rank = parts[-1]
        else:
            prefix = "000"
            taxon_name = fna_file.stem
            taxon_rank = 'unknown'

        sequences = analyze_fna_file(fna_file)

        if not sequences:
            continue

        lengths = [s['length'] for s in sequences]

        # Get taxa info
        info = taxa_info.get(taxon_name, {})
        genera = info.get('genera', [])
        genera_with_scores = info.get('genera_with_scores', [])

        # --- Write individual stats.csv ---
        stats_file = STATS_DIR / f"{prefix}_{taxon_name}_{taxon_rank}_stats.csv"
        with open(stats_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['taxon', 'rank', 'lineage', 'seq_count', 'avg_length', 'min_length', 'max_length', 'priority_score', 'expected_otu', 'genera_count', 'genera_list'])
            writer.writerow([
                taxon_name,
                taxon_rank,
                info.get('path', 'N/A'),
                len(lengths),
                round(sum(lengths) / len(lengths), 1),
                min(lengths),
                max(lengths),
                info.get('priority', 0),
                info.get('otu', 0),
                len(genera),
                '; '.join(genera)
            ])

        # --- Write individual genera.csv ---
        genera_file = STATS_DIR / f"{prefix}_{taxon_name}_{taxon_rank}_genera.csv"
        with open(genera_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['genus', 'priority_score', 'parent_taxon', 'parent_rank', 'lineage'])

            # Use genera_with_scores if available, otherwise just list genera
            if genera_with_scores:
                for g in sorted(genera_with_scores, key=lambda x: -x['score']):
                    writer.writerow([
                        g['name'],
                        g['score'],
                        taxon_name,
                        taxon_rank,
                        info.get('path', 'N/A')
                    ])
            else:
                for genus in genera:
                    writer.writerow([
                        genus,
                        '',  # No score available
                        taxon_name,
                        taxon_rank,
                        info.get('path', 'N/A')
                    ])

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
            'genera_count': len(genera),
        }
        results.append(result)
        total_seqs += len(lengths)

        print(f"  {fna_file.name}: {len(lengths)} seqs, {len(genera)} genera")

    # Write overall summary CSV
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
    print(f"\n📁 Output written to: {STATS_DIR}")
    print(f"   - output_summary.csv (overall summary)")
    print(f"   - *_stats.csv (per-taxon stats)")
    print(f"   - *_genera.csv (per-taxon genera breakdown)")


if __name__ == "__main__":
    main()

