#!/usr/bin/env python3
"""
Generate detailed stats for each umbrella FASTA file.
Shows sequences grouped by source taxon with counts and avg lengths.
Also includes genus-level breakdown.
"""

import csv
import re
import sys
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PROJECT_DIR = Path(__file__).parent.parent.absolute()
PRIORITY_CSV = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "18s" / "nf_abundance_priority_scores.csv"
INPUT_DIR = PROJECT_DIR / "output"  # Individual taxon FNA files
UMBRELLA_DIR = PROJECT_DIR / "umbrella_output"
STATS_DIR = PROJECT_DIR / "umbrella_stats"


def parse_path_to_levels(path_str):
    """Parse taxonomic path into levels."""
    return [p.strip() for p in path_str.split(' > ')]


def get_taxa_info():
    """Load priority data with full lineage info."""
    taxa_info = {}
    with open(PRIORITY_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            path_parts = parse_path_to_levels(row['path'])
            root = path_parts[0] if path_parts else 'Unknown'
            # Parse top_genera to get genus with scores
            top_genera_raw = row.get('top_genera', '')
            genera_with_scores = []
            for g in re.findall(r'(\w+)\(([0-9.]+)\)', top_genera_raw):
                genera_with_scores.append({'name': g[0], 'score': float(g[1])})

            taxa_info[row['node']] = {
                'rank': row['rank'],
                'path': row['path'],
                'root': root,
                'priority': float(row['total_priority_score']),
                'otu': int(row['total_otu']),
                'genera': [g.strip() for g in row.get('all_genera', '').split(';') if g.strip()],
                'genera_with_scores': genera_with_scores
            }
    return taxa_info


def find_source_taxon(header, taxa_info):
    """Try to identify source taxon from sequence header."""
    # Check for REF sequences with full lineage
    if 'REF_SPR_' in header or 'REF_' in header:
        parts = header.split('_')
        # Look for known taxa in header
        for taxon in taxa_info:
            if taxon in header:
                return taxon, taxa_info[taxon]['path']
    return None, None


def analyze_umbrella_file(fna_file, taxa_info):
    """Analyze sequences in an umbrella FASTA file."""
    # Get umbrella name from filename
    umbrella_name = fna_file.stem.split('_', 1)[1].replace('_umbrella', '')
    
    # Track sequences by source file they came from
    source_files = defaultdict(list)
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
        
        # Don't forget last sequence
        if current_header and current_seq:
            seq_len = len(''.join(current_seq))
            all_seqs.append({'header': current_header, 'length': seq_len})
    
    return umbrella_name, all_seqs


def match_seqs_to_source_files(umbrella_fna, input_dir):
    """Match umbrella sequences back to their source FNA files."""
    # Build index of all source sequences
    source_index = {}  # header -> (taxon_name, taxon_rank, file_path)
    
    for fna_file in input_dir.glob("*.fna"):
        # Parse filename like "001_Euamoebida_order.fna"
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
    print("📊 UMBRELLA FASTA STATISTICS GENERATOR")
    print("=" * 80)
    
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Build source index once
    print("\n📁 Building sequence source index...")
    source_index = match_seqs_to_source_files(None, INPUT_DIR)
    print(f"   Indexed {len(source_index):,} sequences from source files")
    
    # Load taxa info for lineage
    taxa_info = get_taxa_info()
    
    # Process each umbrella file
    umbrella_files = sorted(UMBRELLA_DIR.glob("*_umbrella.fna"))
    print(f"\n🔍 Processing {len(umbrella_files)} umbrella files...\n")
    
    for fna_file in umbrella_files:
        umbrella_name, all_seqs = analyze_umbrella_file(fna_file, taxa_info)
        
        # Group by source taxon
        by_taxon = defaultdict(list)
        unmatched = []
        
        # Track high-level unassigned sequences separately
        by_unassigned = defaultdict(list)

        for seq in all_seqs:
            # Check if this is a high-level unassigned sequence
            if '|source=high_level_unassigned' in seq['header']:
                # These are .U. sequences added directly from master FNA
                by_unassigned[('high_level_unassigned', 'phylum/division')].append(seq['length'])
                continue

            # Remove |umbrella=X suffix for matching
            clean_header = seq['header'].split('|umbrella=')[0]

            if clean_header in source_index:
                taxon_name, taxon_rank, source_file = source_index[clean_header]
                by_taxon[(taxon_name, taxon_rank)].append(seq['length'])
            else:
                unmatched.append(seq)

        # Add unassigned to by_taxon for reporting
        for key, lengths in by_unassigned.items():
            by_taxon[key] = lengths
        
        # Write stats file with taxon breakdown
        stats_file = STATS_DIR / f"{fna_file.stem}_stats.csv"
        sorted_taxa = sorted(by_taxon.items(), key=lambda x: -len(x[1]))

        with open(stats_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['taxon', 'rank', 'lineage', 'seq_count', 'avg_length', 'min_length', 'max_length', 'genera'])

            for (taxon_name, taxon_rank), lengths in sorted_taxa:
                # Get lineage and genera from taxa_info
                info = taxa_info.get(taxon_name, {})
                lineage = info.get('path', 'N/A')
                genera = '; '.join(info.get('genera', []))

                writer.writerow([
                    taxon_name,
                    taxon_rank,
                    lineage,
                    len(lengths),
                    round(sum(lengths) / len(lengths), 1),
                    min(lengths),
                    max(lengths),
                    genera
                ])

        # Write detailed genus breakdown file
        genus_file = STATS_DIR / f"{fna_file.stem}_genera.csv"
        all_genera = []

        for (taxon_name, taxon_rank), lengths in sorted_taxa:
            info = taxa_info.get(taxon_name, {})
            lineage = info.get('path', 'N/A')

            for g in info.get('genera_with_scores', []):
                all_genera.append({
                    'genus': g['name'],
                    'priority_score': g['score'],
                    'parent_taxon': taxon_name,
                    'parent_rank': taxon_rank,
                    'lineage': lineage
                })

        # Sort by priority score
        all_genera.sort(key=lambda x: -x['priority_score'])

        with open(genus_file, 'w', newline='') as f:
            if all_genera:
                writer = csv.DictWriter(f, fieldnames=['genus', 'priority_score', 'parent_taxon', 'parent_rank', 'lineage'])
                writer.writeheader()
                writer.writerows(all_genera)

        # Print summary
        total_seqs = len(all_seqs)
        avg_len = sum(s['length'] for s in all_seqs) / total_seqs if total_seqs else 0
        total_genera = len(all_genera)

        print(f"📄 {fna_file.name}")
        print(f"   Total: {total_seqs:,} seqs | Avg length: {avg_len:.0f} bp | Taxa: {len(by_taxon)} | Genera: {total_genera}")
        print(f"   Stats: {stats_file.name} | Genera: {genus_file.name}")
        print()

    print("=" * 80)
    print(f"✅ Stats written to: {STATS_DIR}")
    print("   - *_stats.csv: Taxon breakdown with seq counts and genera list")
    print("   - *_genera.csv: Individual genus breakdown with priority scores")
    print("=" * 80)


if __name__ == "__main__":
    main()

