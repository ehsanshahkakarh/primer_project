#!/usr/bin/env python3
"""
Create combined FASTA files for umbrella primer design (16S version).
Groups sequences from multiple taxa under common umbrella categories.
Handles both Archaea and Bacteria domains.
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PROJECT_DIR = Path(__file__).parent.parent.absolute()
PRIORITY_DIR = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "16s"
INPUT_DIR = PROJECT_DIR / "output"
OUTPUT_DIR = PROJECT_DIR / "umbrella_output"


def parse_path_to_levels(path_str):
    """Parse taxonomic path into levels."""
    return [p.strip() for p in path_str.split(' > ')]


def load_all_priority_data():
    """Load priority data from both Archaea and Bacteria files."""
    umbrella_groups = defaultdict(list)
    
    for domain in ['archaea', 'bacteria']:
        priority_file = PRIORITY_DIR / domain / "nf_abundance_priority_scores.csv"
        if not priority_file.exists():
            print(f"⚠️  {domain} priority file not found")
            continue
            
        with open(priority_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                path_parts = parse_path_to_levels(row['path'])
                # Use first taxonomic level as umbrella (e.g., phylum or class)
                root = path_parts[0] if path_parts else domain.capitalize()
                
                umbrella_groups[root].append({
                    'domain': domain,
                    'node': row['node'],
                    'rank': row['rank'],
                    'total_otu': int(row['total_otu']),
                    'priority_score': float(row['total_priority_score']),
                })
    
    return umbrella_groups


def find_fna_file(node_name, rank):
    """Find the FNA file for a given taxon."""
    for fna_file in INPUT_DIR.glob("*.fna"):
        if node_name in fna_file.stem and rank in fna_file.stem:
            return fna_file
    for fna_file in INPUT_DIR.glob("*.fna"):
        if node_name in fna_file.stem:
            return fna_file
    return None


def combine_fasta_files(fna_files, output_file, umbrella_name):
    """Combine multiple FASTA files into one."""
    total_seqs = 0
    with open(output_file, 'w') as out:
        for fna_file in fna_files:
            if not fna_file.exists():
                continue
            with open(fna_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        header = line.strip()
                        out.write(f"{header}|umbrella={umbrella_name}\n")
                        total_seqs += 1
                    else:
                        out.write(line)
    return total_seqs


def main():
    print("=" * 80)
    print("🎯 CREATING 16S UMBRELLA FASTA FILES FOR PRIMER DESIGN")
    print("   (Archaea + Bacteria combined)")
    print("=" * 80)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    umbrella_groups = load_all_priority_data()

    # Calculate priority scores per umbrella
    umbrella_priority = {}
    for name, members in umbrella_groups.items():
        total_priority = sum(m['priority_score'] for m in members)
        total_otu = sum(m['total_otu'] for m in members)
        umbrella_priority[name] = {
            'priority': total_priority,
            'otu': total_otu,
            'taxa_count': len(members),
            'members': members
        }

    sorted_umbrellas = sorted(umbrella_priority.items(), key=lambda x: -x[1]['priority'])
    print(f"\n📋 Found {len(sorted_umbrellas)} umbrella groups")
    print(f"📁 Input: {INPUT_DIR}")
    print(f"📁 Output: {OUTPUT_DIR}\n")

    results = []
    for i, (name, data) in enumerate(sorted_umbrellas[:30], 1):
        print(f"\n🌳 [{i:2d}] {name}")
        print(f"      Priority: {data['priority']:.1f}, Taxa: {data['taxa_count']}, OTUs: {data['otu']}")

        fna_files = []
        found_taxa = []
        for member in data['members']:
            fna_file = find_fna_file(member['node'], member['rank'])
            if fna_file:
                fna_files.append(fna_file)
                found_taxa.append(member['node'])

        print(f"      Found FNA: {len(fna_files)}/{len(data['members'])} taxa")

        safe_name = name.replace(' ', '_').replace('/', '_')
        output_file = OUTPUT_DIR / f"{i:02d}_{safe_name}_umbrella.fna"
        total_seqs = combine_fasta_files(fna_files, output_file, name)
        print(f"      ✅ Created: {output_file.name} ({total_seqs:,} sequences)")

        results.append({
            'rank': i, 'umbrella': name, 'priority': data['priority'],
            'taxa_in_group': data['taxa_count'], 'taxa_with_fna': len(fna_files),
            'total_sequences': total_seqs, 'output_file': output_file.name,
        })

    summary_file = OUTPUT_DIR / "umbrella_summary.csv"
    with open(summary_file, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)

    print(f"\n{'='*80}")
    print("✅ UMBRELLA FILES CREATED")
    print(f"   Output: {OUTPUT_DIR}")
    print(f"   Files: {len(results)}")
    print("=" * 80)


if __name__ == "__main__":
    main()

