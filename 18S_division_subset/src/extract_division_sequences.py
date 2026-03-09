#!/usr/bin/env python3
"""
Extract 18S sequences for DIVISION-level primer targets from nf_division_priority_scores.csv

For each division, extracts ALL sequences that fall under that taxonomic group
from the EukCensus 18S FNA file.
"""

import csv
import ast
import sys
import re
from pathlib import Path

# Increase CSV field size limit for large member lists
csv.field_size_limit(sys.maxsize)

# Directory paths
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_DIR = SCRIPT_DIR.parent
BRANCH_GAP_DIR = PROJECT_DIR.parent / "branch_gap_analysis"
EUKCENSUS_DIR = PROJECT_DIR.parent.parent / "eukcensus_metadata"

# Input files
PRIORITY_SCORES_CSV = BRANCH_GAP_DIR / "output" / "18s_division" / "nf_division_priority_scores.csv"
CLUSTER_TSV = EUKCENSUS_DIR / "eukcensus_18S.clusters.97.tsv"
FASTA_FILE = EUKCENSUS_DIR / "eukcensus_2025_18S.97.fna"

# Output directory
OUTPUT_DIR = PROJECT_DIR / "output"


def load_division_targets(csv_file: Path) -> list[dict]:
    """Load all division targets from priority scores CSV"""
    print(f"📖 Loading division targets from {csv_file.name}...")
    
    targets = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_genera_str = row.get('all_genera', '')
            all_genera = [g.strip() for g in all_genera_str.split(';') if g.strip()]
            
            targets.append({
                'division': row['division'],
                'total_otu': int(row['total_otu']),
                'all_genera': all_genera
            })
    
    print(f"   ✓ Loaded {len(targets)} division targets")
    return targets


def build_taxonomy_lookup(tsv_file: Path) -> dict:
    """Build lookup: {sequence_id: genus_name}"""
    print(f"\n📖 Building taxonomy lookup from {tsv_file.name}...")
    
    lookup = {}
    cluster_count = 0
    
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            cluster_count += 1
            centroid = row['centroid']
            genus = row['genus']
            members_str = row['members']
            
            lookup[centroid] = genus
            
            try:
                members = ast.literal_eval(members_str)
                for member in members:
                    if member != centroid:
                        lookup[member] = genus
            except (ValueError, SyntaxError):
                pass
    
    print(f"   ✓ Processed {cluster_count} clusters")
    print(f"   ✓ Total sequences indexed: {len(lookup)}")
    return lookup


def find_matching_sequences(target: dict, lookup: dict) -> set:
    """Find all sequence IDs matching genera in the target division"""
    all_genera = set(target.get('all_genera', []))
    
    matching_ids = set()
    for seq_id, genus in lookup.items():
        if genus in all_genera:
            matching_ids.add(seq_id)
    
    return matching_ids


def sanitize_filename(name: str) -> str:
    """Convert a taxonomic name to a safe filename"""
    return re.sub(r'[^\w\-]', '_', name)


def extract_sequences_to_fasta(fasta_file: Path, sequence_ids: set, output_file: Path) -> int:
    """Extract sequences matching the IDs from FASTA file."""
    extracted_count = 0
    current_header = None
    current_seq = []
    write_current = False

    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            line = line.strip()

            if line.startswith('>'):
                if write_current and current_header:
                    fout.write(f"{current_header}\n")
                    fout.write('\n'.join(current_seq) + '\n')
                    extracted_count += 1

                current_header = line
                current_seq = []
                seq_id = line[1:].split()[0]
                write_current = seq_id in sequence_ids
            else:
                current_seq.append(line)

        if write_current and current_header:
            fout.write(f"{current_header}\n")
            fout.write('\n'.join(current_seq) + '\n')
            extracted_count += 1

    return extracted_count


def main():
    print("=" * 70)
    print("🧬 EXTRACTING 18S SEQUENCES FOR DIVISION-LEVEL TARGETS")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load targets and build lookup
    targets = load_division_targets(PRIORITY_SCORES_CSV)
    lookup = build_taxonomy_lookup(CLUSTER_TSV)

    print(f"\n🔬 Processing {len(targets)} division targets...")
    print("-" * 70)

    results = []
    num_digits = len(str(len(targets)))

    for i, target in enumerate(targets, 1):
        division = target['division']
        expected_otu = target['total_otu']

        matching_ids = find_matching_sequences(target, lookup)

        if len(matching_ids) == 0:
            print(f"   ⚠️  [{i:3d}/{len(targets)}] {division}: No sequences found")
            continue

        prefix = str(i).zfill(num_digits)
        safe_name = sanitize_filename(division)
        output_file = OUTPUT_DIR / f"{prefix}_{safe_name}_division.fna"

        extracted = extract_sequences_to_fasta(FASTA_FILE, matching_ids, output_file)

        print(f"   ✓  [{i:3d}/{len(targets)}] {division}: {extracted} sequences → {output_file.name}")

        results.append({
            'division': division,
            'expected_otu': expected_otu,
            'sequences_found': len(matching_ids),
            'sequences_extracted': extracted,
            'output_file': output_file.name
        })

    # Summary
    print("\n" + "=" * 70)
    print("✅ EXTRACTION COMPLETE")
    print("=" * 70)
    print(f"\n📁 Output directory: {OUTPUT_DIR}")
    print(f"   Total divisions processed: {len(results)}")
    print(f"   Total FNA files created: {len([r for r in results if r['sequences_extracted'] > 0])}")

    # Write summary
    summary_file = OUTPUT_DIR / "extraction_summary.csv"
    with open(summary_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['division', 'expected_otu', 'sequences_found', 'sequences_extracted', 'output_file'])
        writer.writeheader()
        writer.writerows(results)
    print(f"   Summary saved to: {summary_file.name}")


if __name__ == "__main__":
    main()

