#!/usr/bin/env python3
"""
Extract 18S sequences for all primer targets from nf_abundance_priority_scores.csv

For each target (family/order/genus/etc.), extracts ALL sequences that fall under
that taxonomic node from the EukCensus 18S FNA file.

Workflow:
1. Read priority scores CSV to get list of targets (node names + ranks)
2. Read cluster TSV to build a lookup: sequence_id -> taxonomy info
3. For each target, find all matching sequences based on rank
4. Extract matching sequences from the FNA file
5. Write separate FNA files for each target
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
PRIORITY_SCORES_CSV = BRANCH_GAP_DIR / "output" / "18s" / "nf_abundance_priority_scores.csv"
CLUSTER_TSV = EUKCENSUS_DIR / "eukcensus_18S.clusters.97.tsv"
FASTA_FILE = EUKCENSUS_DIR / "eukcensus_2025_18S.97.fna"

# Output directory
OUTPUT_DIR = PROJECT_DIR / "output"


def load_priority_targets(csv_file: Path) -> list[dict]:
    """Load all primer targets from priority scores CSV"""
    print(f"📖 Loading primer targets from {csv_file.name}...")

    targets = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse all_genera field (semicolon-separated list of genus names)
            all_genera_str = row.get('all_genera', '')
            all_genera = [g.strip() for g in all_genera_str.split(';') if g.strip()]

            targets.append({
                'node': row['node'],
                'rank': row['rank'],
                'total_otu': int(row['total_otu']),
                'path': row['path'],
                'all_genera': all_genera  # List of all genus names under this node
            })

    print(f"   ✓ Loaded {len(targets)} primer targets")
    return targets


def build_taxonomy_lookup(tsv_file: Path) -> dict:
    """
    Build a lookup from the cluster TSV for matching targets.
    Returns: {sequence_id: {division, family, genus, all taxonomy fields}}
    Also extracts higher taxonomy from the path where available.
    """
    print(f"\n📖 Building taxonomy lookup from {tsv_file.name}...")
    
    lookup = {}
    cluster_count = 0
    
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            cluster_count += 1
            centroid = row['centroid']
            division = row['division']
            family = row['family']
            genus = row['genus']
            members_str = row['members']
            
            # Create taxonomy info dict
            tax_info = {
                'division': division,
                'family': family, 
                'genus': genus,
            }
            
            # Add centroid
            lookup[centroid] = tax_info
            
            # Parse and add all members
            try:
                members = ast.literal_eval(members_str)
                for member in members:
                    if member != centroid:
                        lookup[member] = tax_info
            except (ValueError, SyntaxError):
                pass
    
    print(f"   ✓ Processed {cluster_count} clusters")
    print(f"   ✓ Total sequences indexed: {len(lookup)}")
    return lookup


def find_matching_sequences(target: dict, lookup: dict) -> set:
    """
    Find all sequence IDs that match a given taxonomic target.
    Uses the all_genera list from the priority scores CSV to match against genus names.
    """
    all_genera = set(target.get('all_genera', []))

    if not all_genera:
        # Fallback: try to match by node name directly
        node_name = target['node']
        matching_ids = set()
        for seq_id, tax_info in lookup.items():
            if (tax_info['family'] == node_name or
                tax_info['genus'] == node_name or
                tax_info['division'] == node_name):
                matching_ids.add(seq_id)
        return matching_ids

    # Match by genus names from the all_genera list
    matching_ids = set()
    for seq_id, tax_info in lookup.items():
        if tax_info['genus'] in all_genera:
            matching_ids.add(seq_id)

    return matching_ids


def sanitize_filename(name: str) -> str:
    """Convert a taxonomic name to a safe filename"""
    # Replace problematic characters
    safe_name = re.sub(r'[^\w\-]', '_', name)
    return safe_name


def extract_sequences_to_fasta(fasta_file: Path, sequence_ids: set, output_file: Path) -> int:
    """
    Extract sequences matching the IDs from FASTA file.
    Returns number of sequences extracted.
    """
    extracted_count = 0
    current_header = None
    current_seq = []
    write_current = False

    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            line = line.strip()

            if line.startswith('>'):
                # Write previous sequence if it matched
                if write_current and current_header:
                    fout.write(f"{current_header}\n")
                    fout.write('\n'.join(current_seq) + '\n')
                    extracted_count += 1

                current_header = line
                current_seq = []

                # Extract sequence ID from header (everything after > until first space)
                seq_id = line[1:].split()[0]
                write_current = seq_id in sequence_ids
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if write_current and current_header:
            fout.write(f"{current_header}\n")
            fout.write('\n'.join(current_seq) + '\n')
            extracted_count += 1

    return extracted_count


def main():
    print("=" * 70)
    print("🧬 EXTRACTING 18S SEQUENCES FOR ALL PRIMER TARGETS")
    print("=" * 70)

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load primer targets
    targets = load_priority_targets(PRIORITY_SCORES_CSV)

    # Step 2: Build taxonomy lookup
    lookup = build_taxonomy_lookup(CLUSTER_TSV)

    # Step 3: Process each target
    print(f"\n🔬 Processing {len(targets)} primer targets...")
    print("-" * 70)

    results = []
    for i, target in enumerate(targets, 1):
        node_name = target['node']
        rank = target['rank']
        expected_otu = target['total_otu']

        # Find matching sequences
        matching_ids = find_matching_sequences(target, lookup)

        if len(matching_ids) == 0:
            print(f"   ⚠️  [{i:3d}/{len(targets)}] {node_name} ({rank}): No sequences found")
            continue

        # Create output filename with numeric prefix to preserve order
        # Use zero-padded index based on total number of targets
        num_digits = len(str(len(targets)))
        prefix = str(i).zfill(num_digits)
        safe_name = sanitize_filename(node_name)
        output_file = OUTPUT_DIR / f"{prefix}_{safe_name}_{rank}.fna"

        # Extract sequences
        extracted = extract_sequences_to_fasta(FASTA_FILE, matching_ids, output_file)

        print(f"   ✓  [{i:3d}/{len(targets)}] {node_name} ({rank}): {extracted} sequences → {output_file.name}")

        results.append({
            'node': node_name,
            'rank': rank,
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
    print(f"   Total targets processed: {len(results)}")
    print(f"   Total FNA files created: {len([r for r in results if r['sequences_extracted'] > 0])}")

    # Write summary CSV
    summary_file = OUTPUT_DIR / "extraction_summary.csv"
    with open(summary_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['node', 'rank', 'expected_otu', 'sequences_found', 'sequences_extracted', 'output_file'])
        writer.writeheader()
        writer.writerows(results)
    print(f"   Summary saved to: {summary_file.name}")


if __name__ == "__main__":
    main()

