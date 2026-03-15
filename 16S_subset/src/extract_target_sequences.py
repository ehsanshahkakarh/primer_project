#!/usr/bin/env python3
"""
Extract 16S rRNA sequences for priority taxa from EukCensus master FNA.

Key differences from 18S version:
- Uses UC format cluster file (USEARCH/VSEARCH output)
- Loads priority scores from both Archaea and Bacteria domains
- Same extraction logic: centroid + all cluster members
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

# Paths
PROJECT_DIR = Path(__file__).parent.parent.absolute()
EUKCENSUS_DIR = PROJECT_DIR.parent.parent / "eukcensus_metadata"
UC_FILE = EUKCENSUS_DIR / "eukcensus_2025_16S.97.uc"
FASTA_FILE = EUKCENSUS_DIR / "eukcensus_2025_16S.97.fna"
OUTPUT_DIR = PROJECT_DIR / "output"

# Priority score files for both domains
PRIORITY_DIR = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "16s"
ARCHAEA_PRIORITY = PRIORITY_DIR / "archaea" / "nf_abundance_priority_scores.csv"
BACTERIA_PRIORITY = PRIORITY_DIR / "bacteria" / "nf_abundance_priority_scores.csv"


def load_priority_scores():
    """Load priority scores from both Archaea and Bacteria CSV files."""
    targets = []
    
    for domain, priority_file in [("archaea", ARCHAEA_PRIORITY), ("bacteria", BACTERIA_PRIORITY)]:
        if not priority_file.exists():
            print(f"⚠️  {domain} priority file not found: {priority_file}")
            continue
            
        with open(priority_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                targets.append({
                    'domain': domain,
                    'node': row['node'],
                    'rank': row['rank'],
                    'priority_score': float(row['total_priority_score']),
                    'total_otu': int(row['total_otu']),
                    'path': row['path'],
                    'all_genera': row.get('all_genera', ''),
                })
    
    # Sort by priority score descending
    targets.sort(key=lambda x: -x['priority_score'])
    print(f"📋 Loaded {len(targets)} priority targets (Archaea + Bacteria)")
    return targets


def parse_uc_file():
    """
    Parse UC format file to build cluster membership.
    UC columns: Type, ClusterNum, Length, %Id, Strand, Query, Target, Compressed, QueryLabel, TargetLabel
    
    Returns:
        - genus_to_centroids: {genus_name: set(centroid_ids)}
        - centroid_to_members: {centroid_id: set(member_ids)}
    """
    print(f"📖 Parsing UC file: {UC_FILE.name}...")
    
    genus_to_centroids = defaultdict(set)
    centroid_to_members = defaultdict(set)
    cluster_centroids = {}  # cluster_num -> centroid_id
    
    with open(UC_FILE, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            rec_type = parts[0]  # S=Seed/Centroid, H=Hit/Member
            cluster_num = parts[1]
            query_label = parts[8]  # Sequence ID
            
            if rec_type == 'S':
                # This is a centroid/seed
                cluster_centroids[cluster_num] = query_label
                centroid_to_members[query_label].add(query_label)  # Centroid is also a member
            elif rec_type == 'H':
                # This is a hit/member assigned to a cluster
                target_label = parts[9] if len(parts) > 9 else cluster_centroids.get(cluster_num)
                if target_label:
                    centroid_to_members[target_label].add(query_label)
    
    print(f"   Found {len(cluster_centroids):,} clusters")
    print(f"   Total member assignments: {sum(len(m) for m in centroid_to_members.values()):,}")
    
    return centroid_to_members, cluster_centroids


def extract_genus_from_header(header):
    """Extract genus name from sequence header if possible."""
    if '_' in header:
        parts = header.split('_')
        for part in parts:
            if part and part[0].isupper() and len(part) > 3:
                if part not in ['IMG', 'REF', 'SPR', 'Bacteria', 'Archaea', 'Eukaryota']:
                    return part
    return None


def load_fasta_index(fasta_file):
    """Build index of sequence IDs to file positions for fast extraction."""
    print(f"📖 Indexing FASTA file: {fasta_file.name}...")
    index = {}
    with open(fasta_file, 'r') as f:
        pos = 0
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                index[seq_id] = pos
            pos += len(line)
    print(f"   Indexed {len(index):,} sequences")
    return index


def extract_sequences(fasta_file, seq_ids, output_file):
    """Extract sequences by ID from FASTA file."""
    seq_ids_set = set(seq_ids)
    found = 0
    writing = False

    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                if seq_id in seq_ids_set:
                    fout.write(line)
                    writing = True
                    found += 1
                else:
                    writing = False
            elif writing:
                fout.write(line)

    return found


def find_sequences_for_taxon(taxon_name, all_genera, centroid_to_members, fasta_index):
    """Find all sequence IDs for a given taxon based on genera list."""
    genera_list = [g.strip() for g in all_genera.split(';') if g.strip()]
    matching_ids = set()

    # Search through FASTA index for matching genera in headers
    for seq_id in fasta_index:
        for genus in genera_list:
            if genus in seq_id:
                matching_ids.add(seq_id)
                # Also add cluster members if this is a centroid
                if seq_id in centroid_to_members:
                    matching_ids.update(centroid_to_members[seq_id])
                break

    return matching_ids


def main():
    print("=" * 80)
    print("🧬 16S rRNA SEQUENCE EXTRACTION FOR PRIMER DESIGN")
    print("   (Archaea + Bacteria combined)")
    print("=" * 80)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load priority targets
    targets = load_priority_scores()

    # Parse UC file for cluster membership
    centroid_to_members, cluster_centroids = parse_uc_file()

    # Index FASTA file
    fasta_index = load_fasta_index(FASTA_FILE)

    print(f"\n🎯 Extracting sequences for top priority targets...\n")

    results = []
    for i, target in enumerate(targets[:300], 1):  # Top 300 targets
        taxon_name = target['node']
        rank = target['rank']
        domain = target['domain']

        # Find matching sequences
        matching_ids = find_sequences_for_taxon(
            taxon_name, target['all_genera'], centroid_to_members, fasta_index
        )

        if not matching_ids:
            continue

        # Create output file
        safe_name = taxon_name.replace(' ', '_').replace('/', '_')
        output_file = OUTPUT_DIR / f"{i:03d}_{safe_name}_{rank}.fna"

        # Extract sequences
        found = extract_sequences(FASTA_FILE, matching_ids, output_file)

        if found > 0:
            print(f"  [{i:3d}] {domain:8s} | {taxon_name:30s} | {rank:8s} | {found:5d} seqs")
            results.append({
                'rank_num': i,
                'domain': domain,
                'taxon': taxon_name,
                'taxon_rank': rank,
                'priority': target['priority_score'],
                'sequences': found,
                'file': output_file.name
            })

    # Write summary
    summary_file = OUTPUT_DIR / "extraction_summary.csv"
    with open(summary_file, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)

    print(f"\n{'='*80}")
    print(f"✅ Extraction complete!")
    print(f"   Taxa extracted: {len(results)}")
    print(f"   Total sequences: {sum(r['sequences'] for r in results):,}")
    print(f"   Output directory: {OUTPUT_DIR}")
    print(f"   Summary: {summary_file}")
    print("=" * 80)


if __name__ == "__main__":
    main()

