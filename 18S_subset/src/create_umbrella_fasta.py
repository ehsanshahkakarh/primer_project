#!/usr/bin/env python3
"""
Create combined FASTA files for umbrella primer design.
Groups sequences from multiple taxa under common umbrella categories.

Also includes high-level .U. (unassigned) taxa that don't fall under any
specific order/family target but belong to the umbrella's taxonomic group.
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PROJECT_DIR = Path(__file__).parent.parent.absolute()
PRIORITY_CSV = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "18s" / "nf_abundance_priority_scores.csv"
INPUT_DIR = PROJECT_DIR / "output"  # Individual taxon FNA files
OUTPUT_DIR = PROJECT_DIR / "umbrella_output"

# Source data for extracting .U. sequences directly
EUKCENSUS_DIR = PROJECT_DIR.parent.parent / "eukcensus_metadata"
CLUSTER_TSV = EUKCENSUS_DIR / "eukcensus_18S.clusters.97.tsv"
FASTA_FILE = EUKCENSUS_DIR / "eukcensus_2025_18S.97.fna"

# Mapping of high-level .U. genera to their umbrella groups
# These are taxa classified at phylum/division level that don't roll up to any order/family
# Format: 'genus_name': ['umbrella1', 'umbrella2', ...]  (first match is used)
UNASSIGNED_TO_UMBRELLA = {
    # Amoebozoa-related
    'Amoebozoa.U.genus': ['Amoebozoa', 'Tubulinea', 'Evosea'],

    # Rhizaria-related
    'Cercozoa.U.genus': ['Rhizaria'],
    'Rhizaria.U.genus': ['Rhizaria'],
    'Monothalamids_X.U.genus': ['Rhizaria', 'Retaria'],

    # Stramenopiles-related
    'Stramenopiles.U.genus': ['Stramenopiles', 'Bigyra', 'Sar'],
    'Gyrista.U.genus': ['Stramenopiles', 'Sar'],
    'Ochrophyta.U.genus': ['Stramenopiles', 'Sar'],

    # Alveolata-related
    'Ciliophora.U.genus': ['Intramacronucleata', 'Oligohymenophorea', 'Phyllopharyngea'],
    'Dinophyceae.U.genus': ['Sar'],
    'Dino-Group-II.U.genus': ['Sar'],
    'Alveolata.U.genus': ['Sar'],
    'Apicomplexa.U.genus': ['Apicomplexa', 'Conoidasida'],

    # Discoba/Euglenozoa-related
    'Euglenida.U.genus': ['Euglenozoa', 'Spirocuta'],
    'Heterolobosea.U.genus': ['Euglenozoa'],
    'Kinetoplastea.U.genus': ['Euglenozoa'],
    'Petalomonadida.U.genus': ['Euglenozoa', 'Spirocuta'],
    'Discoba.U.genus': ['Euglenozoa'],

    # SAR supergroup
    'TSAR.U.genus': ['Sar', 'Stramenopiles', 'Rhizaria'],
    'SAR.U.genus': ['Sar'],

    # Chlorophyta
    'Chlorophyta_X.U.genus': ['Sar'],  # No direct umbrella, add to SAR for coverage

    # Fungi-related (under Opisthokonta, but we have sordariomyceta umbrella)
    'Fungi.U.genus': ['sordariomyceta'],
    'Pezizomycotina.U.genus': ['sordariomyceta'],
    'Agaricomycetes.U.genus': ['sordariomyceta'],
    'Dothideomycetes.U.genus': ['sordariomyceta'],
    'Sordariomycetes.U.genus': ['sordariomyceta'],
    'Chytridiomycota.U.genus': ['sordariomyceta'],
}


def parse_path_to_levels(path_str):
    """Parse taxonomic path into levels."""
    return [p.strip() for p in path_str.split(' > ')]


def load_cluster_lookup():
    """
    Load the cluster TSV and create a lookup from centroid ID -> genus.
    Returns dict: {centroid_id: genus_name}
    """
    print(f"📖 Loading cluster taxonomy from {CLUSTER_TSV.name}...")
    lookup = {}

    with open(CLUSTER_TSV, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            centroid = row['centroid']
            genus = row['genus']
            lookup[centroid] = genus

    print(f"   Loaded {len(lookup):,} clusters")
    return lookup


def get_unassigned_sequences_for_umbrella(umbrella_name, cluster_lookup):
    """
    Find all sequence IDs for high-level .U. genera that should be in this umbrella.
    Returns set of centroid IDs.
    """
    matching_ids = set()
    matched_genera = []

    for genus, umbrella_list in UNASSIGNED_TO_UMBRELLA.items():
        if umbrella_name in umbrella_list:
            # Find all centroids with this genus
            for centroid, centroid_genus in cluster_lookup.items():
                if centroid_genus == genus:
                    matching_ids.add(centroid)
            if genus not in matched_genera:
                matched_genera.append(genus)

    return matching_ids, matched_genera


def extract_sequences_from_fasta(fasta_file, sequence_ids, output_handle, umbrella_name, source_tag="unassigned"):
    """
    Extract sequences from a FASTA file and write to output handle.
    Returns count of sequences written.
    """
    count = 0
    writing = False

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]  # Get first part of header
                if seq_id in sequence_ids:
                    # Write modified header
                    header = line.strip()
                    output_handle.write(f"{header}|umbrella={umbrella_name}|source={source_tag}\n")
                    count += 1
                    writing = True
                else:
                    writing = False
            elif writing:
                output_handle.write(line)

    return count


def get_umbrella_groups():
    """Load priority data and group by umbrella category."""
    
    umbrella_groups = defaultdict(list)
    
    with open(PRIORITY_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            path_parts = parse_path_to_levels(row['path'])
            root = path_parts[0] if path_parts else 'Unknown'
            
            umbrella_groups[root].append({
                'node': row['node'],
                'rank': row['rank'],
                'total_otu': int(row['total_otu']),
                'priority_score': float(row['total_priority_score']),
            })
    
    return umbrella_groups


def find_fna_file(node_name, rank):
    """Find the FNA file for a given taxon."""
    # Files are named like: 001_Euamoebida_order.fna
    for fna_file in INPUT_DIR.glob("*.fna"):
        # Check if node name and rank are in filename
        if node_name in fna_file.stem and rank in fna_file.stem:
            return fna_file
    
    # Try partial match
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
                        # Add umbrella group to header
                        header = line.strip()
                        out.write(f"{header}|umbrella={umbrella_name}\n")
                        total_seqs += 1
                    else:
                        out.write(line)

    return total_seqs


def main():
    print("=" * 80)
    print("🎯 CREATING UMBRELLA FASTA FILES FOR PRIMER DESIGN")
    print("   (Now includes high-level .U. unassigned taxa)")
    print("=" * 80)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load cluster lookup for .U. sequence extraction
    cluster_lookup = load_cluster_lookup()

    # Get umbrella groups
    umbrella_groups = get_umbrella_groups()

    # Calculate priority scores
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

    # Sort by priority
    sorted_umbrellas = sorted(umbrella_priority.items(),
                              key=lambda x: -x[1]['priority'])

    print(f"\n📋 Found {len(sorted_umbrellas)} umbrella groups")
    print(f"📁 Input directory: {INPUT_DIR}")
    print(f"📁 Output directory: {OUTPUT_DIR}")
    print(f"📁 Master FASTA: {FASTA_FILE.name}\n")

    # Process top umbrella groups
    results = []

    for i, (name, data) in enumerate(sorted_umbrellas[:20], 1):  # Top 20
        print(f"\n🌳 [{i:2d}] {name}")
        print(f"      Priority: {data['priority']:.1f}, Taxa: {data['taxa_count']}, OTUs: {data['otu']}")

        # Find FNA files for each member taxon
        fna_files = []
        found_taxa = []
        missing_taxa = []

        for member in data['members']:
            fna_file = find_fna_file(member['node'], member['rank'])
            if fna_file:
                fna_files.append(fna_file)
                found_taxa.append(member['node'])
            else:
                missing_taxa.append(member['node'])

        print(f"      Found FNA: {len(fna_files)}/{len(data['members'])} taxa")

        # Create combined FASTA
        safe_name = name.replace(' ', '_').replace('/', '_')
        output_file = OUTPUT_DIR / f"{i:02d}_{safe_name}_umbrella.fna"

        # First combine existing FNA files
        total_seqs = combine_fasta_files(fna_files, output_file, name)

        # Now add high-level .U. sequences for this umbrella
        unassigned_ids, matched_genera = get_unassigned_sequences_for_umbrella(name, cluster_lookup)
        unassigned_count = 0

        if unassigned_ids:
            print(f"      🔍 Found {len(unassigned_ids)} .U. sequences to add from: {', '.join(matched_genera)}")
            with open(output_file, 'a') as out:
                unassigned_count = extract_sequences_from_fasta(
                    FASTA_FILE, unassigned_ids, out, name, "high_level_unassigned"
                )
            print(f"      ➕ Added {unassigned_count} high-level .U. sequences")

        total_seqs += unassigned_count
        print(f"      ✅ Created: {output_file.name} ({total_seqs:,} total sequences)")

        results.append({
            'rank': i,
            'umbrella': name,
            'priority': data['priority'],
            'taxa_in_group': data['taxa_count'],
            'taxa_with_fna': len(fna_files),
            'sequences_from_fna': total_seqs - unassigned_count,
            'sequences_from_unassigned': unassigned_count,
            'total_sequences': total_seqs,
            'output_file': output_file.name,
            'included_taxa': ';'.join(found_taxa),
            'unassigned_genera': ';'.join(matched_genera),
            'missing_taxa': ';'.join(missing_taxa)
        })

    # Write summary
    summary_file = OUTPUT_DIR / "umbrella_summary.csv"
    with open(summary_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\n\n{'='*80}")
    print("✅ UMBRELLA FASTA FILES CREATED")
    print(f"{'='*80}")
    print(f"   Output directory: {OUTPUT_DIR}")
    print(f"   Summary: {summary_file}")
    print(f"   Files created: {len(results)}")

    # Summary of .U. additions
    total_unassigned = sum(r['sequences_from_unassigned'] for r in results)
    print(f"\n   🔬 High-level .U. sequences added: {total_unassigned:,}")

    print("\nNext step: Run ppdesign on these umbrella files with:")
    print("   ppdesign primer main --fasta-input <umbrella.fna> ...")


if __name__ == "__main__":
    main()

