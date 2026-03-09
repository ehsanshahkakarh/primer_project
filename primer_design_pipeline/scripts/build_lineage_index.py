#!/usr/bin/env python3
"""
Build Lineage Index for Taxonomic Queries

This script parses the census CSV files and builds a reverse mapping from
ancestor taxa to their descendant divisions/families/genera.

This enables queries like "find all families under class Discoba" by looking
up which families have "Discoba" in their lineage.

Usage:
    python build_lineage_index.py --gene 18S --output lineage_index_18S.json
    python build_lineage_index.py --gene 16S --output lineage_index_16S.json
"""

import csv
import json
import argparse
from pathlib import Path
from collections import defaultdict

# Census CSV file paths (relative to project root)
CENSUS_CSV_PATHS = {
    '18S': {
        'division': '_project_pulls/otu_assembly_comparative_pipeline-/18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv',
        'family': '_project_pulls/otu_assembly_comparative_pipeline-/18S_censusparse/csv_outputs/eukcensus_18S_by_family.csv',
        'genus': '_project_pulls/otu_assembly_comparative_pipeline-/18S_censusparse/csv_outputs/eukcensus_18S_by_genus.csv',
    },
    '16S': {
        'phylum': '_project_pulls/otu_assembly_comparative_pipeline-/16S_censusparse/csv_outputs/eukcensus_16S_by_phylum.csv',
        'family': '_project_pulls/otu_assembly_comparative_pipeline-/16S_censusparse/csv_outputs/eukcensus_16S_by_family.csv',
        'genus': '_project_pulls/otu_assembly_comparative_pipeline-/16S_censusparse/csv_outputs/eukcensus_16S_by_genus.csv',
    }
}

def parse_lineage(lineage_str: str) -> list:
    """Parse semicolon-separated lineage into list of taxa."""
    if not lineage_str or lineage_str == 'NA':
        return []
    return [t.strip() for t in lineage_str.split(';') if t.strip()]

def build_index_from_csv(csv_file: str, rank: str) -> dict:
    """
    Parse a census CSV file and extract taxon→lineage mappings.
    Returns dict: {taxon_name: {'rank': rank, 'lineage': [...], 'otu_count': N}}
    """
    taxa = {}
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('Name_to_use', '')
            lineage_str = row.get('lineage', '')
            otu_count = int(row.get('otu_count', 0))
            
            if name and not name.endswith('.U.' + rank):  # Skip undefined taxa for now
                lineage = parse_lineage(lineage_str)
                taxa[name] = {
                    'rank': rank,
                    'lineage': lineage,
                    'otu_count': otu_count
                }
    
    return taxa

def build_ancestor_index(all_taxa: dict) -> dict:
    """
    Build reverse index: ancestor_name → list of descendant taxa.
    
    For each taxon, add it to the descendant list of every ancestor in its lineage.
    """
    ancestor_index = defaultdict(lambda: {'divisions': [], 'families': [], 'genera': []})
    
    for taxon_name, info in all_taxa.items():
        rank = info['rank']
        lineage = info['lineage']
        
        # Determine which list to add to based on rank
        if rank == 'division' or rank == 'phylum':
            list_key = 'divisions'
        elif rank == 'family':
            list_key = 'families'
        elif rank == 'genus':
            list_key = 'genera'
        else:
            continue
        
        # Add this taxon as a descendant of each ancestor in its lineage
        for ancestor in lineage:
            # Skip generic ancestors like "cellular organisms", "Eukaryota"
            if ancestor not in ['cellular organisms', 'root']:
                ancestor_index[ancestor][list_key].append({
                    'name': taxon_name,
                    'otu_count': info['otu_count']
                })
        
        # Also index by the taxon's own name (for direct matches)
        ancestor_index[taxon_name][list_key].append({
            'name': taxon_name,
            'otu_count': info['otu_count']
        })
    
    return dict(ancestor_index)

def main():
    parser = argparse.ArgumentParser(description='Build lineage index for taxonomic queries')
    parser.add_argument('--gene', required=True, choices=['18S', '16S'], help='Gene target')
    parser.add_argument('--output', help='Output JSON file (default: lineage_index_{gene}.json)')
    parser.add_argument('--project-root', help='Project root directory')
    
    args = parser.parse_args()
    
    # Determine project root
    if args.project_root:
        project_root = Path(args.project_root)
    else:
        project_root = Path(__file__).resolve().parent.parent.parent
    
    output_file = args.output or f'lineage_index_{args.gene}.json'
    
    print(f"🔍 Building lineage index for {args.gene}")
    print(f"   Project root: {project_root}")
    
    # Parse all census CSV files
    all_taxa = {}
    csv_paths = CENSUS_CSV_PATHS[args.gene]
    
    for rank, rel_path in csv_paths.items():
        csv_file = project_root / rel_path
        if csv_file.exists():
            print(f"   📖 Parsing {rank} CSV: {csv_file.name}")
            taxa = build_index_from_csv(str(csv_file), rank)
            all_taxa.update(taxa)
            print(f"      Found {len(taxa)} {rank}-level taxa")
        else:
            print(f"   ⚠️ File not found: {csv_file}")
    
    print(f"\n   Total taxa parsed: {len(all_taxa)}")

    # Build ancestor index
    print(f"\n🔗 Building ancestor→descendants index...")
    ancestor_index = build_ancestor_index(all_taxa)
    print(f"   Indexed {len(ancestor_index)} unique ancestor names")

    # Save to JSON
    output_path = project_root / 'primer_design_pipeline' / args.gene / output_file
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump({
            'gene': args.gene,
            'total_taxa': len(all_taxa),
            'total_ancestors': len(ancestor_index),
            'ancestor_index': ancestor_index
        }, f, indent=2)

    print(f"\n💾 Saved to: {output_path}")

    # Show some examples
    print(f"\n📋 Example lookups:")
    example_ancestors = ['Discoba', 'Metamonada', 'Amoebozoa', 'Stramenopiles', 'Alveolata']
    for ancestor in example_ancestors:
        if ancestor in ancestor_index:
            info = ancestor_index[ancestor]
            n_div = len(info['divisions'])
            n_fam = len(info['families'])
            n_gen = len(info['genera'])
            print(f"   {ancestor}: {n_div} divisions, {n_fam} families, {n_gen} genera")

if __name__ == "__main__":
    main()

