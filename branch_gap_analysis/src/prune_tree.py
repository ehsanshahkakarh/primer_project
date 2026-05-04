#!/usr/bin/env python3
"""
Prune the 18S genus tree to remove sections with no useful data.

Pruning criteria:
1. Remove genera with inf novelty factor (census-only, no NCBI data)
2. Remove genera with no census data (OTU count = 0)
3. Optionally keep genera if they are at multiple taxonomic ranks (not just genus)
4. Preserve tree structure - only prune terminal leaves unless entire clades are empty

This creates a cleaner tree focused on genera with actual genomic representation.

Usage:
    python prune_tree.py [--aggressive]
    
    Default: Keep genera with finite NF and census data
    --aggressive: Also remove genera with very high NF (>100)
"""

import argparse
import pandas as pd
from pathlib import Path
from Bio import Phylo
from io import StringIO
from collections import defaultdict

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
TREE_FILE = REPO_ROOT / "primer_project/branch_gap_analysis/output/18s/tree/eukcensus_18s_genus_tree.nwk"
OUTPUT_DIR = REPO_ROOT / "primer_project/branch_gap_analysis/output/18s/tree"
GENUS_MERGER_FILE = REPO_ROOT / "00_gaps_taxonomic/00parse_database/final_merger/outputs/18s_ncbi_merged_genus.csv"
GENUS_PARSE_FILE = REPO_ROOT / "00_gaps_taxonomic/00parse_database/eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv"


def load_genus_data():
    """Load genus-level novelty factor and census data."""
    merger_df = pd.read_csv(GENUS_MERGER_FILE)
    census_df = pd.read_csv(GENUS_PARSE_FILE)
    
    genus_info = {}
    for _, row in merger_df.iterrows():
        genus = row['genus']
        nf = row['novelty_factor']
        census_otu = row['census_otu_count']
        ncbi_genomes = row['ncbi_genome_count']
        
        genus_info[genus] = {
            'nf': nf,
            'census_otu': census_otu,
            'ncbi_genomes': ncbi_genomes,
            'has_data': census_otu > 0,
            'has_ncbi': ncbi_genomes > 0
        }
    
    return genus_info


def extract_genus_name(tree_label):
    """Extract genus name from tree label (e.g., 'Acanthamoeba_G' -> 'Acanthamoeba')."""
    # Remove suffix (_G, _ON, etc.)
    if '_' in tree_label:
        parts = tree_label.rsplit('_', 1)
        return parts[0].replace('_', '.')  # Convert underscores back for .U. names
    return tree_label


def should_prune_node(tree_label, genus_info, aggressive=False):
    """
    Determine if a node should be pruned.
    
    Args:
        tree_label: Tree node label (e.g., 'Acanthamoeba_G')
        genus_info: Dictionary of genus data
        aggressive: If True, also prune high NF (>100) genera
    
    Returns:
        (should_prune, reason)
    """
    # Extract genus name from tree label
    genus_name = extract_genus_name(tree_label)
    
    # If not in genus_info, we can't determine - keep it
    if genus_name not in genus_info:
        return False, "Not in genus data"
    
    info = genus_info[genus_name]
    
    # Prune if no census data
    if not info['has_data']:
        return True, "No census data"
    
    # Prune if inf NF (census-only, no NCBI representation)
    if info['nf'] == float('inf'):
        return True, "Census-only (inf NF)"
    
    # Aggressive: prune very high NF
    if aggressive and info['nf'] > 100:
        return True, f"Very high NF ({info['nf']:.1f})"
    
    return False, "Has data"


def prune_tree(tree, genus_info, aggressive=False):
    """
    Prune terminal nodes based on data criteria.
    
    Returns:
        (pruned_tree, stats)
    """
    stats = {
        'total_terminals': 0,
        'pruned': 0,
        'kept': 0,
        'reasons': defaultdict(int)
    }
    
    # Get all terminals
    terminals = tree.get_terminals()
    stats['total_terminals'] = len(terminals)
    
    to_prune = []
    for terminal in terminals:
        if terminal.name:
            should_prune, reason = should_prune_node(terminal.name, genus_info, aggressive)
            if should_prune:
                to_prune.append(terminal)
                stats['pruned'] += 1
                stats['reasons'][reason] += 1
            else:
                stats['kept'] += 1
    
    # Prune the terminals
    for terminal in to_prune:
        tree.prune(terminal)
    
    return tree, stats


def main():
    parser = argparse.ArgumentParser(description='Prune 18S genus tree based on data availability')
    parser.add_argument('--aggressive', action='store_true',
                        help='Also prune genera with very high NF (>100)')
    args = parser.parse_args()
    
    print("=" * 80)
    print("18S Genus Tree Pruning")
    print("=" * 80)
    print(f"Mode: {'AGGRESSIVE' if args.aggressive else 'STANDARD'}")
    print()
    
    # Load genus data
    print("Loading genus data...")
    genus_info = load_genus_data()
    print(f"  Loaded {len(genus_info)} genera")
    
    # Count data availability
    has_ncbi = sum(1 for g in genus_info.values() if g['has_ncbi'])
    census_only = sum(1 for g in genus_info.values() if g['has_data'] and not g['has_ncbi'])
    print(f"  Genera with NCBI data: {has_ncbi}")
    print(f"  Census-only genera (inf NF): {census_only}")
    
    # Load tree
    print(f"\nLoading tree: {TREE_FILE.name}")
    tree_content = TREE_FILE.read_text().strip()
    if not tree_content.endswith(';'):
        tree_content += ';'
    tree = Phylo.read(StringIO(tree_content), 'newick')

    original_terminal_count = len(tree.get_terminals())
    print(f"  Original tree: {original_terminal_count} terminal nodes")

    # Prune tree
    print(f"\nPruning tree...")
    pruned_tree, stats = prune_tree(tree, genus_info, args.aggressive)

    print(f"\nPruning Statistics:")
    print(f"  Total terminals: {stats['total_terminals']}")
    print(f"  Kept: {stats['kept']} ({stats['kept']/stats['total_terminals']*100:.1f}%)")
    print(f"  Pruned: {stats['pruned']} ({stats['pruned']/stats['total_terminals']*100:.1f}%)")
    print(f"\n  Pruning reasons:")
    for reason, count in sorted(stats['reasons'].items(), key=lambda x: -x[1]):
        print(f"    {reason}: {count}")

    # Save pruned tree
    suffix = '_aggressive' if args.aggressive else '_pruned'
    output_file = OUTPUT_DIR / f"eukcensus_18s_genus_tree{suffix}.nwk"

    print(f"\nSaving pruned tree: {output_file.name}")
    with open(output_file, 'w') as f:
        Phylo.write(pruned_tree, f, 'newick')

    print(f"\n✓ Pruned tree saved!")
    print(f"  File: {output_file}")
    print(f"  Terminals: {len(pruned_tree.get_terminals())}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Original tree: {original_terminal_count} genera")
    print(f"Pruned tree: {len(pruned_tree.get_terminals())} genera")
    print(f"Removed: {original_terminal_count - len(pruned_tree.get_terminals())} genera")
    print(f"\nThe pruned tree focuses on genera with actual NCBI genomic data,")
    print(f"removing census-only entries that would show as 'infinite' novelty factor.")


if __name__ == "__main__":
    main()
