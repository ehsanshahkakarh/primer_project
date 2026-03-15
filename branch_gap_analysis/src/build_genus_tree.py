#!/usr/bin/env python3
"""
Build a hierarchical tree from 18S census genus data.

This script:
1. Reads the eukcensus_18S_by_genus.csv file
2. Filters out unassigned (.U.) and environmental lineage taxa
3. Builds a true taxonomic tree using the lineage information
4. Outputs a Newick tree file for iTOL visualization
"""

import csv
import re
from pathlib import Path
from dataclasses import dataclass, field


# Patterns to SKIP (unassigned or environmental/census-only taxa)
SKIP_PATTERNS = [
    r'\.U\.',              # Unassigned: Eukaryota.U.genus, Fungi.U.genus
    r'-lineage',           # Environmental lineages: WIM80-lineage, Rhogostoma-lineage
    r'_X{2,}$',            # Multiple X suffix: _XX, _XXX, _XXXX (unassigned)
    r'_X{2,}_',            # Multiple X in middle: Kinetoplastea_XX_MET10
    r'_MET[0-9]',          # Metagenomic sequences: _MET10, _MET8
    r'^Novel-',            # Novel clades: Novel-Gran-6, Novel-clade-10
    r'^Group-',            # Environmental groups: Group-Te
    r'-Group-',            # Dino-Group-I, Dino-Group-II
    r'^OLIGO[0-9]',        # OLIGO designations: OLIGO2_X, OLIGO5
    r'^ARMOP[0-9]',        # ARMOP designations
    r'-Clade-[0-9]',       # Unnamed numbered clades: Clade-7, Clade-10
    r'^Unclassified_',     # Unclassified entries
    r'^unclassified',      # unclassified sequences/entries
    r'^epibiont$',         # Generic epibiont designation
]

# Patterns to clean from intermediate lineage names
CLEAN_PATTERNS = [
    (r'_X$', ''),          # Remove single _X suffix: Embryophyceae_X -> Embryophyceae
]


def should_skip_taxon(name: str) -> bool:
    """Check if this taxon should be skipped based on naming patterns."""
    for pattern in SKIP_PATTERNS:
        if re.search(pattern, name, re.IGNORECASE):
            return True
    return False


def clean_lineage_name(name: str) -> str:
    """Clean intermediate lineage names (remove _X suffix)."""
    result = name
    for pattern, replacement in CLEAN_PATTERNS:
        result = re.sub(pattern, replacement, result)
    return result


@dataclass
class TaxonNode:
    """Node in the taxonomic tree."""
    name: str
    rank: str
    taxid: str = ""
    otu_count: int = 0
    size_count: int = 0
    is_terminal: bool = False
    children: dict = field(default_factory=dict)

    def get_or_create_child(self, name: str, rank: str, taxid: str = "") -> 'TaxonNode':
        """Get existing child or create new one."""
        if name not in self.children:
            self.children[name] = TaxonNode(name=name, rank=rank, taxid=taxid)
        return self.children[name]


def clean_name(name: str) -> str:
    """Clean name for Newick format."""
    return (str(name)
            .replace(" ", "_")
            .replace("(", "")
            .replace(")", "")
            .replace(",", "")
            .replace(";", "")
            .replace(":", "_")
            .replace("'", "")
            .replace('"', "")
            .replace("[", "")
            .replace("]", ""))


def get_rank_suffix(rank: str) -> str:
    """Get short suffix for rank."""
    suffixes = {
        'domain': 'D', 'kingdom': 'K', 'subkingdom': 'sK',
        'phylum': 'P', 'subphylum': 'sP', 'division': 'Dv',
        'superclass': 'SC', 'class': 'C', 'subclass': 'sC', 'infraclass': 'iC',
        'superorder': 'SO', 'order': 'O', 'suborder': 'sO', 'infraorder': 'iO',
        'superfamily': 'SF', 'family': 'F', 'subfamily': 'sF',
        'tribe': 'T', 'subtribe': 'sT',
        'genus': 'G', 'subgenus': 'sG',
        'species': 'S', 'subspecies': 'sS',
        'clade': 'Cl', 'no rank': 'NR', 'cellular root': 'R',
        'original_name': 'ON'
    }
    return suffixes.get(rank, rank[0].upper() if rank else 'X')


def build_tree_from_lineages(input_file: Path) -> tuple[TaxonNode, dict]:
    """Build hierarchical tree from lineage data.

    Returns:
        Tuple of (root node, stats dict with skip/include counts)
    """
    root = TaxonNode(name="Life", rank="root")
    stats = {'total': 0, 'skipped': 0, 'included': 0}

    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)

        for row in reader:
            stats['total'] += 1
            genus_name = row.get('Name_to_use', '')
            lineage = row.get('lineage', '')
            lineage_ranks = row.get('lineage_ranks', '')
            lineage_taxids = row.get('lineage_taxids', '')
            otu_count = int(row.get('otu_count', 0) or 0)
            size_count = int(row.get('size_count', 0) or 0)

            if not lineage or not lineage_ranks:
                stats['skipped'] += 1
                continue

            # Skip unassigned (.U.) and environmental lineage taxa
            if should_skip_taxon(genus_name):
                stats['skipped'] += 1
                continue

            stats['included'] += 1

            # Parse lineage components
            names = lineage.split(';')
            ranks = lineage_ranks.split(';')
            taxids = lineage_taxids.split(';') if lineage_taxids else [''] * len(names)

            # Build path through tree, cleaning intermediate names
            current = root
            for i, (name, rank, taxid) in enumerate(zip(names, ranks, taxids)):
                name = name.strip()
                rank = rank.strip()
                taxid = taxid.strip() if taxid else ''

                if not name:
                    continue

                # Clean intermediate lineage names (remove _X, _XX)
                # but keep terminal genus name as-is for matching
                if i < len(names) - 1:
                    name = clean_lineage_name(name)

                current = current.get_or_create_child(name, rank, taxid)

            # Mark terminal node with census data
            if current != root:
                current.is_terminal = True
                current.otu_count += otu_count
                current.size_count += size_count

    return root, stats


def tree_to_newick(node: TaxonNode) -> str:
    """Convert tree to Newick format."""
    suffix = get_rank_suffix(node.rank)
    node_label = f"{clean_name(node.name)}_{suffix}"
    
    if not node.children:
        return node_label
    
    children_newick = ",".join(tree_to_newick(child) for child in node.children.values())
    return f"({children_newick}){node_label}"


def count_nodes(node: TaxonNode) -> tuple:
    """Count total and terminal nodes."""
    total = 1
    terminal = 1 if node.is_terminal else 0
    
    for child in node.children.values():
        t, term = count_nodes(child)
        total += t
        terminal += term
    
    return total, terminal


def print_tree_summary(node: TaxonNode, depth: int = 0, max_depth: int = 3):
    """Print tree structure summary."""
    if depth > max_depth:
        return
    
    indent = "  " * depth
    prefix = "├─ " if depth > 0 else ""
    terminal_marker = " *" if node.is_terminal else ""
    
    print(f"{indent}{prefix}{node.name} [{get_rank_suffix(node.rank)}]{terminal_marker} ({len(node.children)} children)")
    
    for child in list(node.children.values())[:5]:
        print_tree_summary(child, depth + 1, max_depth)


def main():
    """Main entry point."""
    input_file = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv")
    output_dir = Path("/Users/ehsankakarh/PROJA/primer_project/branch_gap_analysis/output/18s/tree")

    print("=" * 70)
    print("Building Genus Tree from 18S Census Data")
    print("=" * 70)

    # Build tree
    print(f"\nReading: {input_file}")
    root, stats = build_tree_from_lineages(input_file)

    print(f"\nFiltering stats:")
    print(f"  Total rows: {stats['total']}")
    print(f"  Skipped (unassigned .U. or *-lineage): {stats['skipped']}")
    print(f"  Included (proper genera): {stats['included']}")

    total, terminal = count_nodes(root)
    print(f"\nTree built: {total} total nodes, {terminal} terminal genera")

    # Print summary
    print(f"\nTree structure (top 3 levels):")
    print_tree_summary(root, max_depth=3)

    # Save Newick
    output_dir.mkdir(parents=True, exist_ok=True)
    newick = tree_to_newick(root) + ";"
    newick_file = output_dir / "18s_genus_tree.nwk"
    newick_file.write_text(newick)
    print(f"\nSaved: {newick_file}")


if __name__ == "__main__":
    main()

