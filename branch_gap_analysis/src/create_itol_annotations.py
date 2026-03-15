#!/usr/bin/env python3
"""
Create iTOL annotation files for the 18S genus tree.

Approach:
- Division annotation: Color each GENUS leaf based on its DIVISION's NF value
- Family annotation: Color each GENUS leaf based on its FAMILY's NF value
- Genus annotation: Color each GENUS leaf based on its own NF value

Only genera (leaves) are colored - intermediate nodes are not included.
Filters out unassigned (.U.) and environmental lineage taxa to match tree.
"""

import csv
import re
from pathlib import Path


# File paths
GENUS_PARSE_FILE = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv")
TREE_FILE = Path("/Users/ehsankakarh/PROJA/primer_project/branch_gap_analysis/output/18s/tree/18s_genus_tree.nwk")
OUTPUT_DIR = Path("/Users/ehsankakarh/PROJA/primer_project/branch_gap_analysis/output/18s/tree")

FINAL_MERGER_FILES = {
    'division': Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/final_merger/outputs/18s_ncbi_merged_division.csv"),
    'family': Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/final_merger/outputs/18s_ncbi_merged_family.csv"),
    'genus': Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/final_merger/outputs/18s_ncbi_merged_genus.csv"),
}

# Ranks that count as "division" level (phylum/clade/kingdom)
DIVISION_RANKS = {'phylum', 'division', 'clade', 'kingdom', 'subkingdom'}
# Ranks that count as "family" level
FAMILY_RANKS = {'family', 'subfamily', 'superfamily', 'tribe', 'subtribe'}

# Patterns to SKIP (same as build_genus_tree.py - NOT .U. unassigned)
SKIP_PATTERNS = [
    r'-lineage',           # Environmental lineages: WIM80-lineage, Rhogostoma-lineage
    r'_X{2,}',             # Multiple X suffix: _XX, _XXX, _XXXX (unassigned)
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


def should_skip_taxon(name: str) -> bool:
    """Check if this taxon should be skipped based on naming patterns."""
    for pattern in SKIP_PATTERNS:
        if re.search(pattern, name, re.IGNORECASE):
            return True
    return False


def get_nf_color(nf: float) -> str:
    """Get color based on novelty factor."""
    if nf >= 9999:
        return "#4B0082"  # Indigo - Census only (infinite NF)
    elif nf > 10:
        return "#000080"  # Dark Blue - highly underrepresented
    elif nf > 2:
        return "#DC143C"  # Crimson
    elif nf > 1.1:
        return "#FF8C00"  # Orange
    else:
        return "#228B22"  # Green - well represented


def clean_name(name: str) -> str:
    """Clean name to match Newick format."""
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


def load_nf_lookup(merger_file: Path, level: str) -> dict:
    """Load novelty factor data from final_merger file."""
    nf_data = {}
    with open(merger_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            taxon_name = row.get(level, '')
            nf_str = row.get('novelty_factor', '')
            if taxon_name and nf_str:
                try:
                    nf = float(nf_str)
                    nf_data[taxon_name] = nf
                except ValueError:
                    continue
    return nf_data


def build_genus_data() -> tuple[dict, dict]:
    """
    Build mapping from each genus to its tree node name and ancestors.
    Filters out unassigned (.U.) and environmental lineage taxa.

    Returns:
        Tuple of (genus_data dict, stats dict)
        genus_data: {genus_name: {'tree_node': node_name, 'divisions': [...], 'families': [...]}}
    """
    genus_data = {}
    stats = {'total': 0, 'skipped': 0, 'included': 0}

    with open(GENUS_PARSE_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            stats['total'] += 1
            genus_name = row.get('Name_to_use', '')
            lineage = row.get('lineage', '')
            lineage_ranks = row.get('lineage_ranks', '')
            genus_rank = row.get('rank', 'genus')

            if not genus_name or not lineage or not lineage_ranks:
                stats['skipped'] += 1
                continue

            # Skip unassigned (.U.) and environmental lineage taxa
            if should_skip_taxon(genus_name):
                stats['skipped'] += 1
                continue

            stats['included'] += 1
            names = lineage.split(';')
            ranks = lineage_ranks.split(';')

            # Build the tree node name (same logic as build_genus_tree.py)
            clean = clean_name(genus_name)
            if genus_rank == 'genus':
                tree_node = f"{clean}_G"
            else:
                # Census-only taxa get _ON suffix
                tree_node = f"{clean}_ON"

            # Collect ALL division-level and family-level ancestors
            divisions = []
            families = []

            for name, rank in zip(names, ranks):
                rank = rank.strip().lower()
                name = name.strip()

                if rank in DIVISION_RANKS and name:
                    divisions.append(name)
                if rank in FAMILY_RANKS and name:
                    families.append(name)

            genus_data[genus_name] = {
                'tree_node': tree_node,
                'divisions': divisions,
                'families': families
            }

    return genus_data, stats


def create_node_labels(output_file: Path):
    """
    Create iTOL LABELS file for all nodes in the tree.
    Parses the Newick tree and creates readable labels for all nodes.
    """
    # Read tree to get all node names
    tree_content = TREE_FILE.read_text()

    # Extract all node names with suffixes (e.g., Alveolata_Cl, Homo_G)
    import re
    node_pattern = r'([A-Za-z0-9_.-]+_[A-Za-z]+)'
    all_nodes = set(re.findall(node_pattern, tree_content))

    # Rank suffix to full name mapping
    suffix_to_rank = {
        'D': 'Domain', 'K': 'Kingdom', 'sK': 'Subkingdom',
        'P': 'Phylum', 'sP': 'Subphylum', 'Dv': 'Division',
        'SC': 'Superclass', 'C': 'Class', 'sC': 'Subclass', 'iC': 'Infraclass',
        'SO': 'Superorder', 'O': 'Order', 'sO': 'Suborder', 'iO': 'Infraorder',
        'SF': 'Superfamily', 'F': 'Family', 'sF': 'Subfamily',
        'T': 'Tribe', 'sT': 'Subtribe',
        'G': 'Genus', 'sG': 'Subgenus',
        'S': 'Species', 'sS': 'Subspecies',
        'Cl': 'Clade', 'NR': '', 'R': '',
        'ON': ''  # Original name (census-only)
    }

    lines = [
        "DATASET_TEXT",
        "SEPARATOR TAB",
        "DATASET_LABEL\tNode Labels",
        "COLOR\t#000000",
        "",
        "# Labels for all tree nodes",
        "# Format: node_id <tab> label <tab> position <tab> color <tab> style <tab> size",
        "",
        "DATA"
    ]

    for node in sorted(all_nodes):
        # Split into name and suffix
        parts = node.rsplit('_', 1)
        if len(parts) == 2:
            name, suffix = parts
            rank_name = suffix_to_rank.get(suffix, '')

            # Clean the name (replace underscores with spaces)
            clean_label = name.replace('_', ' ')

            # Add rank in parentheses if it's informative
            if rank_name and rank_name not in ['', 'Genus']:
                label = f"{clean_label} ({rank_name})"
            else:
                label = clean_label
        else:
            label = node.replace('_', ' ')

        # Format: node_id, label, position, color, style, size_factor
        # position: -1=left, 0=center, 1=right
        lines.append(f"{node}\t{label}\t-1\t#000000\tnormal\t1")

    output_file.write_text('\n'.join(lines))
    print(f"  Labels: {len(all_nodes)} nodes → {output_file.name}")


def create_internal_labels(output_file: Path):
    """
    Create iTOL DATASET_TEXT file for INTERNAL nodes only (above genus level).
    Uses SHOW_INTERNAL 1 to display labels at branching points.
    Format matches old full_taxonomy_internal_labels.txt style.
    """
    # Read tree to get all node names
    tree_content = TREE_FILE.read_text()

    # Extract all node names with suffixes
    import re
    node_pattern = r'([A-Za-z0-9_.-]+_[A-Za-z]+)'
    all_nodes = set(re.findall(node_pattern, tree_content))

    lines = [
        "DATASET_TEXT",
        "SEPARATOR TAB",
        "DATASET_LABEL\tClade Labels",
        "COLOR\t#000000",
        "MARGIN\t0",
        "SHOW_INTERNAL\t1",
        "LABELS_BELOW\t0",
        "STRAIGHT_LABELS\t1",
        "SIZE_FACTOR\t1",
        "DATA"
    ]

    internal_count = 0
    for node in sorted(all_nodes):
        # Split into name and suffix
        parts = node.rsplit('_', 1)
        if len(parts) != 2:
            continue

        name, suffix = parts

        # Skip genera - we only want internal nodes
        if suffix == 'G':
            continue

        # Clean the name (replace underscores with spaces)
        clean_label = name.replace('_', ' ')

        # Format: node_id, label, position (0.5=centered), color, style, size, rotation
        lines.append(f"{node}\t{clean_label}\t0.5\t#444444\tbold\t1\t0")
        internal_count += 1

    output_file.write_text('\n'.join(lines))
    print(f"  Internal labels: {internal_count} nodes → {output_file.name}")


def create_colorstrip(level: str, genus_data: dict, nf_lookup: dict, output_file: Path):
    """
    Create iTOL colorstrip - color ONLY genera based on their ancestor's NF.

    genus_data: {genus_name: {'tree_node': ..., 'divisions': [...], 'families': [...]}}
    """

    legend_shapes = "1\t1\t1\t1\t1\t1"
    legend_colors = "#228B22\t#FF8C00\t#DC143C\t#000080\t#4B0082\t#AAAAAA"
    legend_labels = "NF < 1.1\tNF 1.1-2\tNF 2-10\tNF > 10\tCensus Only\tNo Census Data"

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t18S {level.title()} NF",
        "COLOR\t#CCCCCC",
        "STRIP_WIDTH\t25",
        "MARGIN\t0",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t0",  # Only color leaves (genera)
        "",
        f"# All genera colored by their {level}'s NF value",
        "LEGEND_TITLE\tNovelty Factor",
        f"LEGEND_SHAPES\t{legend_shapes}",
        f"LEGEND_COLORS\t{legend_colors}",
        f"LEGEND_LABELS\t{legend_labels}",
        "",
        "DATA"
    ]

    matched = 0
    unmatched = 0

    # Only iterate over genera (not all tree nodes)
    for genus_name, data in sorted(genus_data.items()):
        tree_node = data['tree_node']
        nf = None

        if level == 'genus':
            # For genus level, look up the genus directly
            if genus_name in nf_lookup:
                nf = nf_lookup[genus_name]
        else:
            # For division/family, try ALL ancestors and use first match
            ancestors_key = 'divisions' if level == 'division' else 'families'
            ancestors = data.get(ancestors_key, [])

            # Try each ancestor (most general to most specific)
            for ancestor in ancestors:
                if ancestor in nf_lookup:
                    nf = nf_lookup[ancestor]
                    break  # Use first match found

        if nf is not None:
            color = get_nf_color(nf)
            lines.append(f"{tree_node}\t{color}\t{nf:.1f}")
            matched += 1
        else:
            lines.append(f"{tree_node}\t#AAAAAA\t")
            unmatched += 1

    output_file.write_text('\n'.join(lines))
    print(f"  {level.title()}: {matched} colored, {unmatched} grey → {output_file.name}")



def main():
    print("=" * 70)
    print("Creating iTOL Annotations (Genera Only)")
    print("=" * 70)
    print("\nApproach:")
    print("  - Division: Color ALL genera by their DIVISION's NF")
    print("  - Family: Color ALL genera by their FAMILY's NF")
    print("  - Genus: Color ALL genera by their own NF")
    print("  - Only genus leaves are colored (not intermediate nodes)")
    print("  - Filters out unassigned (.U.) and *-lineage taxa")

    # Build genus data with tree node names and ancestors
    print(f"\nBuilding genus data from: {GENUS_PARSE_FILE.name}")
    genus_data, stats = build_genus_data()
    print(f"  Total rows: {stats['total']}")
    print(f"  Skipped (unassigned .U. or *-lineage): {stats['skipped']}")
    print(f"  Included (proper genera): {stats['included']}")

    # Create node labels for ALL nodes (Order, Family, Clade, etc.)
    print("\nCreating node labels:")
    labels_file = OUTPUT_DIR / "18s_node_labels.txt"
    create_node_labels(labels_file)

    # Create internal labels (above genus level) for iTOL SHOW_INTERNAL display
    print("\nCreating internal clade labels:")
    internal_labels_file = OUTPUT_DIR / "18s_internal_labels.txt"
    create_internal_labels(internal_labels_file)

    # Create colorstrip annotation for each level
    print("\nCreating colorstrip annotations:")
    for level, merger_file in FINAL_MERGER_FILES.items():
        nf_lookup = load_nf_lookup(merger_file, level)
        print(f"  Loaded {len(nf_lookup)} {level} NF values from {merger_file.name}")

        output_file = OUTPUT_DIR / f"18s_{level}_colorstrip.txt"
        create_colorstrip(level, genus_data, nf_lookup, output_file)

    print("\nDone! Upload 18s_genus_tree.nwk to iTOL, then add annotation files.")


if __name__ == "__main__":
    main()

