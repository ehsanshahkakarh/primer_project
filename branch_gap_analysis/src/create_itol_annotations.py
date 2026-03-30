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


# File paths - relative to repo root
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent  # 00data/
GENUS_PARSE_FILE = REPO_ROOT / "00_gaps_taxonomic/00parse_database/eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv"
TREE_FILE = REPO_ROOT / "primer_project/branch_gap_analysis/output/18s/tree/18s_genus_tree.nwk"
OUTPUT_DIR = REPO_ROOT / "primer_project/branch_gap_analysis/output/18s/tree"

FINAL_MERGER_BASE = REPO_ROOT / "00_gaps_taxonomic/00parse_database/final_merger/outputs"
FINAL_MERGER_FILES = {
    'division': FINAL_MERGER_BASE / "18s_ncbi_merged_division.csv",
    'family': FINAL_MERGER_BASE / "18s_ncbi_merged_family.csv",
    'genus': FINAL_MERGER_BASE / "18s_ncbi_merged_genus.csv",
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



def create_multibar(genus_data: dict, nf_genus_file: Path, output_file: Path):
    """
    Create iTOL DATASET_MULTIBAR annotation with OTU counts and species counts
    for each genus leaf.

    Reads census_otu_count and ncbi_species_count from the genus-level NF file.
    """
    # Load OTU counts and species counts from genus NF file
    genus_counts = {}
    with open(nf_genus_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('genus', '')
            try:
                otu = int(row.get('census_otu_count', 0) or 0)
            except (ValueError, TypeError):
                otu = 0
            try:
                species = int(row.get('ncbi_species_count', 0) or 0)
            except (ValueError, TypeError):
                species = 0
            if name:
                genus_counts[name] = {'otu': otu, 'species': species}

    lines = [
        "DATASET_MULTIBAR",
        "SEPARATOR TAB",
        "DATASET_LABEL\tOTU & Species Counts",
        "COLOR\t#555555",
        "",
        "# Bar graph showing OTU counts and NCBI species counts per genus",
        "FIELD_COLORS\t#2196F3\t#FF9800",
        "FIELD_LABELS\tCensus OTU Count\tNCBI Species Count",
        "ALIGN_FIELDS\t0",
        "WIDTH\t200",
        "HEIGHT_FACTOR\t0.8",
        "BAR_SHIFT\t0",
        "MARGIN\t5",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t0",
        "",
        "LEGEND_TITLE\tOTU & Species Counts",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#2196F3\t#FF9800",
        "LEGEND_LABELS\tCensus OTU Count\tNCBI Species Count",
        "",
        "DATA",
    ]

    matched = 0
    for genus_name, data in sorted(genus_data.items()):
        tree_node = data['tree_node']
        counts = genus_counts.get(genus_name, {'otu': 0, 'species': 0})
        otu = counts['otu']
        species = counts['species']
        if otu > 0 or species > 0:
            lines.append(f"{tree_node}\t{otu}\t{species}")
            matched += 1

    output_file.write_text('\n'.join(lines))
    print(f"  Multibar: {matched} genera with data → {output_file.name}")


def create_primer_highlight(genus_data: dict, priority_file: Path,
                            primer_results_dir: Path, output_file: Path):
    """
    Create iTOL DATASET_COLORSTRIP that highlights genera covered by
    successful primer targets in bright/dark green.

    A primer target is successful if its directory contains primer_pairs.csv
    (as opposed to failure_report.txt).
    """
    import os
    import glob

    # Find all successful primer target names
    successful_targets = set()
    for pp_file in glob.glob(str(primer_results_dir / "*" / "primer_pairs.csv")):
        dir_name = os.path.basename(os.path.dirname(pp_file))
        # dir_name like "002_Arcellinida_order" → extract "Arcellinida"
        parts = dir_name.split('_', 1)  # ['002', 'Arcellinida_order']
        if len(parts) == 2:
            # Remove trailing _rank (order, family, subfamily, etc.)
            target_with_rank = parts[1]
            # Split from the right on underscore to strip the rank
            name_parts = target_with_rank.rsplit('_', 1)
            if len(name_parts) == 2:
                successful_targets.add(name_parts[0])
            else:
                successful_targets.add(target_with_rank)

    print(f"  Found {len(successful_targets)} successful primer targets")

    # Load priority scores to map target nodes → genera
    genera_with_primers = set()
    with open(priority_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            node = row.get('node', '')
            if node in successful_targets:
                all_genera_str = row.get('all_genera', '')
                genera = [g.strip() for g in all_genera_str.split(';') if g.strip()]
                genera_with_primers.update(genera)

    print(f"  {len(genera_with_primers)} genera covered by successful primers")

    # Color: bright green with darker hue
    PRIMER_COLOR = "#1B7A2B"  # Dark emerald green (bright but darker hue)
    NO_PRIMER_COLOR = "#E0E0E0"  # Light grey for non-primer genera

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tPrimer Coverage",
        "COLOR\t#1B7A2B",
        "STRIP_WIDTH\t20",
        "MARGIN\t2",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t0",
        "",
        "# Genera highlighted in green = covered by a successful primer design",
        "LEGEND_TITLE\tPrimer Coverage",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#1B7A2B\t#E0E0E0",
        "LEGEND_LABELS\tPrimer Designed\tNo Primer",
        "",
        "DATA",
    ]

    primer_count = 0
    no_primer_count = 0
    for genus_name, data in sorted(genus_data.items()):
        tree_node = data['tree_node']
        if genus_name in genera_with_primers:
            lines.append(f"{tree_node}\t{PRIMER_COLOR}\tPrimer designed")
            primer_count += 1
        else:
            lines.append(f"{tree_node}\t{NO_PRIMER_COLOR}\t")
            no_primer_count += 1

    output_file.write_text('\n'.join(lines))
    print(f"  Primer highlight: {primer_count} with primers, {no_primer_count} without → {output_file.name}")


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

    # Create multibar annotation (OTU counts + species counts)
    print("\nCreating OTU & Species count bar graph:")
    multibar_file = OUTPUT_DIR / "18s_otu_species_bars.txt"
    create_multibar(genus_data, FINAL_MERGER_FILES['genus'], multibar_file)

    # Create primer coverage highlight
    print("\nCreating primer coverage highlight:")
    primer_results_dir = REPO_ROOT / "primer_project" / "18S_subset" / "primer_results"
    priority_file = OUTPUT_DIR.parent / "nf_abundance_priority_scores.csv"
    primer_highlight_file = OUTPUT_DIR / "18s_primer_coverage.txt"
    create_primer_highlight(genus_data, priority_file, primer_results_dir, primer_highlight_file)

    print("\nDone! Upload 18s_genus_tree.nwk to iTOL, then add annotation files.")


if __name__ == "__main__":
    main()

