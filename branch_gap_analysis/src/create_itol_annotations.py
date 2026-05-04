#!/usr/bin/env python3
"""
Create iTOL annotation files for the 18S genus tree.

Approach:
- **Genus colorstrip**: NF from `18s_final_merger_genus_name_novelty.csv` keyed by
  census genus (`Name_to_use`).
- **Family / division colorstrips**: EuKCensus `eukcensus_18S_by_genus.csv`
  lineages → collect family- and division-rank ancestors; for each genus, use the
  **most specific** ancestor name that appears in
  `18s_final_merger_family_name_novelty.csv` or
  `18s_final_merger_division_name_novelty.csv`. Genera sharing that ancestor get
  the same color (per tree node, `max` NF when several census rows collide).
- **Multibar**: census_size_count vs ncbi_genome_count from
  `18s_ncbi_merged_genus.csv`.

Tree node IDs match build_genus_tree.tree_to_newick (walk the built tree along
each census lineage). Multiple census rows mapping to the same node are merged:
  colorstrip / primer: max NF or any-hit; multibar: summed counts.

Uses the same skip rules as build_genus_tree.py (imported should_skip_taxon).
SHOW_INTERNAL 1 on strips/bars so internal nodes receive data when a census row
resolves there (e.g. misaligned lineage_ranks in source CSV).
"""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from collections.abc import Callable
from pathlib import Path

from build_genus_tree import (
    TaxonNode,
    build_tree_from_lineages,
    clean_lineage_name,
    clean_name,
    get_rank_suffix,
    should_skip_taxon,
)

# Layout: this file is primer_project/branch_gap_analysis/src/ → PROJA is four levels up.
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
OTU_PIPELINE_ROOT = REPO_ROOT / "GitHub" / "otu_assembly_comparative_pipeline-"

GENUS_PARSE_FILE = (
    OTU_PIPELINE_ROOT
    / "eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv"
)
TREE_FILE = (
    REPO_ROOT
    / "primer_project/branch_gap_analysis/output/18s/tree/eukcensus_18s_genus_tree.nwk"
)
OUTPUT_DIR = REPO_ROOT / "primer_project/branch_gap_analysis/output/18s/tree"

# Pre-extracted merger NF tables (name, novelty_factor) in the tree output dir
GENUS_NF_EXTRACT = OUTPUT_DIR / "18s_final_merger_genus_name_novelty.csv"
FAMILY_NF_EXTRACT = OUTPUT_DIR / "18s_final_merger_family_name_novelty.csv"
DIVISION_NF_EXTRACT = OUTPUT_DIR / "18s_final_merger_division_name_novelty.csv"

FINAL_MERGER_BASE = OTU_PIPELINE_ROOT / "final_merger" / "outputs"
FINAL_MERGER_FILES = {
    'division': FINAL_MERGER_BASE / "18s_ncbi_merged_division.csv",
    'family': FINAL_MERGER_BASE / "18s_ncbi_merged_family.csv",
    'genus': FINAL_MERGER_BASE / "18s_ncbi_merged_genus.csv",
}

# Ranks that count as "division" level (phylum/clade/kingdom)
DIVISION_RANKS = {'phylum', 'division', 'clade', 'kingdom', 'subkingdom'}
# Ranks that count as "family" level
FAMILY_RANKS = {'family', 'subfamily', 'superfamily', 'tribe', 'subtribe'}


def newick_label_for_node(node: TaxonNode) -> str:
    """Single-node label — must match build_genus_tree.tree_to_newick()."""
    suf = get_rank_suffix(node.rank)
    return f"{clean_name(node.name)}_{suf}"


def find_node_for_lineage(
    root: TaxonNode,
    names: list[str],
    ranks: list[str],
    taxids: list[str],
) -> TaxonNode | None:
    """
    Walk the built tree along one census lineage (same rules as integrate_genus_lineage).
    Returns the taxon node for this row (leaf or internal), or None if path is missing.
    """
    current = root
    for i, (name, rank, taxid) in enumerate(zip(names, ranks, taxids)):
        name = name.strip()
        rank = rank.strip()
        taxid = taxid.strip() if taxid else ''

        if not name:
            continue

        if i < len(names) - 1:
            name = clean_lineage_name(name)

        if name not in current.children:
            return None
        current = current.children[name]

    return current if current is not root else None


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


def _parse_novelty_factor(raw: str) -> float | None:
    s = raw.strip().lower()
    if s in ("inf", "infinity"):
        return float("inf")
    try:
        return float(s)
    except ValueError:
        return None


def load_name_novelty_extract(extract_file: Path, label: str) -> dict[str, float]:
    """Load ``{name: novelty_factor}`` from a two-column name/novelty_factor CSV."""
    if not extract_file.is_file():
        raise FileNotFoundError(
            f"{label} extract missing: {extract_file}\n"
            "Generate it under output/18s/tree first."
        )
    nf_data: dict[str, float] = {}
    with extract_file.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            name = (row.get("name") or "").strip()
            nf = _parse_novelty_factor(row.get("novelty_factor") or "")
            if name and nf is not None:
                nf_data[name] = nf
    return nf_data


def nf_from_lineage_names(
    ancestor_names_root_to_tip: list[str],
    nf_lookup: dict[str, float],
) -> float | None:
    """
    Match the **most tipward** lineage name that exists in ``nf_lookup``.
    ``ancestor_names`` are in EuKCensus order (root → genus); we try subfamily
    before family before superfamily, etc.
    """
    for name in reversed(ancestor_names_root_to_tip):
        name = name.strip()
        if not name:
            continue
        nf = nf_lookup.get(name)
        if nf is not None:
            return nf
    return None


def nf_merger_with_lineage_fallback(
    rank_chain: list[str],
    full_lineage_names: list[str],
    nf_lookup: dict[str, float],
) -> float | None:
    """
    Prefer ancestors at the requested rank(s) (e.g. family tribes), then any
    name along the census lineage so rows like ``Mammalia`` in the family
    merger still apply when ``Hominidae`` is absent.
    """
    nf = nf_from_lineage_names(rank_chain, nf_lookup)
    if nf is not None:
        return nf
    return nf_from_lineage_names(full_lineage_names, nf_lookup)


def build_genus_data(root: TaxonNode) -> tuple[dict, dict]:
    """
    Build mapping from each genus (Name_to_use) to its tree node id and ancestors.
    tree_node matches build_genus_tree.tree_to_newick labels for that row's taxon.

    Returns:
        Tuple of (genus_data dict, stats dict)
        genus_data: {genus_name: {'tree_node', 'divisions', 'families', 'lineage_names'}}
    """
    genus_data = {}
    stats = {'total': 0, 'skipped': 0, 'included': 0, 'path_missing': 0}

    with open(GENUS_PARSE_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            stats['total'] += 1
            genus_name = row.get('Name_to_use', '')
            lineage = row.get('lineage', '')
            lineage_ranks = row.get('lineage_ranks', '')
            lineage_taxids = row.get('lineage_taxids', '')

            if not genus_name or not lineage or not lineage_ranks:
                stats['skipped'] += 1
                continue

            # Skip environmental lineage taxa (keeps .U. unassigned)
            if should_skip_taxon(genus_name):
                stats['skipped'] += 1
                continue

            stats['included'] += 1
            names = [x.strip() for x in lineage.split(';')]
            ranks = [x.strip() for x in lineage_ranks.split(';')]
            taxids = (
                [x.strip() for x in lineage_taxids.split(';')]
                if lineage_taxids
                else [''] * len(names)
            )
            if len(taxids) < len(names):
                taxids = taxids + [''] * (len(names) - len(taxids))

            node = find_node_for_lineage(root, names, ranks, taxids)
            if node is None:
                stats['path_missing'] += 1
                continue
            tree_node = newick_label_for_node(node)

            # Collect ALL division-level and family-level ancestors
            divisions = []
            families = []

            for name, rank in zip(names, ranks):
                rank_l = rank.strip().lower()
                name_s = name.strip()

                if rank_l in DIVISION_RANKS and name_s:
                    divisions.append(name_s)
                if rank_l in FAMILY_RANKS and name_s:
                    families.append(name_s)

            genus_data[genus_name] = {
                'tree_node': tree_node,
                'divisions': divisions,
                'families': families,
                'lineage_names': [n.strip() for n in names if n.strip()],
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
    from Bio import Phylo
    from io import StringIO

    tree_content = TREE_FILE.read_text().rstrip()
    if not tree_content.endswith(';'):
        tree_content += ';'
    tree = Phylo.read(StringIO(tree_content), 'newick')
    terminal_ids = {n.name for n in tree.get_terminals() if n.name}

    # Extract all node names with suffixes
    import re
    node_pattern = r'([A-Za-z0-9_.-]+_[A-Za-z]+)'
    all_nodes = set(re.findall(node_pattern, TREE_FILE.read_text()))

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

        # Skip terminal leaves (all suffix types: G, ON, sG, …)
        if node in terminal_ids:
            continue

        name, _suffix = parts

        # Clean the name (replace underscores with spaces)
        clean_label = name.replace('_', ' ')

        # Format: node_id, label, position (0.5=centered), color, style, size, rotation
        lines.append(f"{node}\t{clean_label}\t0.5\t#444444\tbold\t1\t0")
        internal_count += 1

    output_file.write_text('\n'.join(lines))
    print(f"  Internal labels: {internal_count} nodes → {output_file.name}")


def export_tree_nodes_list(output_file: Path) -> None:
    """
    Write all Newick node IDs (same tokens as iTOL annotation files) to a TSV.
    Columns: node_id, is_leaf (yes=no descendants in tree = terminal clade label).
    """
    from Bio import Phylo
    from io import StringIO

    tree_content = TREE_FILE.read_text()
    node_pattern = r'([A-Za-z0-9_.-]+_[A-Za-z]+)'
    node_ids = sorted(set(re.findall(node_pattern, tree_content)))

    tc = tree_content.rstrip()
    if not tc.endswith(';'):
        tc += ';'
    tree = Phylo.read(StringIO(tc), 'newick')
    terminals = {n.name for n in tree.get_terminals() if n.name}

    with output_file.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['node_id', 'is_leaf'])
        for nid in node_ids:
            w.writerow([nid, 'yes' if nid in terminals else 'no'])

    n_leaf = sum(1 for n in node_ids if n in terminals)
    print(
        f"  Tree nodes: {len(node_ids)} total ({n_leaf} leaves) → {output_file.name}"
    )


def _format_nf_legend_value(nf: float) -> str:
    if math.isinf(nf):
        return "inf"
    return f"{nf:.1f}"


def create_colorstrip(
    genus_data: dict,
    output_file: Path,
    dataset_label: str,
    header_comment: str,
    get_nf_for_genus_row: Callable[[str, dict], float | None],
):
    """
    Create iTOL colorstrip — one row per tree node. Multiple census rows that
    map to the same node use max(NF) so iTOL gets a single color per node.
    """

    legend_shapes = "1\t1\t1\t1\t1\t1"
    legend_colors = "#228B22\t#FF8C00\t#DC143C\t#000080\t#4B0082\t#AAAAAA"
    legend_labels = "NF < 1.1\tNF 1.1-2\tNF 2-10\tNF > 10\tCensus Only\tNo Census Data"

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{dataset_label}",
        "COLOR\t#CCCCCC",
        "STRIP_WIDTH\t25",
        "MARGIN\t0",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t1",  # census rows can map to internal nodes
        "",
        header_comment,
        "LEGEND_TITLE\tNovelty Factor",
        f"LEGEND_SHAPES\t{legend_shapes}",
        f"LEGEND_COLORS\t{legend_colors}",
        f"LEGEND_LABELS\t{legend_labels}",
        "",
        "DATA"
    ]

    by_tree: dict[str, list[tuple[str, dict]]] = defaultdict(list)
    for genus_name, data in genus_data.items():
        by_tree[data['tree_node']].append((genus_name, data))

    matched = 0
    unmatched = 0

    for tree_node in sorted(by_tree):
        nfs: list[float] = []
        for genus_name, data in by_tree[tree_node]:
            nf = get_nf_for_genus_row(genus_name, data)
            if nf is not None:
                nfs.append(nf)

        if nfs:
            nf = max(nfs)
            color = get_nf_color(nf)
            lines.append(f"{tree_node}\t{color}\t{_format_nf_legend_value(nf)}")
            matched += 1
        else:
            lines.append(f"{tree_node}\t#AAAAAA\t")
            unmatched += 1

    output_file.write_text('\n'.join(lines))
    n_nodes = len(by_tree)
    short = output_file.name
    print(f"  {dataset_label}: {n_nodes} tree nodes — {matched} colored, {unmatched} grey → {short}")



def create_multibar(genus_data: dict, nf_genus_file: Path, output_file: Path):
    """
    Create iTOL DATASET_MULTIBAR: EuKCensus cluster size vs NCBI genome count per genus.

    Reads census_size_count and ncbi_genome_count from the genus-level final_merger file.
    Bar colors: EuKCensus orange (#FF8C00), NCBI black (#000000).
    """
    genus_counts = {}
    with open(nf_genus_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('genus', '')
            try:
                census_sz = int(row.get('census_size_count', 0) or 0)
            except (ValueError, TypeError):
                census_sz = 0
            try:
                ncbi_g = int(row.get('ncbi_genome_count', 0) or 0)
            except (ValueError, TypeError):
                ncbi_g = 0
            if name:
                genus_counts[name] = {
                    'census_size': census_sz,
                    'ncbi_genomes': ncbi_g,
                }

    c_euk = "#FF8C00"
    c_ncbi = "#000000"

    lines = [
        "DATASET_MULTIBAR",
        "SEPARATOR TAB",
        "DATASET_LABEL\tEuKCensus vs NCBI",
        "COLOR\t#555555",
        "",
        "# census_size_count (EuKCensus) and ncbi_genome_count (NCBI) per genus",
        f"FIELD_COLORS\t{c_euk}\t{c_ncbi}",
        "FIELD_LABELS\tEuKCensus (census size)\tNCBI (genome count)",
        "ALIGN_FIELDS\t0",
        "WIDTH\t200",
        "HEIGHT_FACTOR\t0.8",
        "BAR_SHIFT\t0",
        "MARGIN\t5",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t1",
        "",
        "LEGEND_TITLE\tCensus size & NCBI genomes",
        "LEGEND_SHAPES\t1\t1",
        f"LEGEND_COLORS\t{c_euk}\t{c_ncbi}",
        "LEGEND_LABELS\tEuKCensus (census size)\tNCBI (genome count)",
        "",
        "DATA",
    ]

    by_tree: dict[str, list[str]] = defaultdict(list)
    for genus_name, data in genus_data.items():
        by_tree[data['tree_node']].append(genus_name)

    matched = 0
    for tree_node in sorted(by_tree):
        v1 = v2 = 0
        for genus_name in by_tree[tree_node]:
            counts = genus_counts.get(
                genus_name, {'census_size': 0, 'ncbi_genomes': 0}
            )
            v1 += counts['census_size']
            v2 += counts['ncbi_genomes']
        if v1 > 0 or v2 > 0:
            lines.append(f"{tree_node}\t{v1}\t{v2}")
            matched += 1

    output_file.write_text('\n'.join(lines))
    print(f"  Multibar: {matched} tree nodes with data (from {len(by_tree)} nodes) → {output_file.name}")


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
        "SHOW_INTERNAL\t1",
        "",
        "# Genera highlighted in green = covered by a successful primer design",
        "LEGEND_TITLE\tPrimer Coverage",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#1B7A2B\t#E0E0E0",
        "LEGEND_LABELS\tPrimer Designed\tNo Primer",
        "",
        "DATA",
    ]

    by_tree: dict[str, list[str]] = defaultdict(list)
    for genus_name, data in genus_data.items():
        by_tree[data['tree_node']].append(genus_name)

    primer_count = 0
    no_primer_count = 0
    for tree_node in sorted(by_tree):
        if any(g in genera_with_primers for g in by_tree[tree_node]):
            lines.append(f"{tree_node}\t{PRIMER_COLOR}\tPrimer designed")
            primer_count += 1
        else:
            lines.append(f"{tree_node}\t{NO_PRIMER_COLOR}\t")
            no_primer_count += 1

    output_file.write_text('\n'.join(lines))
    print(f"  Primer highlight: {primer_count} tree nodes with primers, {no_primer_count} without → {output_file.name}")


def main():
    print("=" * 70)
    print("Creating iTOL Annotations")
    print("=" * 70)
    print("\nApproach:")
    print("  - Node IDs = Newick labels from build_genus_tree (per-row tree walk)")
    print("  - Division / family / genus NF strips (see module docstring)")
    print("  - Duplicate census rows → same node: max NF, summed bar counts")
    print("  - Skip patterns imported from build_genus_tree.should_skip_taxon")

    # Build reference tree (same object build_genus_tree.py uses)
    print(f"\nLoading taxonomy tree from: {GENUS_PARSE_FILE.name}")
    tree_root, tstats = build_tree_from_lineages(GENUS_PARSE_FILE)
    print(
        f"  Tree rows: {tstats['total']} | skipped: {tstats['skipped']} | "
        f"integrated lineages: {tstats['included']}"
    )

    print(f"\nBuilding genus → tree node ids from: {GENUS_PARSE_FILE.name}")
    genus_data, stats = build_genus_data(tree_root)
    print(f"  Total rows: {stats['total']}")
    print(f"  Skipped (*-lineage, _XX, etc.): {stats['skipped']}")
    print(f"  Included (after skip): {stats['included']}")
    print(f"  Rows with missing tree path: {stats['path_missing']}")
    print(f"  Annotation rows (genus_data): {len(genus_data)}")

    # Create node labels for ALL nodes (Order, Family, Clade, etc.)
    print("\nCreating node labels:")
    labels_file = OUTPUT_DIR / "18s_node_labels.txt"
    create_node_labels(labels_file)

    # Create internal labels (above genus level) for iTOL SHOW_INTERNAL display
    print("\nCreating internal clade labels:")
    internal_labels_file = OUTPUT_DIR / "18s_internal_labels.txt"
    create_internal_labels(internal_labels_file)

    print("\nExporting tree node list:")
    export_tree_nodes_list(OUTPUT_DIR / "18s_genus_tree_nodes.tsv")

    print("\nCreating colorstrips (merger NF extracts + EuKCensus lineages):")
    family_nf = load_name_novelty_extract(FAMILY_NF_EXTRACT, "Family NF")
    division_nf = load_name_novelty_extract(DIVISION_NF_EXTRACT, "Division NF")
    genus_nf_lookup = load_name_novelty_extract(GENUS_NF_EXTRACT, "Genus NF")
    print(
        f"  Loaded division={len(division_nf)} family={len(family_nf)} "
        f"genus={len(genus_nf_lookup)} NF rows"
    )

    create_colorstrip(
        genus_data,
        OUTPUT_DIR / "18s_division_colorstrip.txt",
        "18S Division (EuKCensus lineage → merger NF)",
        "# Division-rank chain first, then full lineage; tipward match in division novelty table",
        lambda _gn, d: nf_merger_with_lineage_fallback(
            d["divisions"], d["lineage_names"], division_nf
        ),
    )
    create_colorstrip(
        genus_data,
        OUTPUT_DIR / "18s_family_colorstrip.txt",
        "18S Family (EuKCensus lineage → merger NF)",
        "# Family-rank chain first, then full lineage; tipward match in family novelty table",
        lambda _gn, d: nf_merger_with_lineage_fallback(
            d["families"], d["lineage_names"], family_nf
        ),
    )
    create_colorstrip(
        genus_data,
        OUTPUT_DIR / "18s_genus_colorstrip.txt",
        "18S Genus (merger NF)",
        "# Genus Name_to_use ↔ genus novelty table",
        lambda gn, d: genus_nf_lookup.get(gn),
    )

    # Create multibar (EuKCensus census size vs NCBI genome count)
    print("\nCreating EuKCensus vs NCBI bar graph:")
    multibar_file = OUTPUT_DIR / "18s_eukcensus_ncbi_bars.txt"
    create_multibar(genus_data, FINAL_MERGER_FILES['genus'], multibar_file)

    # Create primer coverage highlight
    print("\nCreating primer coverage highlight:")
    primer_results_dir = REPO_ROOT / "primer_project" / "18S_subset" / "primer_results"
    priority_file = OUTPUT_DIR.parent / "nf_abundance_priority_scores.csv"
    primer_highlight_file = OUTPUT_DIR / "18s_primer_coverage.txt"
    create_primer_highlight(genus_data, priority_file, primer_results_dir, primer_highlight_file)

    print("\nDone! Upload 18s_genus_tree.nwk (or eukcensus_18s_genus_tree.nwk; same file) to iTOL.")


if __name__ == "__main__":
    main()

