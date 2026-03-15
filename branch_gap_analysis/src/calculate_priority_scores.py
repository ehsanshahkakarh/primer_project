#!/usr/bin/env python3
"""
Calculate Primer Target Priority Scores

Identifies the best primer targets by combining:
1. Novelty Factor (NF) - how under-represented a taxon is in NCBI vs census
2. OTU abundance - how common it is in the environment

Formula: (log10(NF) + 1) × (log10(OTU count) + 1)

Higher scores = higher impact targets (both under-represented AND abundant)

Output: nf_abundance_priority_scores.csv with internal nodes (Order level and below)
ranked by total priority score of all descendant genera.
"""

import csv
import math
import re
from pathlib import Path
from dataclasses import dataclass, field

# ============================================================================
# Configuration
# ============================================================================

GENUS_PARSE_FILE = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv")
NF_FILE = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/final_merger/outputs/18s_ncbi_merged_genus.csv")
OUTPUT_DIR = Path("/Users/ehsankakarh/PROJA/primer_project/branch_gap_analysis/output/18s")

TARGET_RANKS = {'order', 'suborder', 'infraorder', 'family', 'subfamily', 'superfamily'}

SKIP_PATTERNS = [
    r'-lineage', r'_X{2,}', r'_MET[0-9]', r'^Novel-', r'^Group-',
    r'-Group-', r'^OLIGO[0-9]', r'^ARMOP[0-9]', r'-Clade-[0-9]',
    r'^Unclassified_', r'^unclassified', r'^epibiont$',
]

# ============================================================================
# Helper Functions
# ============================================================================

def should_skip(name: str) -> bool:
    for pattern in SKIP_PATTERNS:
        if re.search(pattern, name, re.IGNORECASE):
            return True
    return False

def calc_priority(nf: float, otu_count: int) -> float:
    if nf is None or nf <= 0 or otu_count is None or otu_count <= 0:
        return 0
    if math.isinf(nf):
        nf_component = 5
    else:
        nf_component = math.log10(nf) + 1
    otu_component = math.log10(otu_count) + 1
    return nf_component * otu_component

# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class TaxonNode:
    name: str
    rank: str
    children: dict = field(default_factory=dict)
    genera: list = field(default_factory=list)
    
    def get_or_create_child(self, name: str, rank: str) -> 'TaxonNode':
        if name not in self.children:
            self.children[name] = TaxonNode(name=name, rank=rank)
        return self.children[name]

@dataclass
class GenusInfo:
    name: str
    nf: float
    otu: int
    priority_score: float

# ============================================================================
# Main Functions
# ============================================================================

def load_genus_nf() -> dict:
    genus_nf = {}
    with open(NF_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('genus', '')
            try:
                nf = float(row.get('novelty_factor', '')) if row.get('novelty_factor', '') else None
            except:
                nf = None
            try:
                otu = int(row.get('census_otu_count', '')) if row.get('census_otu_count', '') else 0
            except:
                otu = 0
            if nf and nf >= 1.1:
                genus_nf[name] = {'nf': nf, 'otu': otu}
    return genus_nf

def build_tree(genus_nf: dict) -> TaxonNode:
    root = TaxonNode(name="root", rank="root")
    with open(GENUS_PARSE_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            genus_name = row.get('Name_to_use', '').strip()
            if not genus_name or should_skip(genus_name) or genus_name not in genus_nf:
                continue
            lineage = row.get('lineage', '').split(';')
            ranks = row.get('lineage_ranks', '').split(';')
            current = root
            for taxon, rank in zip(lineage, ranks):
                taxon, rank = taxon.strip(), rank.strip().lower()
                if not taxon or should_skip(taxon) or rank == 'original_name':
                    continue
                current = current.get_or_create_child(taxon, rank)
            nf_info = genus_nf[genus_name]
            score = calc_priority(nf_info['nf'], nf_info['otu'])
            current.genera.append(GenusInfo(genus_name, nf_info['nf'], nf_info['otu'], score))
    return root

def collect_all_genera(node: TaxonNode) -> list:
    all_genera = list(node.genera)
    for child in node.children.values():
        all_genera.extend(collect_all_genera(child))
    return all_genera

# Rank hierarchy (lower number = more specific/lower rank)
RANK_SPECIFICITY = {
    'subfamily': 1, 'family': 2, 'superfamily': 3,
    'infraorder': 4, 'suborder': 5, 'order': 6
}

def analyze_tree(root: TaxonNode) -> list:
    results = []

    def analyze_node(node, path=[]):
        current_path = path + [(node.name, node.rank)]
        if node.rank in TARGET_RANKS:
            all_genera = collect_all_genera(node)
            if all_genera:
                total_score = sum(g.priority_score for g in all_genera)
                total_otu = sum(g.otu for g in all_genera)
                all_genera.sort(key=lambda x: x.priority_score, reverse=True)
                path_str = ' > '.join([p[0] for p in current_path[-4:-1]])
                # Create a frozenset of genus names for deduplication
                genera_set = frozenset(g.name for g in all_genera)
                results.append({
                    'node': node.name,
                    'rank': node.rank,
                    'num_genera': len(all_genera),
                    'total_priority_score': round(total_score, 2),
                    'avg_priority_score': round(total_score / len(all_genera), 2),
                    'total_otu': total_otu,
                    'path': path_str,
                    'top_genera': ', '.join([f"{g.name}({g.priority_score:.1f})" for g in all_genera[:5]]),
                    'all_genera': ';'.join([g.name for g in all_genera]),  # All genera for sequence extraction
                    '_genera_set': genera_set  # For deduplication
                })
        for child in node.children.values():
            analyze_node(child, current_path)

    analyze_node(root)

    # Deduplicate: when multiple nodes cover exact same genera, keep the lowest rank
    seen_genera_sets = {}
    deduplicated = []

    for r in results:
        genera_set = r['_genera_set']
        rank_spec = RANK_SPECIFICITY.get(r['rank'], 99)

        if genera_set in seen_genera_sets:
            # We've seen this exact set before - check if this is more specific
            prev_idx, prev_spec = seen_genera_sets[genera_set]
            if rank_spec < prev_spec:
                # This is more specific, replace the previous one
                deduplicated[prev_idx] = r
                seen_genera_sets[genera_set] = (prev_idx, rank_spec)
        else:
            # New genera set
            seen_genera_sets[genera_set] = (len(deduplicated), rank_spec)
            deduplicated.append(r)

    # Remove the internal _genera_set field and sort
    for r in deduplicated:
        del r['_genera_set']

    deduplicated.sort(key=lambda x: x['total_priority_score'], reverse=True)
    return deduplicated

def main():
    print("=" * 70)
    print("Calculating Primer Target Priority Scores")
    print("Formula: (log10(NF) + 1) × (log10(OTU count) + 1)")
    print("=" * 70)

    print("\nLoading NF values...")
    genus_nf = load_genus_nf()
    print(f"  Loaded {len(genus_nf)} genera with NF >= 1.1")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # === Output 1: Individual genus priority scores ===
    print("\nCalculating individual genus priority scores...")
    genus_scores = []
    for name, info in genus_nf.items():
        if should_skip(name):
            continue
        score = calc_priority(info['nf'], info['otu'])
        if score > 0:
            genus_scores.append({
                'genus': name,
                'nf': info['nf'] if not math.isinf(info['nf']) else 'inf',
                'otu_count': info['otu'],
                'priority_score': round(score, 2)
            })

    genus_scores.sort(key=lambda x: x['priority_score'], reverse=True)

    genus_file = OUTPUT_DIR / "genus_priority_scores.csv"
    with open(genus_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['genus', 'nf', 'otu_count', 'priority_score'])
        writer.writeheader()
        writer.writerows(genus_scores)

    print(f"  Wrote {len(genus_scores)} genera to: {genus_file.name}")

    # === Output 2: Grouped primer targets ===
    print("\nBuilding taxonomic tree...")
    root = build_tree(genus_nf)

    print("\nAnalyzing internal nodes (Order level and below)...")
    results = analyze_tree(root)
    print(f"  Found {len(results)} unique primer target candidates (deduplicated)")

    output_file = OUTPUT_DIR / "nf_abundance_priority_scores.csv"

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'node', 'rank', 'num_genera', 'total_priority_score', 'avg_priority_score',
            'total_otu', 'path', 'top_genera', 'all_genera'
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"  Wrote grouped targets to: {output_file.name}")

    # === Summary ===
    print(f"\n" + "=" * 70)
    print("OUTPUT FILES:")
    print(f"  1. {genus_file.name} - Individual genus scores ({len(genus_scores)} genera)")
    print(f"  2. {output_file.name} - Grouped targets ({len(results)} groups)")
    print("=" * 70)

    print(f"\nTop 20 genera by priority score:")
    print(f"{'Genus':<30} {'NF':>10} {'OTU':>8} {'Score':>8}")
    print("-" * 60)
    for g in genus_scores[:20]:
        nf_str = str(g['nf'])[:8]
        print(f"{g['genus']:<30} {nf_str:>10} {g['otu_count']:>8} {g['priority_score']:>8}")

    print("\nDone!")

if __name__ == "__main__":
    main()
