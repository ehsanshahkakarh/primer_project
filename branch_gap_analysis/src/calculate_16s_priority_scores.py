#!/usr/bin/env python3
"""
Calculate Primer Target Priority Scores for 16S (Archaea & Bacteria)

Separates by domain. Outputs to archaea/ and bacteria/ subdirectories.
Formula: (log10(NF) + 1) × (log10(OTU count) + 1)
"""

import csv
import math
import re
import sys
from pathlib import Path
from dataclasses import dataclass, field

# ============================================================================
# Configuration
# ============================================================================

GENUS_PARSE_FILE = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/eukcensus_parse/16S_censusparse/output/eukcensus16S_by_genus.csv")
NF_FILE = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/final_merger/outputs/16s_ncbi_merged_genus.csv")
OUTPUT_BASE = Path("/Users/ehsankakarh/PROJA/primer_project/branch_gap_analysis/output/16s")

TARGET_RANKS = {'order', 'suborder', 'infraorder', 'family', 'subfamily', 'superfamily', 'class', 'phylum'}

SKIP_PATTERNS = [
    r'-lineage', r'_X{2,}', r'_MET[0-9]', r'^Novel-', r'^Group-',
    r'-Group-', r'^OLIGO[0-9]', r'^ARMOP[0-9]', r'-Clade-[0-9]',
    r'^Unclassified_', r'^unclassified', r'^epibiont$',
]

# Domain to process (set via command line or process both)
DOMAIN = None  # Will be set in main()

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

def get_domain(lineage: str) -> str:
    """Extract domain from lineage."""
    if 'Archaea' in lineage:
        return 'archaea'
    elif 'Bacteria' in lineage:
        return 'bacteria'
    return None

def load_genus_nf(target_domain: str) -> dict:
    """Load NF values for a specific domain."""
    # First get domain mapping from census file
    genus_domain = {}
    with open(GENUS_PARSE_FILE, 'r') as f:
        for row in csv.DictReader(f):
            name = row.get('Name_to_use', '').strip()
            domain = get_domain(row.get('lineage', ''))
            if domain and name:
                genus_domain[name] = domain

    # Now load NF values, filtering by domain
    genus_nf = {}
    with open(NF_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('genus', '')
            if name not in genus_domain or genus_domain[name] != target_domain:
                continue
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
                if not taxon or should_skip(taxon) or rank in ('original_name', 'cellular root'):
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

def process_domain(domain: str):
    """Process a single domain (archaea or bacteria)."""
    output_dir = OUTPUT_BASE / domain
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  Processing {domain.upper()}")
    print(f"{'='*60}")

    print(f"\nLoading {domain} NF values...")
    genus_nf = load_genus_nf(domain)
    print(f"  Loaded {len(genus_nf)} genera with NF >= 1.1")

    if not genus_nf:
        print(f"  No data for {domain}, skipping...")
        return

    # === Output 1: Individual genus priority scores ===
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

    genus_file = output_dir / "genus_priority_scores.csv"
    with open(genus_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['genus', 'nf', 'otu_count', 'priority_score'])
        writer.writeheader()
        writer.writerows(genus_scores)
    print(f"  Wrote {len(genus_scores)} genera to: {genus_file.name}")

    # === Output 2: Grouped primer targets ===
    print("  Building taxonomic tree...")
    root = build_tree(genus_nf)
    results = analyze_tree(root)
    print(f"  Found {len(results)} unique primer target candidates")

    output_file = output_dir / "nf_abundance_priority_scores.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'node', 'rank', 'num_genera', 'total_priority_score', 'avg_priority_score',
            'total_otu', 'path', 'top_genera', 'all_genera'
        ])
        writer.writeheader()
        writer.writerows(results)
    print(f"  Wrote grouped targets to: {output_file.name}")

    # Top 10 summary
    print(f"\n  Top 10 {domain} genera:")
    for g in genus_scores[:10]:
        nf_str = str(g['nf'])[:6]
        print(f"    {g['genus']:<30} NF={nf_str:>6} OTU={g['otu_count']:>5} Score={g['priority_score']:.1f}")


def main():
    print("=" * 70)
    print("16S Priority Score Calculation (Archaea & Bacteria)")
    print("Formula: (log10(NF) + 1) × (log10(OTU count) + 1)")
    print("=" * 70)

    for domain in ['archaea', 'bacteria']:
        process_domain(domain)

    print("\n" + "=" * 70)
    print("Done! Output in:")
    print(f"  - {OUTPUT_BASE}/archaea/")
    print(f"  - {OUTPUT_BASE}/bacteria/")
    print("=" * 70)


if __name__ == "__main__":
    main()
