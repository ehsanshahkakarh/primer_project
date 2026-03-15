#!/usr/bin/env python3
"""
Calculate Division-Level Primer Target Priority Scores

Aggregates priority scores at the DIVISION level (e.g., Amoebozoa, Rhizaria, etc.)
for broader taxonomic primer design targets.

Uses the same priority formula as family/order level:
(log10(NF) + 1) × (log10(OTU count) + 1)

Output: nf_division_priority_scores.csv with divisions ranked by total priority score.
"""

import csv
import math
import re
from pathlib import Path
from collections import defaultdict

# ============================================================================
# Configuration
# ============================================================================

GENUS_PARSE_FILE = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv")
NF_FILE = Path("/Users/ehsankakarh/PROJA/GitHub/otu_assembly_comparative_pipeline-/final_merger/outputs/18s_ncbi_merged_genus.csv")
OUTPUT_DIR = Path("/Users/ehsankakarh/PROJA/primer_project/branch_gap_analysis/output/18s_division")

SKIP_PATTERNS = [
    r'-lineage', r'_X{2,}', r'_MET[0-9]', r'^Novel-', r'^Group-',
    r'-Group-', r'^OLIGO[0-9]', r'^ARMOP[0-9]', r'-Clade-[0-9]',
    r'^Unclassified_', r'^unclassified', r'^epibiont$',
]

# Overly broad terms to skip at division level
SKIP_DIVISIONS = {
    'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa',
    'Fungi', 'Viridiplantae', 'SAR', 'Amorphea', 'Obazoa',
    'root', 'life', 'Bilateria', 'Protostomia', 'Ecdysozoa',
}

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
# Main Functions
# ============================================================================

def load_genus_nf() -> dict:
    """Load NF values for genera with NF >= 1.1"""
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

def build_division_aggregates(genus_nf: dict) -> dict:
    """
    Build division-level aggregates from genus data.
    Returns: {division_name: {'genera': [...], 'total_score': X, 'total_otu': Y}}
    """
    divisions = defaultdict(lambda: {'genera': [], 'total_score': 0, 'total_otu': 0})
    
    with open(GENUS_PARSE_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            genus_name = row.get('Name_to_use', '').strip()
            if not genus_name or should_skip(genus_name) or genus_name not in genus_nf:
                continue
            
            # Get lineage and find division (usually early in the lineage)
            lineage = row.get('lineage', '').split(';')
            ranks = row.get('lineage_ranks', '').split(';')
            
            division = None
            for taxon, rank in zip(lineage, ranks):
                taxon, rank = taxon.strip(), rank.strip().lower()
                if rank == 'division' or rank == 'phylum':
                    division = taxon
                    break
            
            # Fallback: use first meaningful taxon after root levels
            if not division:
                for taxon, rank in zip(lineage, ranks):
                    taxon, rank = taxon.strip(), rank.strip().lower()
                    if rank not in ['superkingdom', 'kingdom', 'clade', 'no rank', ''] and not should_skip(taxon):
                        division = taxon
                        break
            
            if division and not should_skip(division) and division not in SKIP_DIVISIONS:
                nf_info = genus_nf[genus_name]
                score = calc_priority(nf_info['nf'], nf_info['otu'])
                divisions[division]['genera'].append({
                    'name': genus_name,
                    'nf': nf_info['nf'],
                    'otu': nf_info['otu'],
                    'score': score
                })
                divisions[division]['total_score'] += score
                divisions[division]['total_otu'] += nf_info['otu']

    return divisions


def main():
    print("=" * 70)
    print("Calculating DIVISION-Level Primer Target Priority Scores")
    print("Formula: (log10(NF) + 1) × (log10(OTU count) + 1)")
    print("=" * 70)

    print("\nLoading NF values...")
    genus_nf = load_genus_nf()
    print(f"  Loaded {len(genus_nf)} genera with NF >= 1.1")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\nBuilding division-level aggregates...")
    divisions = build_division_aggregates(genus_nf)
    print(f"  Found {len(divisions)} divisions with high-priority genera")

    # Convert to sorted list
    results = []
    for div_name, data in divisions.items():
        genera = data['genera']
        genera.sort(key=lambda x: x['score'], reverse=True)

        results.append({
            'division': div_name,
            'num_genera': len(genera),
            'total_priority_score': round(data['total_score'], 2),
            'avg_priority_score': round(data['total_score'] / len(genera), 2) if genera else 0,
            'total_otu': data['total_otu'],
            'top_genera': ', '.join([f"{g['name']}({g['score']:.1f})" for g in genera[:5]]),
            'all_genera': ';'.join([g['name'] for g in genera])
        })

    results.sort(key=lambda x: x['total_priority_score'], reverse=True)

    # Write output
    output_file = OUTPUT_DIR / "nf_division_priority_scores.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'division', 'num_genera', 'total_priority_score', 'avg_priority_score',
            'total_otu', 'top_genera', 'all_genera'
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"\n  Wrote {len(results)} divisions to: {output_file.name}")

    # Summary
    print(f"\n" + "=" * 70)
    print("TOP 20 DIVISIONS BY PRIORITY SCORE:")
    print("=" * 70)
    print(f"{'Division':<30} {'#Genera':>8} {'TotalOTU':>10} {'Score':>10}")
    print("-" * 60)
    for r in results[:20]:
        print(f"{r['division']:<30} {r['num_genera']:>8} {r['total_otu']:>10} {r['total_priority_score']:>10}")

    print("\nDone!")


if __name__ == "__main__":
    main()

