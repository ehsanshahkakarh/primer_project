#!/usr/bin/env python3
"""
Analyze priority taxa and group them into umbrella categories for broader primer design.
"""

import csv
from pathlib import Path
from collections import defaultdict

PROJECT_DIR = Path(__file__).parent.parent.absolute()
PRIORITY_CSV = PROJECT_DIR.parent / "branch_gap_analysis" / "output" / "18s" / "nf_abundance_priority_scores.csv"


def parse_path_to_levels(path_str):
    """Parse taxonomic path into levels."""
    parts = [p.strip() for p in path_str.split(' > ')]
    return parts


def main():
    print("=" * 80)
    print("🎯 UMBRELLA PRIMER STRATEGY ANALYSIS")
    print("=" * 80)
    
    # Load priority data
    taxa = []
    with open(PRIORITY_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            path_parts = parse_path_to_levels(row['path'])
            taxa.append({
                'node': row['node'],
                'rank': row['rank'],
                'total_otu': int(row['total_otu']),
                'priority_score': float(row['total_priority_score']),
                'path': row['path'],
                'path_parts': path_parts,
                'root': path_parts[0] if path_parts else 'Unknown',
                'all_genera': row.get('all_genera', '')
            })
    
    print(f"\n📋 Loaded {len(taxa)} priority taxa\n")
    
    # Group by root taxonomic group
    by_root = defaultdict(list)
    for t in taxa:
        by_root[t['root']].append(t)
    
    print("-" * 80)
    print("📊 UMBRELLA GROUPS (by root taxonomic level)")
    print("-" * 80)
    
    umbrella_summary = []
    
    for root, members in sorted(by_root.items(), key=lambda x: -len(x[1])):
        total_otu = sum(m['total_otu'] for m in members)
        total_priority = sum(m['priority_score'] for m in members)
        all_genera = set()
        for m in members:
            all_genera.update(g.strip() for g in m['all_genera'].split(';') if g.strip())
        
        print(f"\n🌳 {root}")
        print(f"   Taxa count: {len(members)}")
        print(f"   Total OTUs: {total_otu:,}")
        print(f"   Combined priority: {total_priority:.1f}")
        print(f"   Unique genera: {len(all_genera)}")
        print(f"   Members: {', '.join(m['node'] for m in members[:5])}{'...' if len(members) > 5 else ''}")
        
        umbrella_summary.append({
            'umbrella_group': root,
            'child_taxa_count': len(members),
            'total_otu': total_otu,
            'total_priority': round(total_priority, 1),
            'unique_genera': len(all_genera),
            'child_taxa': ';'.join(m['node'] for m in members)
        })
    
    # Sort by priority
    umbrella_summary.sort(key=lambda x: -x['total_priority'])
    
    print("\n" + "=" * 80)
    print("🎯 RECOMMENDED UMBRELLA PRIMER TARGETS (by combined priority)")
    print("=" * 80)
    
    print(f"\n{'Rank':<4} {'Umbrella Group':<25} {'Taxa':<6} {'OTUs':<8} {'Priority':<10} {'Genera':<8}")
    print("-" * 70)
    
    for i, u in enumerate(umbrella_summary[:15], 1):
        print(f"{i:<4} {u['umbrella_group']:<25} {u['child_taxa_count']:<6} {u['total_otu']:<8} {u['total_priority']:<10} {u['unique_genera']:<8}")
    
    # Write umbrella summary
    output_file = PROJECT_DIR / "umbrella_primer_targets.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=umbrella_summary[0].keys())
        writer.writeheader()
        writer.writerows(umbrella_summary)
    
    print(f"\n\n✅ Umbrella targets saved to: {output_file}")
    
    # Strategy recommendations
    print("\n" + "=" * 80)
    print("💡 UMBRELLA PRIMER DESIGN STRATEGY")
    print("=" * 80)
    print("""
    APPROACH 1: Hierarchical Umbrella Design
    ----------------------------------------
    Instead of designing primers for each family/order separately,
    combine ALL sequences from taxa under the same umbrella group.
    
    Example: For "Amoebozoa" umbrella:
    - Combine: Euamoebida + Arcellinida + Hartmannellidae + Schizoplasmodiidae + ...
    - Find conserved regions across ALL Amoebozoa
    - Design primers that amplify the whole group
    
    APPROACH 2: Nested Primer Sets
    ------------------------------
    Design multiple primer pairs per umbrella:
    - Outer primers: Amplify entire umbrella group (broad)
    - Inner primers: Amplify specific sub-groups (specific)
    - Use nested PCR for increased specificity
    
    APPROACH 3: Degenerate Umbrella Primers
    ---------------------------------------
    Allow more degenerate positions (3-4 instead of 2)
    to capture sequence variation across the umbrella group.
""")


if __name__ == "__main__":
    main()

