#!/usr/bin/env python3
"""
Cross-reference primer validation results with EukCensus best-hit data.
Shows which reference sequences (best-hits) were used for primer design
and how they relate to PR2 taxonomy.
"""

import csv
import ast
import sys
from pathlib import Path
from collections import defaultdict

# Increase CSV field size limit
csv.field_size_limit(sys.maxsize)

# Paths
PROJECT_DIR = Path(__file__).parent.parent.absolute()
EUKCENSUS_CLUSTERS = Path("/Users/ehsankakarh/PROJA/eukcensus_metadata/new/eukcensus_2025_18S_clusters.tsv")
VALIDATION_RESULTS = PROJECT_DIR / "validation_results" / "pr2_validation_summary.csv"
OUTPUT_DIR = PROJECT_DIR / "validation_results"


def load_eukcensus_clusters():
    """Load eukcensus cluster data with best-hit information."""
    print("📥 Loading EukCensus cluster data...")
    
    clusters = {}
    with open(EUKCENSUS_CLUSTERS, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cluster_id = row['cluster_id']
            clusters[cluster_id] = {
                'centroid': row['centroid'],
                'best_hit': row.get('best-hit', ''),
                'best_hit_pident': row.get('best-hit-pident', ''),
                'best_hit_length': row.get('best-hit-length', ''),
                'domain': row.get('domain', ''),
                'supergroup': row.get('supergroup', ''),
                'division': row.get('division', ''),
                'order': row.get('order', ''),
                'family': row.get('family', ''),
                'genus': row.get('genus', ''),
            }
    
    print(f"   Loaded {len(clusters):,} clusters")
    return clusters


def find_clusters_for_taxon(clusters, taxon_name):
    """Find all clusters that match a given taxon name."""
    matches = []
    taxon_lower = taxon_name.lower()
    
    for cid, data in clusters.items():
        # Check if taxon appears in any taxonomy field
        for field in ['division', 'order', 'family', 'genus']:
            if taxon_lower in data.get(field, '').lower():
                matches.append((cid, data))
                break
    
    return matches


def analyze_best_hits(matches):
    """Analyze best-hit patterns for a set of clusters."""
    hit_sources = defaultdict(int)
    hit_organisms = defaultdict(int)
    
    for cid, data in matches:
        best_hit = data.get('best_hit', '')
        if best_hit:
            # Parse best-hit to extract source database
            if '_U;' in best_hit or '_U_' in best_hit:
                # SILVA/PR2 format: accession_U;18S_rRNA;...
                parts = best_hit.split(';')
                if len(parts) > 4:
                    # Try to extract organism from lineage
                    organism = parts[-1] if parts[-1] else parts[-2]
                    hit_organisms[organism[:30]] += 1
                    
            # Extract database prefix
            if best_hit.startswith('REF_'):
                hit_sources['REF (EukCensus)'] += 1
            elif '.' in best_hit.split('_')[0]:
                hit_sources['GenBank/SILVA'] += 1
            else:
                hit_sources['Other'] += 1
    
    return hit_sources, hit_organisms


def main():
    print("=" * 80)
    print("🔗 CROSS-REFERENCE: EukCensus Best-Hits vs PR2 Validation")
    print("=" * 80)
    
    # Load data
    clusters = load_eukcensus_clusters()
    
    # Load validation results
    print("\n📥 Loading validation results...")
    validation = []
    with open(VALIDATION_RESULTS, 'r') as f:
        reader = csv.DictReader(f)
        validation = list(reader)
    print(f"   Loaded {len(validation)} primer pairs")
    
    # Group by taxon
    by_taxon = defaultdict(list)
    for row in validation:
        by_taxon[row['target']].append(row)
    
    print(f"\n📋 Analyzing {len(by_taxon)} taxa\n")
    print("-" * 80)
    
    summary_rows = []
    
    for taxon, primers in sorted(by_taxon.items()):
        # Find matching clusters
        matches = find_clusters_for_taxon(clusters, taxon)
        
        print(f"\n🧬 {taxon}")
        print(f"   EukCensus clusters: {len(matches)}")
        
        if matches:
            # Analyze best-hit sources
            hit_sources, hit_organisms = analyze_best_hits(matches)
            
            print(f"   Best-hit sources: {dict(hit_sources)}")
            print(f"   Top organisms in best-hits:")
            for org, count in sorted(hit_organisms.items(), key=lambda x: -x[1])[:3]:
                print(f"      - {org}: {count}")
            
            # Sample best-hits
            print(f"   Sample best-hits:")
            for cid, data in matches[:3]:
                bh = data.get('best_hit', 'N/A')[:60]
                print(f"      - {bh}...")
        
        # Validation stats
        total_hits = sum(int(p['total_hits']) for p in primers)
        target_hits = sum(int(p['target_hits']) for p in primers)
        best_spec = max((float(p['specificity_pct']) for p in primers), default=0)
        
        print(f"   PR2 validation: {total_hits} total hits, {target_hits} target hits, best spec: {best_spec:.1f}%")
        
        summary_rows.append({
            'taxon': taxon,
            'eukcensus_clusters': len(matches),
            'pr2_total_hits': total_hits,
            'pr2_target_hits': target_hits,
            'best_specificity': best_spec,
            'hit_sources': str(dict(hit_sources)) if matches else 'N/A'
        })
    
    # Write summary
    output_file = OUTPUT_DIR / "eukcensus_pr2_crossref.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)
    
    print(f"\n\n✅ Cross-reference saved to: {output_file}")


if __name__ == "__main__":
    main()

