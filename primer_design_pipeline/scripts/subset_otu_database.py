#!/usr/bin/env python3
"""
Subset OTU Database by Novelty Analysis

Creates a filtered subset of the OTU cluster TSV containing only rows
that match taxa from the merger CSV based on novelty criteria.

This creates a smaller, focused database for downstream analysis.

Usage:
    python subset_otu_database.py \
        --cluster-tsv <eukcensus_clusters.tsv> \
        --merger-csv <merger_analysis.csv> \
        --match-status census_only \
        --rank family \
        --output-dir metadata/

Output:
    metadata/
    ├── subset_clusters.tsv          # Filtered cluster data
    ├── subset_taxa_list.txt         # List of included taxa
    └── subset_summary.json          # Summary statistics
"""

import csv
import sys
import json
import argparse
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(sys.maxsize)


def get_taxa_from_merger_csv(merger_csv: str, match_status: str = None,
                              min_nf: float = None, rank: str = 'family') -> dict:
    """
    Get taxa from merger CSV with their novelty info.
    Returns dict: {taxon_name: {census_otu_count, novelty_factor, match_status}}
    """
    taxa = {}
    with open(merger_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Check match_status filter
            if match_status and row.get('match_status') != match_status:
                continue

            # Check novelty factor filter
            if min_nf is not None:
                nf_str = row.get('novelty_factor', '0')
                if nf_str == 'inf':
                    nf = float('inf')
                else:
                    try:
                        nf = float(nf_str)
                    except ValueError:
                        nf = 0
                if nf < min_nf:
                    continue

            # Get taxon name (first column)
            taxon_name = row.get(rank, list(row.values())[0])

            taxa[taxon_name] = {
                'census_otu_count': int(row.get('census_otu_count', 0)),
                'ncbi_genome_count': int(row.get('ncbi_genome_count', 0)),
                'novelty_factor': row.get('novelty_factor', 'N/A'),
                'match_status': row.get('match_status', 'unknown')
            }

    return taxa


def subset_cluster_tsv(cluster_tsv: str, target_taxa: set, rank: str) -> tuple:
    """
    Extract rows from cluster TSV where the rank column matches target taxa.
    Returns (matching_rows, stats_dict)
    """
    matching_rows = []
    taxa_found = defaultdict(int)
    total_sequences = 0

    # Map rank to column name
    rank_column = rank if rank != 'phylum' else 'division'

    with open(cluster_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames

        for row in reader:
            taxon_value = row.get(rank_column, '')

            if taxon_value in target_taxa:
                matching_rows.append(row)
                taxa_found[taxon_value] += 1
                total_sequences += int(row.get('size', 1))

    stats = {
        'total_rows': len(matching_rows),
        'total_sequences': total_sequences,
        'taxa_found': len(taxa_found),
        'taxa_counts': dict(taxa_found)
    }

    return matching_rows, fieldnames, stats


def write_subset_tsv(rows: list, fieldnames: list, output_file: str):
    """Write subset rows to TSV file."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description='Create subset OTU database based on novelty analysis')
    parser.add_argument('--cluster-tsv', required=True,
                        help='Full cluster TSV file (eukcensus_*.clusters.*.tsv)')
    parser.add_argument('--merger-csv', required=True,
                        help='Merger CSV with novelty analysis results')
    parser.add_argument('--rank', required=True,
                        choices=['division', 'phylum', 'family', 'genus'],
                        help='Taxonomic rank to filter by')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for subset files')
    parser.add_argument('--match-status',
                        choices=['census_only', 'matched', 'ncbi_only'],
                        help='Filter by match_status')
    parser.add_argument('--min-nf', type=float,
                        help='Minimum novelty factor threshold')
    parser.add_argument('--prefix', default='subset',
                        help='Prefix for output files (default: subset)')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Get target taxa from merger CSV
    print(f"📊 Reading merger CSV: {args.merger_csv}")
    print(f"   Filters: match_status={args.match_status}, min_nf={args.min_nf}")

    taxa_info = get_taxa_from_merger_csv(
        args.merger_csv,
        match_status=args.match_status,
        min_nf=args.min_nf,
        rank=args.rank
    )

    print(f"   ✓ Found {len(taxa_info)} taxa matching criteria")

    if len(taxa_info) == 0:
        print("   ⚠️ No taxa found. Exiting.")
        sys.exit(0)

    # Step 2: Subset the cluster TSV
    print(f"\n📖 Subsetting cluster TSV: {args.cluster_tsv}")
    target_taxa = set(taxa_info.keys())

    rows, fieldnames, stats = subset_cluster_tsv(
        args.cluster_tsv, target_taxa, args.rank)

    print(f"   ✓ Found {stats['total_rows']} matching rows")
    print(f"   ✓ {stats['taxa_found']} taxa represented")
    print(f"   ✓ {stats['total_sequences']} total sequences (including cluster members)")

    # Step 3: Write output files
    # 3a: Subset TSV
    subset_tsv_file = output_dir / f"{args.prefix}_clusters.tsv"
    write_subset_tsv(rows, fieldnames, str(subset_tsv_file))
    print(f"\n💾 Saved: {subset_tsv_file}")

    # 3b: Summary JSON
    summary = {
        'source_cluster_tsv': args.cluster_tsv,
        'source_merger_csv': args.merger_csv,
        'filters': {
            'rank': args.rank,
            'match_status': args.match_status,
            'min_nf': args.min_nf
        },
        'results': {
            'taxa_in_merger': len(taxa_info),
            'taxa_found_in_clusters': stats['taxa_found'],
            'total_cluster_rows': stats['total_rows'],
            'total_sequences': stats['total_sequences']
        },
        'taxa_details': {
            taxon: {
                **taxa_info[taxon],
                'clusters_found': stats['taxa_counts'].get(taxon, 0)
            }
            for taxon in taxa_info
        },
        'output_files': {
            'subset_clusters': str(subset_tsv_file)
        }
    }

    summary_file = output_dir / f"{args.prefix}_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"💾 Saved: {summary_file}")

    # Print summary table
    print(f"\n{'='*60}")
    print(f"📊 SUBSET DATABASE SUMMARY")
    print(f"{'='*60}")
    print(f"   Filter: {args.rank} with match_status={args.match_status}")
    print(f"   Taxa in merger CSV: {len(taxa_info)}")
    print(f"   Taxa found in clusters: {stats['taxa_found']}")
    print(f"   Total cluster rows: {stats['total_rows']}")
    print(f"   Total sequences: {stats['total_sequences']}")
    print(f"{'='*60}")

    # Show top taxa by cluster count
    print(f"\n📈 Top 10 taxa by cluster count:")
    sorted_taxa = sorted(stats['taxa_counts'].items(), key=lambda x: x[1], reverse=True)[:10]
    for i, (taxon, count) in enumerate(sorted_taxa, 1):
        nf = taxa_info[taxon]['novelty_factor']
        print(f"   {i:2}. {taxon}: {count} clusters (NF={nf})")

    print(f"\n✅ Subset database created successfully!")


if __name__ == "__main__":
    main()

