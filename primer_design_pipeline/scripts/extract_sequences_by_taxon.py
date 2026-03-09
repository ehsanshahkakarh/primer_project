#!/usr/bin/env python3
"""
Extract Sequences by Taxon from Cluster TSV

Generalized extraction script that can query by:
- division/phylum (direct match)
- family (direct match)
- genus (direct match or partial match for species)
- class/order (lineage-based extraction using ancestor index)

New: Can also filter by novelty status using merger CSV:
- census_only: Taxa only in census data (NF = infinity)
- matched: Taxa in both census and NCBI
- Filter by minimum novelty factor threshold

Usage:
    python extract_sequences_by_taxon.py --tsv <cluster_tsv> --taxon <taxon_name>
           --rank <rank> --output-dir <output_dir> [--fasta <main_fasta>]
           [--lineage-index <lineage_index.json>]

    # Extract all census-only taxa from merger CSV:
    python extract_sequences_by_taxon.py --tsv <cluster_tsv> --merger-csv <merger.csv>
           --match-status census_only --rank family --output-dir <output_dir>
"""

import csv
import ast
import sys
import json
import argparse
from pathlib import Path

csv.field_size_limit(sys.maxsize)


def get_taxa_from_merger_csv(merger_csv: str, match_status: str = None,
                              min_nf: float = None, rank: str = 'family') -> list:
    """
    Get list of taxa from merger CSV, optionally filtered by match_status or min novelty factor.

    Args:
        merger_csv: Path to the merger CSV file
        match_status: Filter by match_status ('census_only', 'matched', 'ncbi_only')
        min_nf: Minimum novelty factor threshold (use float('inf') for infinity)
        rank: The rank column to extract (family, genus, etc.)

    Returns:
        List of dicts with taxon info: {'name': ..., 'census_otu_count': ..., 'novelty_factor': ...}
    """
    taxa = []
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

            # Get taxon name (first column is usually the rank name)
            taxon_name = row.get(rank, list(row.values())[0])

            taxa.append({
                'name': taxon_name,
                'census_otu_count': int(row.get('census_otu_count', 0)),
                'ncbi_genome_count': int(row.get('ncbi_genome_count', 0)),
                'novelty_factor': row.get('novelty_factor', 'N/A'),
                'match_status': row.get('match_status', 'unknown')
            })

    return taxa

def extract_by_division(tsv_file: str, target_name: str) -> list:
    """Extract all sequences where division matches target."""
    sequences = []
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['division'] == target_name or row['division'].startswith(target_name):
                sequences.extend(process_cluster_row(row))
    return sequences

def extract_by_family(tsv_file: str, target_name: str) -> list:
    """Extract all sequences where family matches target."""
    sequences = []
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['family'] == target_name:
                sequences.extend(process_cluster_row(row))
    return sequences

def extract_by_genus(tsv_file: str, target_name: str, partial: bool = False) -> list:
    """Extract all sequences where genus matches target."""
    sequences = []
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if partial:
                if row['genus'].startswith(target_name):
                    sequences.extend(process_cluster_row(row))
            else:
                if row['genus'] == target_name:
                    sequences.extend(process_cluster_row(row))
    return sequences


def load_lineage_index(index_file: str) -> dict:
    """Load the ancestor→descendants index from JSON."""
    with open(index_file, 'r') as f:
        data = json.load(f)
    return data['ancestor_index']


def extract_by_lineage(tsv_file: str, ancestor_name: str, lineage_index: dict) -> list:
    """
    Extract sequences for all descendants of an ancestor taxon.

    This is used for class/order level queries where we need to find all
    divisions/families/genera that have this ancestor in their lineage.
    """
    # Get all descendant taxa from the lineage index
    descendants = lineage_index.get(ancestor_name, {})

    # Build sets of target names for each rank
    target_divisions = {d['name'] for d in descendants.get('divisions', [])}
    target_families = {f['name'] for f in descendants.get('families', [])}
    target_genera = {g['name'] for g in descendants.get('genera', [])}

    if not target_divisions and not target_families and not target_genera:
        print(f"   ⚠️ No descendants found for '{ancestor_name}' in lineage index")
        return []

    print(f"   📊 Lineage search for '{ancestor_name}':")
    print(f"      - {len(target_divisions)} divisions")
    print(f"      - {len(target_families)} families")
    print(f"      - {len(target_genera)} genera")

    # Extract sequences matching any descendant
    sequences = []
    matched_divisions = set()
    matched_families = set()
    matched_genera = set()

    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            match = False

            # Check division match
            if row['division'] in target_divisions:
                match = True
                matched_divisions.add(row['division'])

            # Check family match
            if row['family'] in target_families:
                match = True
                matched_families.add(row['family'])

            # Check genus match
            if row['genus'] in target_genera:
                match = True
                matched_genera.add(row['genus'])

            if match:
                sequences.extend(process_cluster_row(row))

    print(f"   ✓ Matched: {len(matched_divisions)} divisions, {len(matched_families)} families, {len(matched_genera)} genera")

    return sequences

def process_cluster_row(row: dict) -> list:
    """Process a cluster row and return list of sequence records."""
    sequences = []
    centroid = row['centroid']
    division = row['division']
    family = row['family']
    genus = row['genus']
    cluster_size = row['size']
    
    try:
        members = ast.literal_eval(row['members'])
    except (ValueError, SyntaxError):
        members = [centroid]
    
    # Add centroid
    sequences.append({
        'sequence_id': centroid,
        'division': division,
        'family': family,
        'genus': genus,
        'cluster_centroid': centroid,
        'cluster_size': cluster_size,
        'type': 'centroid'
    })
    
    # Add members
    for member in members:
        if member != centroid:
            sequences.append({
                'sequence_id': member,
                'division': division,
                'family': family,
                'genus': genus,
                'cluster_centroid': centroid,
                'cluster_size': cluster_size,
                'type': 'member'
            })
    
    return sequences

def write_sequence_tsv(sequences: list, output_file: str):
    """Write sequences to TSV file."""
    fieldnames = ['sequence_id', 'division', 'family', 'genus', 'cluster_centroid', 'cluster_size', 'type']
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(sequences)

def extract_fasta_sequences(fasta_file: str, sequence_ids: set, output_file: str):
    """Extract FASTA sequences matching the given IDs."""
    found = 0
    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout:
        write_current = False
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                write_current = seq_id in sequence_ids
                if write_current:
                    found += 1
            if write_current:
                fout.write(line)
    return found

def extract_single_taxon(args, taxon_name: str, output_dir: Path, taxon_label: str = None):
    """Extract sequences for a single taxon. Returns extraction stats dict."""
    taxon_label = taxon_label or taxon_name
    print(f"📖 Extracting sequences for: {taxon_label} (rank: {args.rank})")

    # Extract based on rank
    if args.rank == 'division':
        sequences = extract_by_division(args.tsv, taxon_name)
    elif args.rank == 'family':
        sequences = extract_by_family(args.tsv, taxon_name)
    elif args.rank == 'genus':
        sequences = extract_by_genus(args.tsv, taxon_name)
    elif args.rank == 'species':
        sequences = extract_by_genus(args.tsv, taxon_name, partial=True)
    elif args.rank in ['class', 'order', 'infraclass']:
        # Lineage-based extraction for class/order level taxa
        if not args.lineage_index:
            print(f"   ❌ Lineage index required for {args.rank} level queries")
            print(f"      Use: --lineage-index <lineage_index.json>")
            return None

        if not Path(args.lineage_index).exists():
            print(f"   ❌ Lineage index not found: {args.lineage_index}")
            print(f"      Run: python build_lineage_index.py --gene {args.gene}")
            return None

        lineage_index = load_lineage_index(args.lineage_index)
        sequences = extract_by_lineage(args.tsv, taxon_name, lineage_index)
    else:
        print(f"   ❌ Unsupported rank: {args.rank}")
        return None

    print(f"   ✓ Found {len(sequences)} sequences")

    if len(sequences) == 0:
        print(f"   ⚠️ No sequences found for {taxon_label}")
        return None

    # Create taxon-specific output directory
    taxon_output_dir = output_dir / taxon_name
    taxon_output_dir.mkdir(parents=True, exist_ok=True)

    # Write TSV
    tsv_output = taxon_output_dir / f"{taxon_name}_sequences.tsv"
    write_sequence_tsv(sequences, str(tsv_output))
    print(f"   💾 TSV saved: {tsv_output}")

    # Calculate centroid/member stats
    centroids = sum(1 for s in sequences if s['type'] == 'centroid')
    members = len(sequences) - centroids

    # Extract FASTA if provided
    found = 0
    if args.fasta:
        sequence_ids = {s['sequence_id'] for s in sequences}
        fasta_output = taxon_output_dir / f"{taxon_name}_{args.gene}.fasta"
        found = extract_fasta_sequences(args.fasta, sequence_ids, str(fasta_output))
        print(f"   💾 FASTA saved: {fasta_output} ({found} sequences)")

    # Summary
    print(f"   📊 Centroids: {centroids}, Members: {members}")
    print(f"   📊 Total in TSV: {len(sequences)}, Extracted to FASTA: {found}")

    # Save extraction stats to JSON for pipeline integration
    extraction_stats = {
        'taxon': taxon_label,
        'rank': args.rank,
        'gene': args.gene,
        'tsv_total': len(sequences),
        'centroids': centroids,
        'cluster_members': members,
        'fasta_sequences': found,
        'note': 'FASTA contains only centroid sequences; members are represented by their centroid at 97% OTU clustering'
    }

    stats_file = taxon_output_dir / 'extraction_stats.json'
    with open(stats_file, 'w') as f:
        json.dump(extraction_stats, f, indent=2)
    print(f"   💾 Stats saved: {stats_file}")

    return extraction_stats


def main():
    parser = argparse.ArgumentParser(description='Extract sequences by taxonomic group')
    parser.add_argument('--tsv', required=True, help='Cluster TSV file')
    parser.add_argument('--taxon', help='Taxon name (base name without rank suffix)')
    parser.add_argument('--rank', required=True,
                        choices=['division', 'family', 'genus', 'species', 'class', 'order', 'infraclass'],
                        help='Taxonomic rank to query')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--fasta', help='Main FASTA file to extract sequences from')
    parser.add_argument('--taxon-original', help='Original taxon name with suffix for labeling')
    parser.add_argument('--lineage-index', help='Lineage index JSON file for class/order queries')
    parser.add_argument('--gene', default='18S', help='Gene name for output file naming')

    # New merger CSV filtering options
    parser.add_argument('--merger-csv', help='Merger CSV file with novelty analysis results')
    parser.add_argument('--match-status', choices=['census_only', 'matched', 'ncbi_only'],
                        help='Filter taxa by match_status from merger CSV')
    parser.add_argument('--min-nf', type=float,
                        help='Minimum novelty factor threshold (use "inf" for infinity)')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Mode 1: Extract from merger CSV (batch mode for census_only, etc.)
    if args.merger_csv:
        if not Path(args.merger_csv).exists():
            print(f"❌ Merger CSV not found: {args.merger_csv}")
            sys.exit(1)

        print(f"📊 Reading taxa from merger CSV: {args.merger_csv}")
        print(f"   Filters: match_status={args.match_status}, min_nf={args.min_nf}")

        taxa = get_taxa_from_merger_csv(
            args.merger_csv,
            match_status=args.match_status,
            min_nf=args.min_nf,
            rank=args.rank
        )

        print(f"   ✓ Found {len(taxa)} taxa matching filters")

        if len(taxa) == 0:
            print("   ⚠️ No taxa found matching criteria")
            sys.exit(0)

        # Print taxa list
        print(f"\n{'='*60}")
        print(f"Taxa to process ({len(taxa)}):")
        print(f"{'='*60}")
        for i, t in enumerate(taxa[:20], 1):
            print(f"   {i:3}. {t['name']} (OTUs: {t['census_otu_count']}, NF: {t['novelty_factor']})")
        if len(taxa) > 20:
            print(f"   ... and {len(taxa) - 20} more")
        print(f"{'='*60}\n")

        # Process each taxon
        results = []
        for i, taxon_info in enumerate(taxa, 1):
            print(f"\n[{i}/{len(taxa)}] Processing: {taxon_info['name']}")
            print("-" * 40)

            stats = extract_single_taxon(
                args,
                taxon_name=taxon_info['name'],
                output_dir=output_dir,
                taxon_label=taxon_info['name']
            )

            if stats:
                stats['novelty_factor'] = taxon_info['novelty_factor']
                stats['match_status'] = taxon_info['match_status']
                stats['census_otu_count'] = taxon_info['census_otu_count']
                results.append(stats)

        # Save batch summary
        summary = {
            'merger_csv': args.merger_csv,
            'match_status_filter': args.match_status,
            'min_nf_filter': args.min_nf,
            'total_taxa_found': len(taxa),
            'taxa_with_sequences': len(results),
            'results': results
        }

        summary_file = output_dir / 'batch_extraction_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        print(f"\n{'='*60}")
        print(f"✅ Batch extraction complete!")
        print(f"   Taxa processed: {len(taxa)}")
        print(f"   Taxa with sequences: {len(results)}")
        print(f"   Summary saved: {summary_file}")
        print(f"{'='*60}")

    # Mode 2: Single taxon extraction (original behavior)
    else:
        if not args.taxon:
            print("❌ Either --taxon or --merger-csv is required")
            sys.exit(1)

        taxon_label = args.taxon_original or args.taxon

        # For single taxon, output directly to output_dir (original behavior)
        stats = extract_single_taxon(args, args.taxon, output_dir.parent, taxon_label)

        if stats is None:
            sys.exit(1)

if __name__ == "__main__":
    main()

