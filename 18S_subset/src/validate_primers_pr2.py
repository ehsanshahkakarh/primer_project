#!/usr/bin/env python3
"""
Validate primers against PR2 database using in-silico PCR.

Tests designed primers against the PR2 18S rRNA database to check:
1. Specificity: How many sequences match each primer pair
2. Coverage: What taxonomic groups are amplified
3. Off-target hits: What non-target taxa would be amplified

Usage:
    python src/validate_primers_pr2.py --taxon 001_Euamoebida_order --top-n 10
    python src/validate_primers_pr2.py --all --top-n 5
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict
from Bio import SeqIO
from Bio.Seq import Seq

# Paths
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_DIR = SCRIPT_DIR.parent
RESULTS_DIR = PROJECT_DIR / "primer_results"
PR2_FASTA = Path("/Users/ehsankakarh/PROJA/eukcensus_metadata/primer_amoebozoa/PR2/pr2_version_5.1.1_SSU_taxo_long.fasta")


def reverse_complement(seq: str) -> str:
    """Get reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def find_primer_binding_exact(seq: str, primer: str) -> List[int]:
    """Find all exact match positions where primer binds (fast)."""
    positions = []
    primer_upper = primer.upper()
    seq_upper = seq.upper()
    start = 0
    while True:
        pos = seq_upper.find(primer_upper, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions


def find_primer_binding(seq: str, primer: str, max_mismatches: int = 2) -> List[int]:
    """Find all positions where primer binds (with allowed mismatches)."""
    # Fast path: try exact match first
    if max_mismatches == 0:
        return find_primer_binding_exact(seq, primer)

    positions = []
    primer_len = len(primer)
    seq_upper = seq.upper()
    primer_upper = primer.upper()

    # Optimization: use first 8bp as anchor, then check full primer
    anchor = primer_upper[:8]
    start = 0
    while True:
        pos = seq_upper.find(anchor, start)
        if pos == -1:
            break
        if pos + primer_len <= len(seq_upper):
            region = seq_upper[pos:pos + primer_len]
            mismatches = sum(1 for a, b in zip(region, primer_upper) if a != b and a != 'N' and b != 'N')
            if mismatches <= max_mismatches:
                positions.append(pos)
        start = pos + 1

    return positions


def check_amplicon(seq: str, fwd_primer: str, rev_primer: str, 
                   min_size: int = 100, max_size: int = 1000,
                   max_mismatches: int = 2) -> Tuple[bool, int]:
    """
    Check if a sequence would be amplified by the primer pair.
    
    Returns: (would_amplify, amplicon_size)
    """
    # Find forward primer binding sites
    fwd_sites = find_primer_binding(seq, fwd_primer, max_mismatches)
    
    # Find reverse primer binding sites (need reverse complement on opposite strand)
    rev_rc = reverse_complement(rev_primer)
    rev_sites = find_primer_binding(seq, rev_rc, max_mismatches)
    
    # Check for valid amplicons
    for fwd_pos in fwd_sites:
        for rev_pos in rev_sites:
            amplicon_size = rev_pos + len(rev_primer) - fwd_pos
            if min_size <= amplicon_size <= max_size:
                return True, amplicon_size
    
    return False, 0


def parse_pr2_header(header: str) -> Dict[str, str]:
    """Parse PR2 FASTA header to extract taxonomy."""
    # Format: >AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|...
    parts = header.split('|')
    taxonomy = {
        'accession': parts[0] if len(parts) > 0 else '',
        'gene': parts[1] if len(parts) > 1 else '',
        'compartment': parts[2] if len(parts) > 2 else '',
    }
    # Taxonomy starts at position 4 (after empty field)
    if len(parts) > 4:
        tax_parts = parts[4:]
        taxonomy['domain'] = tax_parts[0] if len(tax_parts) > 0 else ''
        taxonomy['supergroup'] = tax_parts[1] if len(tax_parts) > 1 else ''
        taxonomy['phylum'] = tax_parts[2] if len(tax_parts) > 2 else ''
        taxonomy['class'] = tax_parts[3] if len(tax_parts) > 3 else ''
        taxonomy['order'] = tax_parts[4] if len(tax_parts) > 4 else ''
        taxonomy['family'] = tax_parts[5] if len(tax_parts) > 5 else ''
        taxonomy['genus'] = tax_parts[6] if len(tax_parts) > 6 else ''
        taxonomy['species'] = tax_parts[7] if len(tax_parts) > 7 else ''
        taxonomy['full_lineage'] = '|'.join(tax_parts)
    return taxonomy


def load_top_primers(taxon_dir: Path, top_n: int = 10) -> List[Tuple[str, str, str]]:
    """Load top N primer pairs from a taxon results directory."""
    fasta_file = taxon_dir / "primer_pairs.fasta"
    if not fasta_file.exists():
        return []
    
    primers = []
    current_pair = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        pair_id = record.id.rsplit('_', 1)[0]  # e.g., "pair_001"
        direction = record.id.rsplit('_', 1)[1]  # "forward" or "reverse"
        
        if pair_id not in current_pair:
            current_pair[pair_id] = {}
        current_pair[pair_id][direction] = str(record.seq)
        
        if len(current_pair) > top_n:
            break
    
    # Convert to list of tuples
    for pair_id in sorted(current_pair.keys())[:top_n]:
        if 'forward' in current_pair[pair_id] and 'reverse' in current_pair[pair_id]:
            primers.append((
                pair_id,
                current_pair[pair_id]['forward'],
                current_pair[pair_id]['reverse']
            ))
    
    return primers


def validate_primers(taxon: str, primers: List[Tuple[str, str, str]], 
                    pr2_seqs: Dict[str, Tuple[str, Dict]], 
                    max_mismatches: int = 2) -> Dict:
    """Validate primers against PR2 sequences."""
    results = {
        'taxon': taxon,
        'primers_tested': len(primers),
        'primer_results': []
    }
    
    for pair_id, fwd, rev in primers:
        hits = []
        taxonomy_hits = defaultdict(int)
        
        for seq_id, (sequence, taxonomy) in pr2_seqs.items():
            amplifies, amp_size = check_amplicon(
                sequence, fwd, rev, 
                min_size=100, max_size=1000,
                max_mismatches=max_mismatches
            )
            if amplifies:
                hits.append({
                    'seq_id': seq_id,
                    'amplicon_size': amp_size,
                    'taxonomy': taxonomy.get('full_lineage', '')
                })
                # Count by phylum/class
                key = f"{taxonomy.get('phylum', 'Unknown')}|{taxonomy.get('class', 'Unknown')}"
                taxonomy_hits[key] += 1
        
        results['primer_results'].append({
            'pair_id': pair_id,
            'forward': fwd,
            'reverse': rev,
            'total_hits': len(hits),
            'taxonomy_breakdown': dict(taxonomy_hits),
            'sample_hits': hits[:5]  # First 5 hits as examples
        })
    
    return results


def load_pr2_by_taxon(pr2_fasta: Path, target_taxon: str, limit: int = 5000) -> Dict[str, Tuple[str, Dict]]:
    """Load PR2 sequences filtered by target taxon (case-insensitive)."""
    pr2_seqs = {}
    target_lower = target_taxon.lower()
    count = 0

    for record in SeqIO.parse(pr2_fasta, "fasta"):
        taxonomy = parse_pr2_header(record.description)
        lineage = taxonomy.get('full_lineage', '').lower()

        if target_lower in lineage:
            pr2_seqs[record.id] = (str(record.seq), taxonomy)
            count += 1
            if count >= limit:
                break

    return pr2_seqs


def main():
    parser = argparse.ArgumentParser(description="Validate primers against PR2 database")
    parser.add_argument("--taxon", help="Specific taxon folder to validate (e.g., 001_Euamoebida_order)")
    parser.add_argument("--all", action="store_true", help="Validate all taxa")
    parser.add_argument("--top-n", type=int, default=5, help="Number of top primers to test per taxon")
    parser.add_argument("--max-mismatches", type=int, default=2, help="Max mismatches allowed for primer binding")
    parser.add_argument("--limit-seqs", type=int, default=50000, help="Limit PR2 sequences to load (for speed)")
    parser.add_argument("--filter-taxon", action="store_true", help="Filter PR2 to target taxon only")
    args = parser.parse_args()

    if not args.taxon and not args.all:
        parser.error("Must specify --taxon or --all")

    print("=" * 70)
    print("🧬 PRIMER VALIDATION AGAINST PR2 DATABASE")
    print("=" * 70)

    # Get taxa to process
    if args.all:
        taxa_dirs = sorted([d for d in RESULTS_DIR.iterdir() if d.is_dir()])
    else:
        taxa_dirs = [RESULTS_DIR / args.taxon]

    # Load PR2 sequences (either filtered or full)
    pr2_seqs = None
    if not args.filter_taxon:
        print(f"\n📥 Loading PR2 database (limit: {args.limit_seqs} sequences)...")
        pr2_seqs = {}
        for i, record in enumerate(SeqIO.parse(PR2_FASTA, "fasta")):
            if i >= args.limit_seqs:
                break
            taxonomy = parse_pr2_header(record.description)
            pr2_seqs[record.id] = (str(record.seq), taxonomy)
        print(f"   Loaded {len(pr2_seqs):,} sequences")

    print(f"\n📋 Validating {len(taxa_dirs)} taxa, top {args.top_n} primers each\n")

    all_results = []
    for taxon_dir in taxa_dirs:
        if not taxon_dir.exists():
            print(f"⚠ Taxon not found: {taxon_dir.name}")
            continue

        # Extract taxon name from folder (e.g., "001_Euamoebida_order" -> "Euamoebida")
        taxon_name = taxon_dir.name.split('_')[1] if '_' in taxon_dir.name else taxon_dir.name

        # Load taxon-filtered PR2 if requested
        if args.filter_taxon:
            print(f"📥 Loading PR2 sequences for {taxon_name}...")
            pr2_seqs_local = load_pr2_by_taxon(PR2_FASTA, taxon_name, args.limit_seqs)
            print(f"   Found {len(pr2_seqs_local):,} {taxon_name} sequences")
        else:
            pr2_seqs_local = pr2_seqs

        print(f"🔍 {taxon_dir.name}...")
        primers = load_top_primers(taxon_dir, args.top_n)

        if not primers:
            print(f"   ⚠ No primers found")
            continue

        results = validate_primers(
            taxon_dir.name, primers, pr2_seqs_local,
            max_mismatches=args.max_mismatches
        )
        all_results.append(results)

        # Print summary for this taxon
        for pr in results['primer_results']:
            top_taxa = sorted(pr['taxonomy_breakdown'].items(), key=lambda x: -x[1])[:3]
            top_taxa_str = ", ".join([f"{t[0]}:{t[1]}" for t in top_taxa])
            print(f"   {pr['pair_id']}: {pr['total_hits']} hits | Top: {top_taxa_str[:60]}")

    # Summary
    print("\n" + "=" * 70)
    print("📊 SUMMARY")
    print("=" * 70)
    for res in all_results:
        total_hits = sum(pr['total_hits'] for pr in res['primer_results'])
        avg_hits = total_hits / len(res['primer_results']) if res['primer_results'] else 0
        print(f"   {res['taxon']}: avg {avg_hits:.0f} hits/primer")


if __name__ == "__main__":
    main()

