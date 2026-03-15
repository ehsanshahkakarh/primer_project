#!/usr/bin/env python3
"""
Comprehensive PR2 validation for all primer candidates.
Outputs a summary CSV and detailed report.
"""

import csv
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime

# Paths
PROJECT_DIR = Path(__file__).parent.parent.absolute()
RESULTS_DIR = PROJECT_DIR / "primer_results"
PR2_FASTA = Path("/Users/ehsankakarh/PROJA/eukcensus_metadata/primer_amoebozoa/PR2/pr2_version_5.1.1_SSU_taxo_long.fasta")
OUTPUT_DIR = PROJECT_DIR / "validation_results"


def load_pr2_sequences(limit=None):
    """Load all PR2 sequences into memory for fast searching."""
    print("📥 Loading PR2 database...")
    sequences = []
    for i, record in enumerate(SeqIO.parse(PR2_FASTA, "fasta")):
        if limit and i >= limit:
            break
        # Parse taxonomy from header
        parts = record.description.split('|')
        taxonomy = {
            'domain': parts[4] if len(parts) > 4 else '',
            'supergroup': parts[5] if len(parts) > 5 else '',
            'phylum': parts[6] if len(parts) > 6 else '',
            'class': parts[7] if len(parts) > 7 else '',
            'order': parts[8] if len(parts) > 8 else '',
            'family': parts[9] if len(parts) > 9 else '',
            'genus': parts[10] if len(parts) > 10 else '',
        }
        sequences.append({
            'id': record.id,
            'seq': str(record.seq).upper(),
            'taxonomy': taxonomy
        })
    print(f"   Loaded {len(sequences):,} sequences")
    return sequences


def check_primer_hits(pr2_seqs, fwd_primer, rev_primer, anchor_len=10, min_amp=100, max_amp=1000):
    """Check how many PR2 sequences would be amplified by this primer pair."""
    rev_rc = str(Seq(rev_primer).reverse_complement())
    fwd_anchor = fwd_primer[:anchor_len].upper()
    rev_anchor = rev_rc[:anchor_len].upper()
    
    hits = []
    for seq_data in pr2_seqs:
        seq = seq_data['seq']
        fwd_pos = seq.find(fwd_anchor)
        rev_pos = seq.find(rev_anchor)
        
        if fwd_pos != -1 and rev_pos != -1:
            amp_size = rev_pos + len(rev_rc) - fwd_pos
            if min_amp <= amp_size <= max_amp:
                hits.append({
                    'id': seq_data['id'],
                    'amplicon_size': amp_size,
                    'taxonomy': seq_data['taxonomy']
                })
    return hits


def load_primers(taxon_dir, top_n=10):
    """Load top N primers from a taxon directory."""
    fasta_file = taxon_dir / "primer_pairs.fasta"
    if not fasta_file.exists():
        return []
    
    primers = []
    current = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        pair_id = record.id.rsplit('_', 1)[0]
        direction = record.id.rsplit('_', 1)[1]
        if pair_id not in current:
            current[pair_id] = {}
        current[pair_id][direction] = str(record.seq)
    
    for pair_id in sorted(current.keys())[:top_n]:
        if 'forward' in current[pair_id] and 'reverse' in current[pair_id]:
            primers.append((pair_id, current[pair_id]['forward'], current[pair_id]['reverse']))
    
    return primers


def get_target_taxon(folder_name):
    """Extract target taxon name from folder (e.g., '001_Euamoebida_order' -> 'Euamoebida')."""
    parts = folder_name.split('_')
    if len(parts) >= 2:
        return parts[1]
    return folder_name


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    print("=" * 80)
    print("🧬 COMPREHENSIVE PR2 PRIMER VALIDATION")
    print(f"   Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 80)
    
    # Load PR2 database
    pr2_seqs = load_pr2_sequences()
    
    # Get all taxa directories
    taxa_dirs = sorted([d for d in RESULTS_DIR.iterdir() if d.is_dir()])
    print(f"\n📋 Found {len(taxa_dirs)} taxa to validate\n")
    
    # Prepare output
    summary_rows = []
    detailed_report = []
    
    for taxon_dir in taxa_dirs:
        taxon_name = taxon_dir.name
        target_taxon = get_target_taxon(taxon_name)
        
        print(f"🔍 {taxon_name}...")
        primers = load_primers(taxon_dir, top_n=10)
        
        if not primers:
            print(f"   ⚠️ No primers found")
            continue
        
        for pair_id, fwd, rev in primers:
            hits = check_primer_hits(pr2_seqs, fwd, rev)
            
            # Count by supergroup
            supergroup_counts = defaultdict(int)
            target_hits = 0
            for h in hits:
                sg = h['taxonomy'].get('supergroup', 'Unknown')
                supergroup_counts[sg] += 1
                # Check if hit matches target taxon
                tax_str = '|'.join(h['taxonomy'].values()).lower()
                if target_taxon.lower() in tax_str:
                    target_hits += 1
            
            top_groups = sorted(supergroup_counts.items(), key=lambda x: -x[1])[:3]
            top_groups_str = "; ".join([f"{g}:{c}" for g, c in top_groups])
            
            # Specificity score
            specificity = (target_hits / len(hits) * 100) if hits else 0
            
            row = {
                'taxon': taxon_name,
                'target': target_taxon,
                'pair_id': pair_id,
                'fwd_primer': fwd,
                'rev_primer': rev,
                'total_hits': len(hits),
                'target_hits': target_hits,
                'specificity_pct': round(specificity, 1),
                'top_groups': top_groups_str,
                'amplicon_sizes': f"{min(h['amplicon_size'] for h in hits)}-{max(h['amplicon_size'] for h in hits)}" if hits else "N/A"
            }
            summary_rows.append(row)
            
            # Print progress
            status = "✓" if specificity > 50 else "⚠️" if hits else "❌"
            print(f"   {status} {pair_id}: {len(hits)} hits, {target_hits} target, {specificity:.0f}% spec | {top_groups_str[:50]}")
    
    # Write summary CSV
    csv_file = OUTPUT_DIR / "pr2_validation_summary.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)
    
    print(f"\n✅ Summary saved to: {csv_file}")
    print(f"   Total primer pairs validated: {len(summary_rows)}")


if __name__ == "__main__":
    main()

