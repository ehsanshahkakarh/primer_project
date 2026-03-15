#!/usr/bin/env python3
"""
Primer Design Pipeline Orchestrator

This script runs the complete primer design pipeline for all dark_blue priority taxa
from the iTOL tree annotations.

Pipeline Steps:
1. Extract sequence IDs from cluster TSV by taxonomic group
2. Extract FASTA sequences from main census file
3. Filter sequences by length (≥1200bp)
4. Select representatives (1 per genus)
5. MAFFT alignment
6. Consensus generation

Usage:
    python run_primer_pipeline.py --gene 18S [--dry-run] [--taxon <specific_taxon>]
    python run_primer_pipeline.py --gene 16S [--dry-run]
"""

import subprocess
import argparse
import json
import sys
from pathlib import Path
from datetime import datetime

# Import the lineage index builder
sys.path.insert(0, str(Path(__file__).parent))
from build_lineage_index import build_index_from_csv, build_ancestor_index, CENSUS_CSV_PATHS

# Configuration paths (relative to project root)
CONFIG = {
    '18S': {
        'cluster_tsv': 'primer_design_pipeline/metadata/eukcensus_2025_18S_clusters.tsv',
        'main_fasta': 'primer_design_pipeline/metadata/eukcensus_2025_18S.fna',
        'priority_list': 'primer_design_pipeline/metadata/18s_priority_prune_list.txt',
        'output_base': 'primer_design_pipeline/18S',
        'min_length': 1200,
    },
    '16S': {
        'cluster_tsv': 'primer_design_pipeline/metadata/eukcensus_2025_16S_clusters.tsv',
        'main_fasta': 'primer_design_pipeline/metadata/eukcensus_2025_16S.fna',
        'priority_list': 'primer_design_pipeline/metadata/16s_priority_prune_list.txt',
        'output_base': 'primer_design_pipeline/16S',
        'min_length': 1200,
    }
}

# Map rank suffix to query column
RANK_MAP = {
    '_P': 'division',
    '_C': 'class',  # Needs special handling
    '_O': 'order',  # Needs special handling
    '_F': 'family',
    '_G': 'genus',
    '_S': 'species',
    '_I': 'infraclass',
}

def parse_prune_list(prune_list_file: str) -> list:
    """Parse iTOL prune list file and extract taxa."""
    taxa = []
    with open(prune_list_file, 'r') as f:
        in_data = False
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line == 'DATA':
                in_data = True
                continue
            if not in_data:
                continue

            # Parse taxon name and rank
            for suffix, rank in RANK_MAP.items():
                if line.endswith(suffix):
                    base_name = line[:-len(suffix)]
                    taxa.append({
                        'original': line,
                        'base_name': base_name,
                        'rank': rank,
                        'suffix': suffix
                    })
                    break
    return taxa


def ensure_lineage_index(gene: str, project_root: Path) -> Path:
    """Build lineage index if it doesn't exist, return path to index file."""
    output_dir = project_root / 'primer_design_pipeline' / gene
    output_dir.mkdir(parents=True, exist_ok=True)
    index_file = output_dir / f"lineage_index_{gene}.json"

    if index_file.exists():
        return index_file

    print(f"\n🔗 Building lineage index for {gene}...")

    # Parse all census CSV files
    all_taxa = {}
    csv_paths = CENSUS_CSV_PATHS.get(gene, {})

    for rank, rel_path in csv_paths.items():
        csv_file = project_root / rel_path
        if csv_file.exists():
            print(f"   📖 Parsing {rank} CSV: {csv_file.name}")
            taxa = build_index_from_csv(str(csv_file), rank)
            all_taxa.update(taxa)

    # Build ancestor index
    ancestor_index = build_ancestor_index(all_taxa)

    # Save to JSON
    with open(index_file, 'w') as f:
        json.dump({
            'gene': gene,
            'total_taxa': len(all_taxa),
            'total_ancestors': len(ancestor_index),
            'ancestor_index': ancestor_index
        }, f, indent=2)

    print(f"   ✓ Indexed {len(ancestor_index)} ancestor names")
    return index_file


def run_step(cmd: list, step_name: str, dry_run: bool = False) -> bool:
    """Run a pipeline step."""
    print(f"   → {step_name}")
    if dry_run:
        print(f"     [DRY-RUN] {' '.join(cmd)}")
        return True

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            print(f"     ❌ Error: {result.stderr}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"     ⏰ Timeout after 1 hour")
        return False
    except Exception as e:
        print(f"     ❌ Exception: {e}")
        return False

def count_fasta_sequences(fasta_file: Path) -> int:
    """Count number of sequences in a FASTA file."""
    if not fasta_file.exists():
        return 0
    count = 0
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def get_fasta_stats(fasta_file: Path) -> dict:
    """Get detailed stats from a FASTA file."""
    if not fasta_file.exists():
        return {'count': 0, 'total_bp': 0, 'min_len': 0, 'max_len': 0, 'avg_len': 0}

    lengths = []
    current_len = 0
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_len > 0:
            lengths.append(current_len)

    if not lengths:
        return {'count': 0, 'total_bp': 0, 'min_len': 0, 'max_len': 0, 'avg_len': 0}

    return {
        'count': len(lengths),
        'total_bp': sum(lengths),
        'min_len': min(lengths),
        'max_len': max(lengths),
        'avg_len': round(sum(lengths) / len(lengths), 1)
    }


def select_representatives(sequences_tsv: Path, fasta_input: Path,
                           reps_output: Path, unknown_output: Path) -> tuple:
    """
    Select one representative per unique genus from extracted sequences.

    For each unique genus:
    - Pick the centroid with the largest cluster size
    - Only considers sequences that are actually in the input FASTA
    - Sequences with '.U.' in genus go to unknown_output (for --add later)

    Args:
        sequences_tsv: TSV file with sequence metadata (has genus, cluster_size columns)
        fasta_input: Input FASTA file with all sequences
        reps_output: Output FASTA for representative sequences
        unknown_output: Output FASTA for unknown (.U.) sequences

    Returns:
        (num_representatives, num_unknown)
    """
    import csv

    # First, read FASTA to get available sequence IDs
    fasta_seqs = {}
    current_id = None
    current_seq = []

    with open(fasta_input, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    fasta_seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            fasta_seqs[current_id] = ''.join(current_seq)

    available_ids = set(fasta_seqs.keys())

    # Read sequence metadata from TSV - only consider sequences in FASTA
    genus_best = {}  # genus -> (seq_id, cluster_size)
    unknown_ids = set()

    with open(sequences_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_id = row.get('sequence_id', row.get('centroid', ''))

            # Skip if this sequence is not in the FASTA
            if seq_id not in available_ids:
                continue

            genus = row.get('genus', '')

            # Handle cluster_size - might be string
            try:
                cluster_size = int(row.get('cluster_size', row.get('size', 1)))
            except (ValueError, TypeError):
                cluster_size = 1

            # Check if this is an unknown genus
            if '.U.' in genus or genus.endswith('.U'):
                unknown_ids.add(seq_id)
            else:
                # Track best representative per genus (largest cluster)
                if genus not in genus_best or cluster_size > genus_best[genus][1]:
                    genus_best[genus] = (seq_id, cluster_size)

    # Get representative IDs
    rep_ids = {info[0] for info in genus_best.values()}

    # Write representatives
    reps_written = 0
    with open(reps_output, 'w') as f:
        for seq_id in rep_ids:
            f.write(f">{seq_id}\n{fasta_seqs[seq_id]}\n")
            reps_written += 1

    # Write unknowns
    unknowns_written = 0
    with open(unknown_output, 'w') as f:
        for seq_id in unknown_ids:
            f.write(f">{seq_id}\n{fasta_seqs[seq_id]}\n")
            unknowns_written += 1

    return reps_written, unknowns_written


def run_pipeline_for_taxon(taxon: dict, config: dict, project_root: Path, dry_run: bool = False) -> dict:
    """Run the complete pipeline for a single taxon."""
    taxon_name = taxon['base_name']
    taxon_orig = taxon['original']
    rank = taxon['rank']

    # Create output directory
    output_dir = project_root / config['output_base'] / taxon_name
    output_dir.mkdir(parents=True, exist_ok=True)

    result = {
        'taxon': taxon_orig,
        'status': 'pending',
        'steps_completed': [],
        'errors': [],
        'stats': {
            'extraction': {},
            'filtering': {},
            'alignment': {},
            'primers': {}
        }
    }

    print(f"\n{'='*60}")
    print(f"🧬 Processing: {taxon_orig}")
    print(f"   Rank: {rank}, Base name: {taxon_name}")
    print(f"   Output: {output_dir}")
    print(f"{'='*60}")

    # Check if already completed (alignment file exists and is non-empty)
    gene = config['output_base'].split('/')[-1]  # Extract gene from path (18S or 16S)
    align_dir = output_dir / 'align'
    aligned_fasta = align_dir / f"{taxon_name}_aligned.fasta"
    if not dry_run and aligned_fasta.exists() and aligned_fasta.stat().st_size > 0:
        print(f"   ✓ Already completed - skipping")
        result['status'] = 'complete'
        result['steps_completed'] = ['extract', 'filter', 'align']
        return result

    # STEP 1 & 2: Extract sequences
    scripts_dir = project_root / 'primer_design_pipeline' / 'scripts'
    extract_script = scripts_dir / 'extract_sequences_by_taxon.py'

    cmd = [
        'python3', str(extract_script),
        '--tsv', str(project_root / config['cluster_tsv']),
        '--taxon', taxon_name,
        '--rank', rank if rank != 'species' else 'species',
        '--output-dir', str(output_dir),
        '--fasta', str(project_root / config['main_fasta']),
        '--taxon-original', taxon_orig,
        '--gene', gene
    ]

    # Add lineage index for class/order/infraclass level taxa (auto-built if needed)
    if rank in ['class', 'order', 'infraclass']:
        lineage_index_file = project_root / config['output_base'] / f"lineage_index_{gene}.json"
        cmd.extend(['--lineage-index', str(lineage_index_file)])

    if not run_step(cmd, "Step 1-2: Extract sequences", dry_run):
        result['status'] = 'failed'
        result['errors'].append('Sequence extraction failed')
        return result
    result['steps_completed'].append('extract')

    # Check if FASTA was created and get stats
    fasta_file = output_dir / f"{taxon_name}_{gene}.fasta"
    if not dry_run and not fasta_file.exists():
        print(f"   ⚠️ No FASTA file created - taxon may not exist in cluster TSV")
        result['status'] = 'no_sequences'
        return result

    # Get extraction stats (from extraction_stats.json if available, plus FASTA stats)
    if not dry_run:
        extract_stats = get_fasta_stats(fasta_file)

        # Try to load detailed extraction stats from JSON
        extraction_json = output_dir / 'extraction_stats.json'
        if extraction_json.exists():
            with open(extraction_json) as f:
                detailed_stats = json.load(f)
            extract_stats['centroids'] = detailed_stats.get('centroids', extract_stats['count'])
            extract_stats['cluster_members'] = detailed_stats.get('cluster_members', 0)
            extract_stats['tsv_total'] = detailed_stats.get('tsv_total', extract_stats['count'])
        else:
            # Fallback: all sequences are centroids
            extract_stats['centroids'] = extract_stats['count']
            extract_stats['cluster_members'] = 0
            extract_stats['tsv_total'] = extract_stats['count']

        result['stats']['extraction'] = extract_stats
        print(f"     📊 Extracted: {extract_stats['count']} sequences (centroids)")
        if extract_stats['cluster_members'] > 0:
            print(f"        Cluster members in TSV: {extract_stats['cluster_members']} (represented by centroids at 97% OTU)")
        print(f"        Length range: {extract_stats['min_len']}-{extract_stats['max_len']} bp (avg: {extract_stats['avg_len']})")

    # STEP 3: Length filter
    filtered_dir = output_dir / 'filtered'
    filtered_dir.mkdir(exist_ok=True)
    filtered_fasta = filtered_dir / f"{taxon_name}_1200bp.fasta"

    cmd = ['seqkit', 'seq', '-m', str(config['min_length']), str(fasta_file)]
    if not dry_run:
        try:
            with open(filtered_fasta, 'w') as f:
                subprocess.run(cmd, stdout=f, check=True, timeout=300)
            result['steps_completed'].append('filter')

            # Get filter stats
            filter_stats = get_fasta_stats(filtered_fasta)
            result['stats']['filtering'] = filter_stats

            # Calculate retention rate
            if extract_stats['count'] > 0:
                retention = round(100 * filter_stats['count'] / extract_stats['count'], 1)
            else:
                retention = 0
            result['stats']['filtering']['retention_pct'] = retention

            print(f"   → Step 3: Length filter (≥{config['min_length']}bp)")
            print(f"     📊 After filter: {filter_stats['count']} sequences ({retention}% retained)")
            if filter_stats['count'] > 0:
                print(f"        Length range: {filter_stats['min_len']}-{filter_stats['max_len']} bp (avg: {filter_stats['avg_len']})")
        except Exception as e:
            result['errors'].append(f'Filter failed: {e}')
            print(f"   → Step 3: Length filter - FAILED (seqkit not available)")
    else:
        print(f"   → Step 3: Length filter [DRY-RUN]")
        result['steps_completed'].append('filter')

    # Save stats to file
    if not dry_run:
        stats_file = filtered_dir / 'filter_stats.txt'
        try:
            stats_result = subprocess.run(
                ['seqkit', 'stats', str(fasta_file), str(filtered_fasta)],
                capture_output=True, text=True
            )
            with open(stats_file, 'w') as f:
                f.write(stats_result.stdout)
        except Exception:
            pass

    # STEP 4: Select representatives (1 per unique genus)
    reps_dir = output_dir / 'representatives'
    reps_dir.mkdir(exist_ok=True)
    reps_fasta = reps_dir / f"{taxon_name}_representatives.fasta"
    unknown_fasta = reps_dir / f"{taxon_name}_unknown.fasta"

    # Determine input file - prefer filtered, fall back to original
    if dry_run:
        rep_input = filtered_fasta
    elif filtered_fasta.exists() and filtered_fasta.stat().st_size > 0:
        rep_input = filtered_fasta
    else:
        rep_input = fasta_file

    print(f"   → Step 4: Select representatives (1 per genus) {'[DRY-RUN]' if dry_run else ''}")
    if not dry_run:
        # Read the sequences TSV to get genus info
        sequences_tsv = output_dir / f"{taxon_name}_sequences.tsv"
        if sequences_tsv.exists():
            reps_count, unknown_count = select_representatives(
                sequences_tsv, rep_input, reps_fasta, unknown_fasta
            )
            result['steps_completed'].append('select_reps')
            result['stats']['representatives'] = {
                'unique_genera': reps_count,
                'unknown_sequences': unknown_count
            }
            print(f"     📊 Representatives: {reps_count} unique genera")
            if unknown_count > 0:
                print(f"     📊 Unknown (.U.): {unknown_count} sequences (will add via --add)")
        else:
            print(f"     ⚠️ No sequences TSV found, using all sequences")
            reps_fasta = rep_input
            unknown_count = 0
    else:
        result['steps_completed'].append('select_reps')

    # STEP 5: MAFFT alignment using Q-INS-i (considers RNA secondary structure)
    # Q-INS-i is slower but critical for rRNA to identify structurally accessible regions
    align_dir = output_dir / 'align'
    align_dir.mkdir(exist_ok=True)
    aligned_fasta = align_dir / f"{taxon_name}_aligned.fasta"

    print(f"   → Step 5: MAFFT Q-INS-i alignment (secondary structure aware) {'[DRY-RUN]' if dry_run else ''}")
    if not dry_run:
        try:
            import time
            start_time = time.time()

            # Check if we have representatives
            if reps_fasta.exists() and reps_fasta.stat().st_size > 0:
                reps_stats = get_fasta_stats(reps_fasta)
                result['stats']['alignment'] = {'input_representatives': reps_stats['count']}
                result['stats']['alignment']['algorithm'] = 'Q-INS-i'
                print(f"     📊 Aligning {reps_stats['count']} representative sequences")

                # Step 5a: Align representatives using Q-INS-i
                # Q-INS-i considers secondary structure via four-way consistency
                ref_aligned = align_dir / f"{taxon_name}_ref_aligned.fasta"
                cmd = ['mafft-qinsi', '--thread', '4', str(reps_fasta)]
                with open(ref_aligned, 'w') as f:
                    subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True, timeout=14400)

                # Step 5b: Add unknown sequences if any
                if unknown_fasta.exists() and unknown_fasta.stat().st_size > 0:
                    unknown_stats = get_fasta_stats(unknown_fasta)
                    print(f"     📊 Adding {unknown_stats['count']} unknown sequences via --add")
                    result['stats']['alignment']['unknown_added'] = unknown_stats['count']

                    # Use --add with the Q-INS-i aligned reference
                    cmd_add = ['mafft', '--add', str(unknown_fasta), '--thread', '4', str(ref_aligned)]
                    with open(aligned_fasta, 'w') as f:
                        subprocess.run(cmd_add, stdout=f, stderr=subprocess.DEVNULL, check=True, timeout=7200)
                else:
                    # No unknowns, just use the reference alignment
                    import shutil
                    shutil.copy(ref_aligned, aligned_fasta)
            else:
                # Fallback: align all sequences directly with Q-INS-i
                mafft_input = rep_input
                input_stats = get_fasta_stats(mafft_input)
                result['stats']['alignment'] = {'input_sequences': input_stats['count'], 'algorithm': 'Q-INS-i'}
                print(f"     📊 Input for alignment: {input_stats['count']} sequences")

                cmd = ['mafft-qinsi', '--thread', '4', str(mafft_input)]
                with open(aligned_fasta, 'w') as f:
                    subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True, timeout=14400)

            elapsed = round(time.time() - start_time, 1)
            result['steps_completed'].append('align')

            # Get alignment stats
            align_stats = get_fasta_stats(aligned_fasta)
            result['stats']['alignment']['output_sequences'] = align_stats['count']
            result['stats']['alignment']['aligned_length'] = align_stats['max_len']
            result['stats']['alignment']['time_seconds'] = elapsed

            print(f"     ✓ Alignment complete in {elapsed}s")
            print(f"     📊 Aligned: {align_stats['count']} sequences, length: {align_stats['max_len']} bp")
        except subprocess.TimeoutExpired:
            result['errors'].append('MAFFT timeout (2h)')
        except Exception as e:
            result['errors'].append(f'MAFFT failed: {e}')
    else:
        result['steps_completed'].append('align')

    # =========================================================================
    # STEP 6: Design primers using consensus and Primer3
    # =========================================================================
    print(f"   → Step 6: Primer design (consensus + primer3) {'[DRY-RUN]' if dry_run else ''}")
    primers_dir = output_dir / 'primers'
    primers_dir.mkdir(exist_ok=True)

    if not dry_run and aligned_fasta.exists() and aligned_fasta.stat().st_size > 0:
        try:
            # Import design_primers module
            import sys
            scripts_dir = Path(__file__).parent
            if str(scripts_dir) not in sys.path:
                sys.path.insert(0, str(scripts_dir))
            from design_primers import design_primers_for_taxon

            primer_result = design_primers_for_taxon(
                alignment_file=aligned_fasta,
                output_dir=primers_dir,
                taxon_name=taxon_name,
                max_seqs=500
            )

            result['steps_completed'].append('primers')
            result['stats']['primers'] = {
                'consensus_length': primer_result.get('consensus_length', 0),
                'conserved_regions': primer_result.get('conserved_regions', 0),
                'primers_designed': {}
            }

            # Count primers by size
            if 'primers_by_size' in primer_result:
                for size_label, size_data in primer_result['primers_by_size'].items():
                    primers = size_data.get('primers', [])
                    result['stats']['primers']['primers_designed'][size_label] = len(primers)

                total_primers = sum(result['stats']['primers']['primers_designed'].values())
                print(f"     ✓ Designed {total_primers} primer pairs")
                for size_label, count in result['stats']['primers']['primers_designed'].items():
                    print(f"       • {size_label}: {count} pairs")
            else:
                print(f"     ⚠️ No primers designed (check alignment quality)")

        except ImportError as e:
            result['errors'].append(f'Primer design import failed: {e}')
            print(f"     ❌ Could not import design_primers module: {e}")
        except Exception as e:
            result['errors'].append(f'Primer design failed: {e}')
            print(f"     ❌ Primer design error: {e}")
    else:
        if dry_run:
            result['steps_completed'].append('primers')
        else:
            result['errors'].append('No alignment file for primer design')

    result['status'] = 'complete' if len(result['errors']) == 0 else 'partial'

    # Save detailed stats to JSON file
    if not dry_run:
        stats_json = output_dir / 'pipeline_stats.json'
        with open(stats_json, 'w') as f:
            json.dump(result['stats'], f, indent=2)

    return result


def main():
    parser = argparse.ArgumentParser(description='Run primer design pipeline for priority taxa')
    parser.add_argument('--gene', required=True, choices=['18S', '16S'], help='Gene target (18S or 16S)')
    parser.add_argument('--dry-run', action='store_true', help='Print commands without executing')
    parser.add_argument('--taxon', help='Process only a specific taxon (by original name with suffix)')
    parser.add_argument('--limit', type=int, help='Limit number of taxa to process')

    args = parser.parse_args()

    # Determine project root (2 levels up from scripts dir)
    project_root = Path(__file__).resolve().parent.parent.parent
    config = CONFIG[args.gene]

    print(f"🧬 Primer Design Pipeline - {args.gene}")
    print(f"{'='*60}")
    print(f"Project root: {project_root}")
    print(f"Priority list: {config['priority_list']}")
    print(f"Cluster TSV: {config['cluster_tsv']}")
    print(f"Main FASTA: {config['main_fasta']}")
    print(f"{'='*60}")

    # Parse priority taxa
    priority_list = project_root / config['priority_list']
    taxa = parse_prune_list(str(priority_list))
    print(f"\n📋 Found {len(taxa)} priority taxa")

    # Ensure lineage index exists (auto-build if needed)
    ensure_lineage_index(args.gene, project_root)

    # Filter to specific taxon if requested
    if args.taxon:
        taxa = [t for t in taxa if t['original'] == args.taxon]
        if not taxa:
            print(f"❌ Taxon not found: {args.taxon}")
            sys.exit(1)

    # Apply limit
    if args.limit:
        taxa = taxa[:args.limit]
        print(f"   Limiting to first {args.limit} taxa")

    # Track results
    results = []
    start_time = datetime.now()

    for i, taxon in enumerate(taxa, 1):
        print(f"\n[{i}/{len(taxa)}] Processing {taxon['original']}...")
        result = run_pipeline_for_taxon(taxon, config, project_root, args.dry_run)
        results.append(result)

    # Summary
    elapsed = datetime.now() - start_time
    print(f"\n{'='*60}")
    print(f"🏁 PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"Elapsed time: {elapsed}")

    # Count by status
    by_status = {}
    for r in results:
        by_status[r['status']] = by_status.get(r['status'], 0) + 1

    print(f"\n📊 Results by status:")
    for status, count in sorted(by_status.items()):
        print(f"   {status}: {count}")

    # Save results
    results_file = project_root / config['output_base'] / 'pipeline_results.json'
    with open(results_file, 'w') as f:
        json.dump({
            'gene': args.gene,
            'timestamp': datetime.now().isoformat(),
            'elapsed_seconds': elapsed.total_seconds(),
            'total_taxa': len(results),
            'results': results
        }, f, indent=2)
    print(f"\n💾 Results saved: {results_file}")

if __name__ == "__main__":
    main()

