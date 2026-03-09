# Primer Design Pipeline

Automated pipeline for designing PCR primers targeting underrepresented microbial lineages from 18S/16S rRNA gene sequences.

## Overview

This pipeline extracts sequences for priority taxa, performs multiple sequence alignment, and designs primers using primer3.

## Requirements

- Python 3.8+
- MAFFT (for sequence alignment)
- primer3 (for primer design)
- Optional: seqkit (for length filtering)

## Installation

```bash
# Install MAFFT
brew install mafft  # macOS
# or: apt-get install mafft  # Ubuntu

# Install primer3
brew install primer3  # macOS
# or: apt-get install primer3  # Ubuntu
```

## Usage

### Run the full pipeline for all priority taxa:

```bash
python scripts/run_primer_pipeline.py --gene 18S
```

### Run for a specific taxon:

```bash
python scripts/run_primer_pipeline.py --gene 18S --taxon Hexamita_G
```

### Design primers from an existing alignment:

```bash
python scripts/design_primers.py --alignment path/to/aligned.fasta --output output_dir/
```

### Process all alignments:

```bash
python scripts/design_primers.py --gene 18S --all
```

## Pipeline Steps

1. **Extract sequences** - Extract sequences for the target taxon from census data
2. **Length filter** - Filter sequences ≥1200bp (requires seqkit)
3. **MAFFT alignment** - Multiple sequence alignment
4. **Consensus generation** - Generate consensus sequence from alignment
5. **Primer design** - Design primers using primer3

## Output Structure

```
{gene}/{taxon}/
├── {taxon}_{gene}.fasta          # Extracted sequences
├── {taxon}_sequences.tsv         # Sequence metadata
├── pipeline_stats.json           # Pipeline statistics
├── filtered/
│   └── {taxon}_1200bp.fasta      # Length-filtered sequences
├── align/
│   └── {taxon}_aligned.fasta           # MAFFT alignment
└── primers/
    ├── {taxon}_consensus.fasta         # Consensus sequence
    ├── {taxon}_primers_all_sizes.json  # All primers (500bp, 1000bp, 1500bp)
    ├── {taxon}_primers_500bp.json      # 500bp product primers only
    ├── {taxon}_primers_1000bp.json     # 1000bp product primers only
    └── {taxon}_primers_1500bp.json     # 1500bp product primers only
```

## Test Data

The `test_data/` directory contains example outputs for 6 taxa (3 family-genus pairs):

| Family | Genus | Sequences |
|--------|-------|-----------|
| Spizellomycetaceae | Spizellomyces | 52 |
| Vannellidae | Vannella | 56 |
| Paramoebidae | Paramoeba | 20-54 |

## Statistics Tracked

The pipeline tracks detailed statistics at each step:

- **Extraction**: sequence count, length range, average length
- **Filtering**: filtered count, retention percentage
- **Alignment**: input/output counts, aligned length, time
- **Primers**: conservation analysis, conserved regions, primer pairs

## Scripts

- `run_primer_pipeline.py` - Main pipeline orchestrator
- `design_primers.py` - Primer design from alignments
- `extract_sequences_by_taxon.py` - Sequence extraction
- `build_lineage_index.py` - Build taxonomic lineage index
- `parse_priority_taxa.py` - Parse priority taxa list

## Primer Design Parameters

The pipeline now designs primers for **three different product sizes**:
- **500bp** (range: 450-550bp)
- **1000bp** (range: 950-1050bp)
- **1500bp** (range: 1450-1550bp)

### Primer3 Settings (Applied to All Sizes)

- Optimal Tm: 60°C (range: 55-65°C)
- Optimal GC: 50% (range: 40-60%)
- Primer length: 18-25bp (optimal: 20bp)
- Number of primer pairs returned: 5 per product size

