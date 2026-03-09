# Multi-Size Primer Design

## Overview

The primer design pipeline has been updated to design primers for **three different product sizes** simultaneously:

- **500bp** (range: 450-550bp)
- **1000bp** (range: 950-1050bp)
- **1500bp** (range: 1450-1550bp)

This allows you to choose the most appropriate amplicon size for your downstream application.

## Why Multiple Sizes?

Different applications require different amplicon sizes:

| Product Size | Best For |
|--------------|----------|
| **500bp** | • Degraded DNA samples<br>• High-throughput sequencing<br>• qPCR applications<br>• Environmental samples |
| **1000bp** | • Standard Sanger sequencing<br>• Balanced between specificity and coverage<br>• Good for most applications |
| **1500bp** | • Maximum sequence coverage<br>• Phylogenetic analysis<br>• High-quality DNA samples<br>• Detailed taxonomic resolution |

## How It Works

### 1. Representative Selection & Alignment

The pipeline ensures taxonomic diversity and structural awareness:

1. **Select representatives**: One sequence per unique genus (largest cluster wins)
2. **Separate unknowns**: `.U.` sequences (unknown genus) are set aside
3. **Q-INS-i alignment**: Representatives are aligned using MAFFT's Q-INS-i algorithm, which considers **RNA secondary structure**
4. **Add unknowns**: `.U.` sequences are added to the reference MSA via `mafft --add`

**Why Q-INS-i?**
- rRNA (18S/16S) has significant secondary structure
- Q-INS-i identifies regions that are both **conserved** AND **structurally accessible**
- Primers in accessible regions hybridize better to the target

### 2. Single Consensus Sequence

The pipeline generates **one consensus sequence** from the alignment, which is used for all three product sizes.

### 2. Three Primer3 Runs

Primer3 is run **three times** with different `PRIMER_PRODUCT_SIZE_RANGE` settings:

```
500bp:  PRIMER_PRODUCT_SIZE_RANGE=450-550
1000bp: PRIMER_PRODUCT_SIZE_RANGE=950-1050
1500bp: PRIMER_PRODUCT_SIZE_RANGE=1450-1550
```

### 3. Same Primer Settings

All other primer3 parameters remain **identical** across all sizes:

- Optimal Tm: 60°C (range: 55-65°C)
- GC content: 40-60%
- Primer length: 18-25bp (optimal: 20bp)
- Number returned: 5 primer pairs per size

## Output Files

For each taxon, the following files are generated:

```
{taxon}/primers/
├── {taxon}_consensus.fasta           # Consensus sequence
├── {taxon}_primers_all_sizes.json    # Combined results for all sizes
├── {taxon}_primers_500bp.json        # Only 500bp primers
├── {taxon}_primers_1000bp.json       # Only 1000bp primers
└── {taxon}_primers_1500bp.json       # Only 1500bp primers
```

### File Structure

**`{taxon}_primers_all_sizes.json`** contains:
```json
{
  "500bp": {
    "target_size": 500,
    "size_range": [450, 550],
    "primers": [...]
  },
  "1000bp": {
    "target_size": 1000,
    "size_range": [950, 1050],
    "primers": [...]
  },
  "1500bp": {
    "target_size": 1500,
    "size_range": [1450, 1550],
    "primers": [...]
  },
  "stats": {...}
}
```

Each individual size file (e.g., `{taxon}_primers_500bp.json`) contains only the primers for that specific size.

## Usage

### Design primers for a single taxon

```bash
python scripts/design_primers.py \
  --alignment test_data/18S/Spizellomyces/align/Spizellomyces_aligned.fasta \
  --output test_data/18S/Spizellomyces/primers
```

### Process all alignments

```bash
python scripts/design_primers.py --gene 18S --all
```

## Example Output

When running the pipeline, you'll see output like:

```
🧬 Designing primers for: Spizellomyces
============================================================
   → Reading alignment...
     📊 Found 52 sequences in alignment
   → Generating consensus sequence...
     📊 Consensus length: 1789 bp
   → Analyzing conservation...
     📊 Average conservation: 94.2%
   → Finding conserved regions...
     📊 Found 12 conserved regions suitable for primers
   → Running primer3 for multiple product sizes...
     Target sizes: 500bp, 1000bp, 1500bp
     ✓ 500bp: Designed 5 primer pairs
        Best pair: 495bp product
        Forward: CTCATCCGAAACGCTGCATG (Tm=60.2°C, GC=50%)
        Reverse: TTGGCTGCGGATTACCCATT (Tm=59.8°C, GC=55%)
     ✓ 1000bp: Designed 5 primer pairs
        Best pair: 1005bp product
     ✓ 1500bp: Designed 5 primer pairs
        Best pair: 1498bp product

   📊 Summary: Successfully designed primers for 3/3 product sizes
```

## Handling Short Sequences

If the consensus sequence is shorter than a target product size, the pipeline will:

1. Skip that product size
2. Report an error in the output JSON
3. Continue with other product sizes

Example:
```json
{
  "500bp": {
    "primers": [...]
  },
  "1000bp": {
    "error": "Sequence too short (850 bp) for 1000bp product",
    "primers": []
  },
  "1500bp": {
    "error": "Sequence too short (850 bp) for 1500bp product",
    "primers": []
  }
}
```

## Testing

Run the demonstration script to see the configuration:

```bash
python scripts/test_multi_size_primers.py
```

This shows the settings and expected output structure without requiring primer3 to be installed.

