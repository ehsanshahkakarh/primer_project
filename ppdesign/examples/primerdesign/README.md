# Primer Design Example: Burkholderia 16S rRNA

Design PCR primers targeting conserved regions of the 16S rRNA gene across 10 Burkholderia species.

## Input

`burkholderia_16s.fna` — 10 sequences (1401-1550 bp) from NCBI:

| Accession   | Species              |
|-------------|----------------------|
| M22518.1    | B. cepacia           |
| X67038.1    | B. gladioli          |
| PX997116.1  | B. gladioli          |
| PX997090.1  | B. gladioli          |
| PX997089.1  | B. gladioli          |
| PX997088.1  | B. oklahomensis      |
| PX997084.1  | B. cenocepacia       |
| PX985651.1  | B. contaminans       |
| PX974699.1  | B. ambifaria         |
| PX969661.1  | B. gladioli          |

## Run

From the repository root:

```bash
# Short amplicons (200-400 bp, e.g. V3-V5 region)
PYTHONPATH=src pixi run -e dev python -m ppdesign.cli primer main \
  --fasta-input examples/primerdesign/burkholderia_16s.fna \
  --output-dir burkholderia_300bp \
  --amplicon-min 200 \
  --amplicon-max 400 \
  --conservation 0.8

# Near-full-length amplicons (1300-1600 bp)
PYTHONPATH=src pixi run -e dev python -m ppdesign.cli primer main \
  --fasta-input examples/primerdesign/burkholderia_16s.fna \
  --output-dir burkholderia_1500bp \
  --amplicon-min 1300 \
  --amplicon-max 1600 \
  --conservation 0.8

# Fast mode (k-mer, no MAFFT required)
PYTHONPATH=src pixi run -e dev python -m ppdesign.cli primer main \
  --fasta-input examples/primerdesign/burkholderia_16s.fna \
  --output-dir burkholderia_kmer \
  --method kmer \
  --amplicon-min 200 \
  --amplicon-max 400 \
  --conservation 0.8
```

The default method is `msa` (MAFFT L-INS-i alignment). Use `--method kmer` if MAFFT is not installed.

## Output Location

Relative paths are placed under `results/`:

```
results/burkholderia_300bp/
  primer_pairs.csv              # All primer pairs with full stats
  primer_pairs.fasta            # Sequences ready for ordering
  amplicon_predictions.tsv      # Genomic coordinates per pair
  summary.txt                   # Ranked overview of top 10 pairs
```

## Reading the Output

### summary.txt

Quick overview. Shows total pair count, amplicon/Tm/conservation statistics, and the top 10 pairs ranked by quality score (0-100):

```
1. Score: 99.3/100  Amplicon: 346bp  Conservation: 100%
   Fwd: TACGTAGGGTGCAAGCGTTAATCG  (24bp, Tm=60.8C, GC=50%, hairpin=-0.4 kcal/mol, pos=529)
   Rev: GGTCAACTTCACGCGTTAGATACGT  (25bp, Tm=60.7C, GC=48%, hairpin=-1.2 kcal/mol, pos=875)
   dTm=0.1C  cross-dimer=0
```

### primer_pairs.csv

One row per pair, sorted by quality score (best first). Key columns:

| Column | Meaning |
|--------|---------|
| `forward_seq` / `reverse_seq` | Primer sequences (5'->3') |
| `forward_length` / `reverse_length` | Primer length in bp |
| `forward_tm` / `reverse_tm` | Melting temperature (nearest-neighbor, C) |
| `forward_gc` / `reverse_gc` | GC content (%) |
| `forward_hairpin_dg` / `reverse_hairpin_dg` | Hairpin delta-G from primer3 (kcal/mol); values > -2 = no concern |
| `forward_pos` / `reverse_pos` | Binding position on reference alignment |
| `amplicon_size` | Distance between forward and reverse positions (bp) |
| `tm_difference` | abs(forward_tm - reverse_tm); lower is better |
| `cross_dimer_score` | Heterodimer potential; 0 = none detected |
| `quality_score` | Composite 0-100 score (conservation 35%, Tm match 25%, GC 15%, hairpin 15%, cross-dimer 10%) |
| `conservation` | Fraction of input sequences where both primers bind (1.0 = all) |

### primer_pairs.fasta

Each pair produces two FASTA entries (`pair_001_forward`, `pair_001_reverse`). Use directly for primer ordering.

### amplicon_predictions.tsv

Tab-separated. Shows amplicon start/end positions and which input sequences each pair targets:

```
pair_id  amplicon_start  amplicon_end  amplicon_size  target_sequences
pair_001 529             875           346            M22518.1,PX969661.1,...
```
