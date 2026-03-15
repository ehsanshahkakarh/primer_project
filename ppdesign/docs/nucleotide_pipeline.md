# Direct Nucleotide Pipeline (`ppdesign nucleotide`)

Design probes straight from nucleotide contigs when annotated genes are unavailable.

## Available Methods

| Method | Engine | When to use |
| --- | --- | --- |
| `kmer` (default) | `KmerBasedFinder` (`conserved_finder.py`) | Fast scan for conserved windows in moderately sized contig sets. |
| `msa` | `MSABasedFinder` | Uses MUSCLE-based MSAs; appropriate for shorter amplicons or when structural conservation matters. |
| `minimap2` | `Minimap2BasedFinder` | Aligns long contigs/genomes using minimap2 and extracts conserved regions at ≥90% identity. |

## CLI Essentials

```bash
ppdesign nucleotide main \
  --fasta-input path/to/contigs \
  --method kmer \
  --min-length 18 --max-length 30 \
  --min-conservation 0.85 \
  --output-dir my_probes
```

Key parameters: `--kmer-size`, `--window-size`, `--gc-range`, `--tm-range`, `--threads`, `--allow-mismatches` (k-mer tolerance).

## Filtering & Scoring

After conserved blocks are found, the pipeline applies:

- GC, Tm, and nearest-neighbor percentage thresholds.
- Optional hairpin/dimer checks (`--no-hairpin`, `--no-dimer` to disable).
- A composite quality score defined in `NucleotideProbeDesign._calculate_quality_score`.

## Outputs

- `probes.csv` – Columns include Sequence, Length, GC%, Tm, Conservation, ND%, hairpin/dimer flags, quality score.
- `candidate_regions.csv` – Optional dump of all conserved segments before filtering.

## Implementation Observations

- K-mer method now honours the `allow_mismatches` argument when extending seeds; adjust the CLI flag to tune tolerance.
- The minimap2 branch assumes the `minimap2` binary is available; wrapping the call with a pre-flight check would improve UX.
- Consider caching intermediate alignments (especially for the MSA mode) to avoid recomputation when tweaking filters.
