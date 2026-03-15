# CRISPR gRNA Design Workflow

The gRNA pipeline (`ppdesign grna`) discovers SpCas9-compatible guides in one or more FASTA files, optionally enforcing coverage constraints across hundreds of contigs. The underlying engine lives in `src/ppdesign/guide_rna_finder.py` and was recently hardened in the following ways:

- **Overlap-aware PAM search** – NGG/NAG motifs are now scanned with a sliding window, ensuring overlapping PAMs (for example `GGGG`) each yield a candidate guide.
- **Canonical base enforcement** – Candidate guides are discarded unless their 20 nt payload consists solely of `A/C/G/T`. Ambiguous IUPAC bases are no longer permitted in emitted guides.
- **Perfect-coverage mode** – When `--perfect-coverage` is active, a greedy set-cover optimizer guarantees a minimum number of exact guides per sequence (`--min-coverage-per-target`) while respecting a global cap (`--max-total-grnas`).

## Typical invocation

```bash
ppdesign grna main \
  --fasta-input path/to/cluster.fna \
  --output-dir cluster_grna \
  --conservation 0.7 \
  --mismatches 2 \
  --min-guides 50 \
  --max-guides 400 \
  --min-gc 30 --max-gc 70 \
  --perfect-coverage \
  --min-coverage-per-target 3 \
  --max-total-grnas 400
```

### Key parameters

| Flag | Description |
| --- | --- |
| `--fasta-input` | Single FASTA file or directory of FASTA files to scan. |
| `--conservation` | Minimum fraction of sequences containing the guide (applied after clustering). Lower values allow more lineage-specific guides. |
| `--mismatches` | Mismatches tolerated while clustering guides prior to consensus building. |
| `--no-degenerate` | Emit only exact 20-mers; otherwise a consensus with IUPAC codes may be returned. |
| `--perfect-coverage` | Enable set-cover routine to guarantee coverage per target (requires `--min-coverage-per-target`). |
| `--min-gc` / `--max-gc` | GC bounds for candidate guides. |

### Output files

Each run produces a structured result folder (see `docs/results.md` for the full layout). Important artifacts include:

- `guide_rnas.csv` – quantitative summary for every surviving guide (coverage parameters reported when perfect coverage mode is active).
- `guide_rnas.fasta` – FASTA-formatted guides suitable for synthesis or downstream BLAST.
- `summary.txt` – human-readable metrics (counts, conservation, PAM breakdown, coverage statistics).
- `target_mapping.tsv` – mapping between guide IDs and contig coordinates (20 nt window plus strand) for the entire candidate set.
- Optional matrices (`coverage_matrix.tsv`, `target_coverage_details.tsv`) when perfect coverage is requested.

### Post-processing tips

- Use `target_mapping.tsv` to collapse redundant guides or to downselect to a non-overlapping panel.
- Combine the `coverage_matrix` with external metadata to prioritize guides that hit critical samples or strains.
- For large inputs, consider pre-dereplicating contigs or adjusting `--mismatches` to keep runtime manageable.
