# Result Files Cheat Sheet

Every PPDesign CLI writes its outputs to `results/<run_name>/`. While individual pipelines add mode-specific tables, the following patterns are shared across tools.

## Core artifacts

| File | Present in | Purpose |
| --- | --- | --- |
| `summary.txt` | all CLIs | Human-readable overview (inputs, parameter values, counts, highlight statistics). |
| `*.csv` | all CLIs | Machine-friendly tabular output (e.g., `guide_rnas.csv`, `probes.csv`). |
| `*.fasta` | all CLIs | Final oligo sequences suitable for synthesis or downstream BLAST. |
| `target_mapping.tsv` | gRNA + probe modules | Maps each candidate to contig coordinates and strands. |
| `coverage_matrix.tsv` | gRNA (perfect coverage) | Binary matrix showing which guides cover which sequences. |
| `target_coverage_details.tsv` | gRNA (perfect coverage) | Per-contig coverage counts and the guides that contribute to coverage. |

## Mode-specific notes

### Gene-based (`ppdesign unified`)

- `orthogroups.tsv` – orthogroup assignments from ProteinOrtho.
- `alignment/` – intermediate protein and codon alignments.
- `probes_with_metrics.csv` – per-probe GC/Tm/structure metrics.

### Nucleotide (`ppdesign nucleotide`)

- `kmer_hits.tsv` – seed positions used to anchor candidate design blocks.
- `msa/` – multi-sequence alignment snapshots when the `--method msa` path is selected.

### gRNA (`ppdesign grna`)

- `guide_rnas.csv` – includes conservation, degeneracy count, quality score, and target list.
- `guide_rnas.fasta` – sequence payload; headers match the Guide_ID column.
- `coverage_matrix.tsv` / `target_coverage_details.tsv` – only when `--perfect-coverage` is active.
- `*.json` – optional machine-readable summaries created during post-processing.

## Naming conventions

- Guide/Probe IDs follow zero-padded numeric ordering (`gRNA_0001`, `probe_0001`, etc.).
- Contig identifiers are not altered; multi-part headers retain the original pipe-delimited tokens from the FASTA.
- All counts in summaries are exact hits (no mismatch tolerance) unless explicitly labelled otherwise.

Refer to the project README for instructions on running individual pipelines and to `docs/grna_pipeline.md` for gRNA-specific command guidance.
