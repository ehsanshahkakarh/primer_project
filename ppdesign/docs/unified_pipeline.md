# Gene-Based Probe Design (`ppdesign unified`)

This pipeline targets conserved genes across genomes by chaining several well-known tools.

## Workflow Summary

1. **Genome acquisition** – Either accepts user FASTA/GBK files or queries the NeLLi database via taxonomic strings (`--db-query --taxonomy`).
2. **Gene calling** – Uses `prodigal-gv` to predict CDS features for each genome.
3. **Orthogrouping** – Runs ProteinOrtho to cluster proteins into orthologous groups across genomes.
4. **Alignments** – Aligns proteins with MAFFT, then converts them to codon-aware alignments (pal2nal).
5. **Probe selection** – Passes conserved codon blocks into the shared filtering engine (`probedesign_seqselection.py`) to apply GC/Tm/structure checks.

## Key CLI Flags

| Flag | Description |
| --- | --- |
| `--taxonomy` / `--db-query` | Query NeLLi when genomes are not supplied locally. |
| `--min-conservation` | Fraction of genomes that must contain the probe target. |
| `--threads` | Parallelism passed to ProteinOrtho and MAFFT. |
| `--max-probes` | Limit number of probes reported per orthogroup. |

## Outputs

- `results/<run>/summary.txt` – Run configuration and probe counts.
- `results/<run>/probes.csv` – Final probe panel with GC/Tm/structure metrics.
- `results/<run>/orthogroups.tsv` – ProteinOrtho cluster assignments.
- `results/<run>/alignment/` – Protein and codon alignments for reproducibility.

## Implementation Notes

- `probedesign_unified.py` now checks for required binaries (prodigal-gv, proteinortho, mafft) before any work begins and raises a descriptive error if one is absent.
- Alignment/parsing helpers live in `conserved_finder.py` and `thermodynamics.py`.
- Suggestion: wrap external command invocations with richer error handling (capture non-zero exit codes and surface stdout/stderr). Currently the pipeline raises generic `RuntimeError` on failure.
