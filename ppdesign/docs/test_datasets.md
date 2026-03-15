# Test Data Sets

To validate the updated workflows we ship three lightweight FASTA collections (all under 100 KB):

| Dataset | Location | Purpose | Notes |
| --- | --- | --- | --- |
| `resources/test_data/grna/cluster45_subset.fna` | 5 Pepper mild mottle virus contigs (cluster45) | Exercises the gRNA finder, including PAM overlap detection and perfect coverage mode | Extracted from the cluster45 assembly; retains both Illumina contigs and the NCBI reference accession. |
| `resources/test_data/nucleotide/toy_contigs.fna` | 3 synthetic contigs | Quick smoke test for the nucleotide probe pipeline (k-mer mode) | Sequences differ at 1–2 positions to test mismatch handling. |
| `resources/test_data/unified/toy_genomes.fna` | 3 synthetic genomes with two ORFs each | Minimal input for the unified gene-based pipeline | Short genomes with clear start/stop codons suitable for prodigal-gv. |

These files are used in the new automated tests (`tests/test_nucleotide_pipeline.py`, `tests/test_unified_dependencies.py`) and for manual sanity checks described below.

## Manual Sanity Checks

1. **gRNA pipeline**
   ```bash
   ppdesign grna main \
     --fasta-input resources/test_data/grna/cluster45_subset.fna \
     --output-dir demo_cluster45 \
     --conservation 0.7 --mismatches 2 \
     --perfect-coverage --min-coverage-per-target 3 --max-total-grnas 50
   ```
   Confirms the revised clustering and coverage reporting. The run completes in under a minute.

2. **Nucleotide pipeline**
   ```bash
   ppdesign nucleotide main \
     --fasta-input resources/test_data/nucleotide/toy_contigs.fna \
     --method kmer --kmer-size 4 \
     --min-length 6 --max-length 8 \
     --output-dir demo_nucleotide
   ```
   Produces a small probe set and exercises the `--allow-mismatches` handling.

3. **Unified pipeline**
   ```bash
   ppdesign unified main \
     --fasta-input resources/test_data/unified/toy_genomes.fna \
     --output-dir demo_unified --threads 1
   ```
   This command requires external binaries (`prodigal-gv`, `proteinortho`, `mafft`). The pipeline now checks for these dependencies and exits with a clear error message if they are absent, making it straightforward to enable once the tools are installed.
