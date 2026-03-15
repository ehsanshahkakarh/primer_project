# PPDesign Overview

PPDesign delivers a collection of pipelines for probe, primer, and CRISPR gRNA discovery. All command line tools share a common project layout (`src/ppdesign/`) and write their outputs beneath `results/<run_name>`. The table below highlights the primary front-ends and when to choose each one.

| CLI | Location | Primary use case | Key dependencies |
| --- | --- | --- | --- |
| `ppdesign unified` | `src/ppdesign/probedesign_unified.py` | Bacterial/archaeal probe design starting from genomes; performs gene calling, orthogrouping, protein & codon alignments | prodigal-gv, ProteinOrtho, MAFFT, pal2nal |
| `ppdesign nucleotide` | `src/ppdesign/probedesign_nucleotide.py` | Probe/primer design directly on nucleotide contigs or amplicons (k-mer, MSA, minimap2 workflows) | minimap2, MUSCLE |
| `ppdesign primer` | `src/ppdesign/probedesign_primer.py` | PCR primer pair design for amplifying conserved regions; generates forward/reverse pairs with amplicon size control (100bp-2kb) | Biopython, primer3-py (optional) |
| `ppdesign grna` | `src/ppdesign/probedesign_grna.py` | Conserved SpCas9 gRNA discovery on a FASTA collection; supports perfect-coverage mode | Biopython |
| `ppdesign select` | `src/ppdesign/probedesign_seqselection.py` | Filter & rank codon-aligned oligos by thermodynamics and structure | ViennaRNA (optional) |
| `ppdesign rank` | `src/ppdesign/probedesign_rank.py` | Post-hoc scoring/visualization of probe panels against a reference genome | matplotlib |

Common features:

- Pixi environment management (`pixi.toml`) with Python 3.11 and required bioinformatics packages.
- Uniform CLI flag style (`--fasta-input`, `--output-dir`, `--threads`).
- Consistent results layout with CSV + FASTA + summary text files for downstream review.
- Reusable Python API surface inside `src/ppdesign/` for notebook-driven workflows.

See `docs/grna_pipeline.md` for the CRISPR workflow, `docs/primer_pipeline.md` for PCR primer pair design, and `docs/results.md` for an explanation of the output artifacts that each run produces.
