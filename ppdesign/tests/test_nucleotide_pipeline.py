import shutil
from pathlib import Path

from ppdesign.probedesign_nucleotide import NucleotideProbeDesign


def test_kmer_pipeline_smoke(tmp_path):
    sequences = {
        "seq1": "ATGCGTACGTAGCTAGCTAG",
        "seq2": "ATGCGTACGTAGCTAGCTAA",
        "seq3": "ATGCGTACGTGGCTAGCTAG",
    }

    output_tag = "test_nuc_smoke"
    pipeline = NucleotideProbeDesign(output_dir=output_tag, threads=1)

    regions = pipeline.find_conserved_regions(
        sequences=sequences,
        method="kmer",
        min_length=6,
        max_length=8,
        min_conservation=0.66,
        kmer_size=4,
    )

    assert regions, "Expected at least one conserved region from k-mer finder"

    df = pipeline.filter_oligos(
        regions,
        gc_range=(0, 100),
        tm_range=(0, 120),
        check_hairpins=False,
        check_dimers=False,
    )

    assert not df.empty, "Filtered oligos should not be empty"

    results_dir = Path("results") / output_tag
    if results_dir.exists():
        shutil.rmtree(results_dir)
