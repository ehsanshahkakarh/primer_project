"""Pytest configuration and shared fixtures."""

import pytest
import tempfile
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    temp_path = tempfile.mkdtemp()
    yield Path(temp_path)
    shutil.rmtree(temp_path)


@pytest.fixture
def sample_fasta(temp_dir):
    """Create a sample FASTA file for testing."""
    sequences = [
        SeqRecord(Seq("ATCGATCGATCG"), id="seq1", description="Test sequence 1"),
        SeqRecord(Seq("GCTAGCTAGCTA"), id="seq2", description="Test sequence 2"),
        SeqRecord(Seq("ATCNATCNATCN"), id="seq3", description="Test sequence with Ns"),
    ]

    fasta_path = temp_dir / "test.fasta"
    SeqIO.write(sequences, fasta_path, "fasta")
    return fasta_path


@pytest.fixture
def sample_genomad_tsv(temp_dir):
    """Create a sample GeNomad taxonomy TSV file."""
    tsv_path = temp_dir / "test_virus_summary.tsv"

    content = """seq_name\ttaxonomy\tconfidence
contig1|provirus_1_1000\tViruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes\t0.95
contig2|provirus_1_2000\tViruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Faserviricetes\t0.90
contig3|provirus_1_3000\tViruses;Monodnaviria;Sangervirae;Phixviricota\t0.85
"""

    tsv_path.write_text(content)
    return tsv_path


@pytest.fixture
def sample_alignment(temp_dir):
    """Create a sample alignment file for testing."""
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import AlignIO

    # Create aligned sequences
    alignment = MultipleSeqAlignment(
        [
            SeqRecord(Seq("ATGATGATG"), id="seq1"),
            SeqRecord(Seq("ATGATGATG"), id="seq2"),
            SeqRecord(Seq("ATGCTGATG"), id="seq3"),  # One mismatch
        ]
    )

    alignment_path = temp_dir / "test_alignment.fasta"
    AlignIO.write(alignment, alignment_path, "fasta")
    return alignment_path
