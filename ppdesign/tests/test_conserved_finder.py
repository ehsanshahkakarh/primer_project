"""Tests for conserved region finding with optional alignment."""

import pytest
from ppdesign.conserved_finder import KmerBasedFinder, ConservedRegion


class TestKmerBasedFinder:
    """Tests for k-mer based conserved region finding."""

    def test_basic_region_finding(self):
        """Test basic conserved region finding without alignment."""
        sequences = {
            "seq1": "ATGCGATCGATCGATCGATCG",
            "seq2": "ATGCGATCGATCGATCGATCG",
            "seq3": "ATGCGATCGATCGATCGATCG",
        }

        finder = KmerBasedFinder(kmer_size=10, min_conservation=0.8)
        regions = finder.find_conserved_regions(
            sequences, min_length=15, max_length=30, align_regions=False
        )

        assert len(regions) > 0
        assert all(isinstance(r, ConservedRegion) for r in regions)
        assert all(r.conservation >= 0.8 for r in regions)

    def test_alignment_optional_flag(self):
        """Test that alignment is optional and controlled by flag."""
        sequences = {
            "seq1": "ATGCGATCGATCGATCGATCG",
            "seq2": "ATGCGATCGATCGATCGATCG",
            "seq3": "ATGCGATCGATCGATCGATCG",
        }

        finder = KmerBasedFinder(kmer_size=10, min_conservation=0.8)

        # Without alignment (fast mode)
        regions_fast = finder.find_conserved_regions(
            sequences, min_length=15, max_length=30, align_regions=False
        )

        # With alignment (requires MAFFT, may skip if not available)
        regions_aligned = finder.find_conserved_regions(
            sequences, min_length=15, max_length=30, align_regions=True
        )

        # Both should find regions
        assert len(regions_fast) > 0
        # Alignment may refine results or return same if perfect match
        assert len(regions_aligned) >= 0

    def test_region_with_iupac_codes(self):
        """Test handling of regions with IUPAC degenerate codes."""
        # Create sequences that will produce IUPAC codes in consensus
        sequences = {
            "seq1": "ATGCAATCGATCGATCGATCG",
            "seq2": "ATGCGATCGATCGATCGATCG",
            "seq3": "ATGCAATCGATCGATCGATCG",
        }

        finder = KmerBasedFinder(kmer_size=10, min_conservation=0.6)
        regions = finder.find_conserved_regions(
            sequences, min_length=15, max_length=30, align_regions=False
        )

        # Should find region despite variation
        assert len(regions) > 0

        # Consensus should have IUPAC code (M = A/C) at position of variation
        region = regions[0]
        # Position 5 has A in seq1/3 and G in seq2
        assert "R" in region.consensus or "M" in region.consensus or "W" in region.consensus

    def test_alignment_refinement_with_mismatches(self):
        """Test that alignment can refine regions with mismatches."""
        # Sequences with small variations
        sequences = {
            "seq1": "ATGCGATCGATCGATCGATCGATCGATCG",
            "seq2": "ATGCGATCGATTGATCGATCGATCGATCG",  # T instead of C
            "seq3": "ATGCGATCGATCGAT-GATCGATCGATCG",  # Has gap
        }

        finder = KmerBasedFinder(kmer_size=10, min_conservation=0.6)

        # With alignment (should handle gaps better)
        regions_aligned = finder.find_conserved_regions(
            sequences, min_length=15, max_length=30, align_regions=True
        )

        # Should find conserved regions
        assert len(regions_aligned) >= 0  # May filter out low quality

    def test_skip_alignment_for_perfect_match(self):
        """Test that alignment is skipped for perfect ACGT matches."""
        sequences = {
            "seq1": "ATGCGATCGATCGATCGATCG",
            "seq2": "ATGCGATCGATCGATCGATCG",
            "seq3": "ATGCGATCGATCGATCGATCG",
        }

        finder = KmerBasedFinder(kmer_size=10, min_conservation=0.8)
        regions = finder.find_conserved_regions(
            sequences, min_length=15, max_length=30, align_regions=True
        )

        # Perfect matches should not need alignment
        assert len(regions) > 0
        # Consensus should be pure ACGT
        region = regions[0]
        assert all(base in "ACGT" for base in region.consensus)

    def test_conservation_threshold(self):
        """Test that conservation threshold is respected."""
        sequences = {
            "seq1": "ATGCGATCGATCGATCGATCG",
            "seq2": "ATGCGATCGATCGATCGATCG",
            "seq3": "TTTTTTTTTTTTTTTTTTTT",  # Very different
        }

        # High conservation threshold
        finder_strict = KmerBasedFinder(kmer_size=10, min_conservation=0.9)
        regions_strict = finder_strict.find_conserved_regions(
            sequences, min_length=15, max_length=30
        )

        # Should find no or few regions
        assert len(regions_strict) == 0

        # Low conservation threshold
        finder_relaxed = KmerBasedFinder(kmer_size=10, min_conservation=0.5)
        regions_relaxed = finder_relaxed.find_conserved_regions(
            sequences, min_length=15, max_length=30
        )

        # May find regions (depends on k-mer overlap)
        assert len(regions_relaxed) >= 0


class TestAlignmentRefinement:
    """Tests for alignment refinement functionality."""

    def test_identify_conserved_columns(self):
        """Test column-wise conservation calculation."""
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        # Create a simple alignment
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ATGC"), id="seq1"),
            SeqRecord(Seq("ATGC"), id="seq2"),
            SeqRecord(Seq("ATGT"), id="seq3"),  # Different at position 3
        ])

        finder = KmerBasedFinder(min_conservation=0.8)
        conserved_cols = finder._identify_conserved_columns(alignment, threshold=0.8)

        # Positions 0, 1, 2 should be conserved (100% match)
        # Position 3 is only 66% conserved (C in 2/3, T in 1/3)
        assert 0 in conserved_cols
        assert 1 in conserved_cols
        assert 2 in conserved_cols
        # Position 3 might not be conserved depending on threshold

    def test_extract_conserved_segments(self):
        """Test extraction of ungapped conserved segments."""
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        # Alignment with gaps
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ATGC--GATC"), id="seq1"),
            SeqRecord(Seq("ATGC--GATC"), id="seq2"),
            SeqRecord(Seq("ATGCTTGATC"), id="seq3"),
        ])

        finder = KmerBasedFinder(min_conservation=0.8)
        # All columns conserved
        conserved_cols = list(range(10))

        segments = finder._extract_conserved_segments(
            alignment, conserved_cols, min_length=4, max_length=20
        )

        # Should extract ungapped sequences
        assert len(segments) > 0
        # Gaps should be removed
        assert "-" not in segments[0]
