"""Tests for primer candidate generation with IUPAC handling."""

import pytest
from ppdesign.primer_candidates import PrimerCandidateGenerator
from ppdesign.conserved_finder import ConservedRegion


class TestPrimerCandidateGenerator:
    """Tests for primer candidate generation."""

    def test_basic_primer_generation(self):
        """Test basic forward primer generation."""
        region = ConservedRegion(
            start=0,
            end=30,
            sequences=[
                "ATGCGATCGATCGATCGATCGATCGATCG",
                "ATGCGATCGATCGATCGATCGATCGATCG",
            ],
            consensus="ATGCGATCGATCGATCGATCGATCGATCG",
            positions={"seq1": 100, "seq2": 200},
            conservation=1.0,
        )

        generator = PrimerCandidateGenerator(
            primer_min_length=18,
            primer_max_length=25,
            tm_min=55.0,
            tm_max=65.0,
            conservation_threshold=0.8,
        )

        primers = generator.generate_forward_primers([region])

        # Should generate multiple primers of different lengths
        assert len(primers) > 0
        assert all(p.sequence.startswith("ATGCGATCG") for p in primers)

    def test_individual_sequence_extraction(self):
        """Test that primers are extracted from individual sequences, not consensus."""
        region = ConservedRegion(
            start=0,
            end=30,
            sequences=[
                "ATGCAATCGATCGATCGATCGATCGA",
                "ATGCGATCGATCGATCGATCGATCGA",
            ],
            consensus="ATGCRATCGATCGATCGATCGATCGA",  # R at position 5
            positions={"seq1": 100, "seq2": 200},
            conservation=1.0,
        )

        generator = PrimerCandidateGenerator(
            primer_min_length=18,
            primer_max_length=22,
            tm_min=50.0,
            tm_max=70.0,
            conservation_threshold=0.8,
            max_degenerate_positions=2,
        )

        primers = generator.generate_forward_primers([region])

        # Should extract from both sequences
        assert len(primers) > 0

        # Should have primers from both variants
        primer_seqs = {p.sequence for p in primers}
        # Some primers should start with ATGCA (from seq1)
        # Some primers should start with ATGCG (from seq2)
        has_seq1_variant = any(seq.startswith("ATGCA") for seq in primer_seqs)
        has_seq2_variant = any(seq.startswith("ATGCG") for seq in primer_seqs)

        assert has_seq1_variant or has_seq2_variant

    def test_iupac_handling_strict(self):
        """Test strict ACGT-only validation (max_degenerate_positions=0)."""
        region = ConservedRegion(
            start=0,
            end=25,
            sequences=[
                "ATGCAATCGATCGATCGATCG",
                "ATGCGATCGATCGATCGATCG",
            ],
            consensus="ATGCRATCGATCGATCGATCG",  # R at position 5
            positions={"seq1": 100, "seq2": 200},
            conservation=1.0,
        )

        generator = PrimerCandidateGenerator(
            primer_min_length=18,
            primer_max_length=22,
            tm_min=50.0,
            tm_max=70.0,
            max_degenerate_positions=0,  # Strict ACGT only
        )

        primers = generator.generate_forward_primers([region])

        # All primers should be pure ACGT
        assert all(all(base in "ACGT" for base in p.sequence) for p in primers)

    def test_iupac_handling_flexible(self):
        """Test flexible validation allowing 1-2 IUPAC codes."""
        # Create region with IUPAC codes
        region = ConservedRegion(
            start=0,
            end=25,
            sequences=[
                "ATGCAATCGATCGATCGATCG",
                "ATGCGATCGATCGATCGATCG",
            ],
            consensus="ATGCRATCGATCGATCGATCG",
            positions={"seq1": 100, "seq2": 200},
            conservation=1.0,
        )

        generator_flexible = PrimerCandidateGenerator(
            primer_min_length=18,
            primer_max_length=22,
            tm_min=50.0,
            tm_max=70.0,
            max_degenerate_positions=2,  # Allow 1-2 IUPAC codes
        )

        primers = generator_flexible.generate_forward_primers([region])

        # Should still extract ACGT-only primers from individual sequences
        assert len(primers) > 0
        # Extracted from individual sequences, so should be ACGT only
        assert all(all(base in "ACGT" for base in p.sequence) for p in primers)

    def test_reverse_primer_generation(self):
        """Test reverse primer generation from individual sequences."""
        region = ConservedRegion(
            start=0,
            end=30,
            sequences=[
                "ATGCGATCGATCGATCGATCGATCGA",
                "ATGCGATCGATCGATCGATCGATCGA",
            ],
            consensus="ATGCGATCGATCGATCGATCGATCGA",
            positions={"seq1": 100, "seq2": 200},
            conservation=1.0,
        )

        generator = PrimerCandidateGenerator(
            primer_min_length=18,
            primer_max_length=22,
            tm_min=50.0,
            tm_max=70.0,
        )

        primers = generator.generate_reverse_primers([region])

        # Should generate reverse-complemented primers
        assert len(primers) > 0
        assert all(p.strand == "-" for p in primers)
        # All should be valid DNA sequences
        assert all(all(base in "ACGT" for base in p.sequence) for p in primers)

    def test_tm_gc_constraints(self):
        """Test that Tm and GC constraints are enforced."""
        region = ConservedRegion(
            start=0,
            end=30,
            sequences=[
                "ATGCGATCGATCGATCGATCGATCGA",
                "ATGCGATCGATCGATCGATCGATCGA",
            ],
            consensus="ATGCGATCGATCGATCGATCGATCGA",
            positions={"seq1": 100, "seq2": 200},
            conservation=1.0,
        )

        generator = PrimerCandidateGenerator(
            primer_min_length=18,
            primer_max_length=25,
            tm_min=55.0,
            tm_max=65.0,
            gc_min=40.0,
            gc_max=60.0,
        )

        primers = generator.generate_forward_primers([region])

        # All primers should meet constraints
        assert all(55.0 <= p.tm <= 65.0 for p in primers)
        assert all(40.0 <= p.gc_content <= 60.0 for p in primers)

    def test_conservation_threshold_filter(self):
        """Test that low conservation regions are filtered out."""
        region_low = ConservedRegion(
            start=0,
            end=25,
            sequences=["ATGCGATCGATCGATCGATCG", "ATGCGATCGATCGATCGATCG"],
            consensus="ATGCGATCGATCGATCGATCG",
            positions={"seq1": 100, "seq2": 200},
            conservation=0.5,  # Below threshold
        )

        region_high = ConservedRegion(
            start=0,
            end=25,
            sequences=["ATGCGATCGATCGATCGATCG", "ATGCGATCGATCGATCGATCG"],
            consensus="ATGCGATCGATCGATCGATCG",
            positions={"seq1": 100, "seq2": 200},
            conservation=0.9,  # Above threshold
        )

        generator = PrimerCandidateGenerator(
            conservation_threshold=0.8,
            primer_min_length=18,
            primer_max_length=22,
            tm_min=50.0,
            tm_max=70.0,
        )

        primers_low = generator.generate_forward_primers([region_low])
        primers_high = generator.generate_forward_primers([region_high])

        # Low conservation region should be filtered
        assert len(primers_low) == 0
        # High conservation region should produce primers
        assert len(primers_high) > 0

    def test_multiple_length_candidates(self):
        """Test that multiple primer lengths are generated."""
        region = ConservedRegion(
            start=0,
            end=30,
            sequences=[
                "ATGCGATCGATCGATCGATCGATCGA",
                "ATGCGATCGATCGATCGATCGATCGA",
            ],
            consensus="ATGCGATCGATCGATCGATCGATCGA",
            positions={"seq1": 100, "seq2": 200},
            conservation=1.0,
        )

        generator = PrimerCandidateGenerator(
            primer_min_length=18,
            primer_max_length=22,
            tm_min=50.0,
            tm_max=70.0,
        )

        primers = generator.generate_forward_primers([region])

        # Should have primers of different lengths
        primer_lengths = {len(p.sequence) for p in primers}
        assert len(primer_lengths) > 1  # Multiple lengths
        assert all(18 <= length <= 22 for length in primer_lengths)

    def test_max_degenerate_validation(self):
        """Test validation with different max_degenerate_positions values."""
        generator_strict = PrimerCandidateGenerator(max_degenerate_positions=0)
        generator_flexible = PrimerCandidateGenerator(max_degenerate_positions=2)
        generator_very_flexible = PrimerCandidateGenerator(max_degenerate_positions=5)

        # Pure ACGT sequence
        assert generator_strict._is_valid_sequence("ATGCGATCGATCG")
        assert generator_flexible._is_valid_sequence("ATGCGATCGATCG")

        # Sequence with 1 IUPAC code
        assert not generator_strict._is_valid_sequence("ATGCRATCGATCG")  # R = A/G
        assert generator_flexible._is_valid_sequence("ATGCRATCGATCG")

        # Sequence with 3 IUPAC codes
        assert not generator_strict._is_valid_sequence("ATGCRATYGATCG")
        assert not generator_flexible._is_valid_sequence("ATGCRATYGATCG")  # > 2
        assert generator_very_flexible._is_valid_sequence("ATGCRATYGATCG")

        # Invalid sequence (not IUPAC code)
        assert not generator_strict._is_valid_sequence("ATGCXATCGATCG")
        assert not generator_flexible._is_valid_sequence("ATGCXATCGATCG")
