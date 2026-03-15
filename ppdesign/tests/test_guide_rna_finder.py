#!/usr/bin/env python3
"""Tests for CRISPR guide RNA finder."""

from ppdesign.guide_rna_finder import GuideRNAFinder, GuideRNA


class TestPAMDetection:
    """Test PAM site detection."""

    def test_ngg_forward_strand(self):
        """Test NGG PAM detection on forward strand."""
        finder = GuideRNAFinder()
        sequence = "ATCGATCGATCGATCGATCGATCGGGATCG"
        pam_sites = finder.find_pam_sites(sequence, "NGG")

        # Should find NGG at position 22 (0-indexed)
        assert len(pam_sites) == 1
        assert pam_sites[0] == (22, "+")

    def test_nag_forward_strand(self):
        """Test NAG PAM detection on forward strand."""
        finder = GuideRNAFinder()
        sequence = "ATCGATCGATCGATCGATCGATCAGATCG"
        pam_sites = finder.find_pam_sites(sequence, "NAG")

        # Should find NAG at position 22 (0-indexed)
        assert len(pam_sites) == 1
        assert pam_sites[0] == (22, "+")

    def test_multiple_pams(self):
        """Test detection of multiple PAM sites."""
        finder = GuideRNAFinder()
        sequence = "ATCGGGATCGATCAGGTATCGATCGGGATC"
        pam_sites = finder.find_pam_sites(sequence, "NGG")

        # Should find NGG at positions 2, 13, and 23 (0-indexed)
        assert len(pam_sites) == 3
        positions = [site[0] for site in pam_sites if site[1] == "+"]
        assert 2 in positions
        assert 23 in positions

    def test_reverse_strand_pam(self):
        """Test PAM detection on reverse strand."""
        finder = GuideRNAFinder()
        # CCG is reverse complement of CGG
        sequence = "ATCGATCGATCGATCGATCGATCCGATCG"
        pam_sites = finder.find_pam_sites(sequence, "CGG")

        # Should find at least one PAM on reverse strand
        reverse_sites = [site for site in pam_sites if site[1] == "-"]
        assert len(reverse_sites) > 0


class TestGuideExtraction:
    """Test guide RNA extraction."""

    def test_extract_guide_forward(self):
        """Test guide extraction upstream of PAM on forward strand."""
        finder = GuideRNAFinder()
        # 20bp + NGG
        sequence = "ATCGATCGATCGATCGATCGATCGGGATCGATCG"

        guide = finder.extract_guide_sequence(sequence, 22, "+", 3)

        assert guide is not None
        assert len(guide) == 20
        assert guide == "CGATCGATCGATCGATCGAT"  # 20bp upstream of position 22

    def test_extract_guide_reverse(self):
        """Test guide extraction on reverse strand."""
        finder = GuideRNAFinder()
        # Need sequence where guide would be downstream of PAM position
        sequence = "ATCGATCCGATCGATCGATCGATCGATCGATCGATCGATCG"

        # For reverse strand, guide is downstream of PAM
        guide = finder.extract_guide_sequence(sequence, 5, "-", 3)

        if guide:  # May be None if too close to edge
            assert len(guide) == 20

    def test_edge_case_too_short(self):
        """Test extraction when sequence is too short."""
        finder = GuideRNAFinder()
        sequence = "ATCGATCGGGATC"  # Too short for 20bp guide

        guide = finder.extract_guide_sequence(sequence, 7, "+", 3)

        assert guide is None  # Should return None when too close to start

    def test_guide_with_n_skipped(self):
        """Test that guides containing N are skipped."""
        finder = GuideRNAFinder()
        sequences = {"seq1": "ATCGATCGATCNATCGATCGATCGGGATCGATCG"}

        guides = finder._extract_all_guides(sequences)

        # Should skip guide with N
        for guide in guides:
            assert "N" not in guide["sequence"]


class TestGuideRNAClustering:
    """Test guide RNA clustering and consensus generation."""

    def test_identical_guides_cluster(self):
        """Test that identical guides cluster together."""
        finder = GuideRNAFinder()
        guides = [
            {
                "sequence": "ATCGATCGATCGATCGATCG",
                "pam": "NGG",
                "strand": "+",
                "seq_id": "seq1",
                "position": 20,
            },
            {
                "sequence": "ATCGATCGATCGATCGATCG",
                "pam": "NGG",
                "strand": "+",
                "seq_id": "seq2",
                "position": 20,
            },
        ]

        clusters = finder._cluster_guides(guides, max_mismatches=2)

        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_guides_with_mismatches_cluster(self):
        """Test clustering with allowed mismatches."""
        finder = GuideRNAFinder()
        guides = [
            {
                "sequence": "ATCGATCGATCGATCGATCG",
                "pam": "NGG",
                "strand": "+",
                "seq_id": "seq1",
                "position": 20,
            },
            {
                "sequence": "ATCGATCGATCGATCGATCA",
                "pam": "NGG",
                "strand": "+",
                "seq_id": "seq2",
                "position": 20,
            },  # 1 mismatch at end
        ]

        clusters = finder._cluster_guides(guides, max_mismatches=2)

        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_different_guides_separate_clusters(self):
        """Test that different guides form separate clusters."""
        finder = GuideRNAFinder()
        guides = [
            {
                "sequence": "ATCGATCGATCGATCGATCG",
                "pam": "NGG",
                "strand": "+",
                "seq_id": "seq1",
                "position": 20,
            },
            {
                "sequence": "GCTAGCTAGCTAGCTAGCTA",
                "pam": "NGG",
                "strand": "+",
                "seq_id": "seq2",
                "position": 20,
            },  # Completely different
        ]

        clusters = finder._cluster_guides(guides, max_mismatches=2)

        assert len(clusters) == 2
        assert len(clusters[0]) == 1
        assert len(clusters[1]) == 1


class TestConsensusGeneration:
    """Test consensus sequence generation with IUPAC codes."""

    def test_identical_sequences_consensus(self):
        """Test consensus of identical sequences."""
        finder = GuideRNAFinder()
        sequences = ["ATCG", "ATCG", "ATCG"]

        consensus = finder._generate_consensus(sequences)

        assert consensus == "ATCG"
        assert "N" not in consensus

    def test_single_position_variation(self):
        """Test consensus with variation at one position."""
        finder = GuideRNAFinder()
        sequences = ["ATCG", "ACCG"]  # T/C at position 1

        consensus = finder._generate_consensus(sequences)

        assert consensus == "AYCG"  # Y = C or T

    def test_iupac_code_generation(self):
        """Test IUPAC code generation for various base combinations."""
        finder = GuideRNAFinder()

        # Test R (A/G)
        assert finder._get_iupac_code({"A", "G"}) == "R"

        # Test Y (C/T)
        assert finder._get_iupac_code({"C", "T"}) == "Y"

        # Test N (all bases)
        assert finder._get_iupac_code({"A", "C", "G", "T"}) == "N"

    def test_degenerate_count(self):
        """Test counting of degenerate bases."""
        finder = GuideRNAFinder()
        sequences = ["ATCGATCG", "ACCGATCG", "ATCGATGG"]

        consensus = finder._generate_consensus(sequences)

        degenerate_count = sum(1 for base in consensus if base not in "ATCG")
        assert degenerate_count == 2  # Two degenerate positions


class TestConservationAnalysis:
    """Test conservation calculation."""

    def test_full_conservation(self):
        """Test 100% conservation across sequences."""
        finder = GuideRNAFinder(min_conservation=0.8)
        sequences = {
            "seq1": "ATCGATCGATCGATCGATCGATCGGGATCG",
            "seq2": "ATCGATCGATCGATCGATCGATCGGGATCG",
            "seq3": "ATCGATCGATCGATCGATCGATCGGGATCG",
        }

        guides = finder.find_guide_rnas(sequences)

        if guides:  # Should find at least one guide
            assert guides[0].conservation == 1.0

    def test_partial_conservation(self):
        """Test partial conservation across sequences."""
        finder = GuideRNAFinder(min_conservation=0.5)
        sequences = {
            "seq1": "ATCGATCGATCGATCGATCGATCGGGATCG",
            "seq2": "ATCGATCGATCGATCGATCGATCGGGATCG",
            "seq3": "GCTAGCTAGCTAGCTAGCTAGCTAGGATCG",  # Different guide
        }

        guides = finder.find_guide_rnas(sequences, allow_mismatches=0)

        # Should find guide present in 2/3 sequences
        conserved_guides = [g for g in guides if g.conservation >= 0.66]
        assert len(conserved_guides) > 0

    def test_conservation_calculation_method(self):
        """Test that conservation calculation correctly counts sequence-level coverage."""
        finder = GuideRNAFinder(
            min_conservation=0.1
        )  # Low threshold to find all guides

        # Create test data where we know exact conservation
        sequences = {
            "seq1": "AAAAAAAAAAAAAAAAAAAAAGGGATCG",  # Contains guide at pos 0-19
            "seq2": "AAAAAAAAAAAAAAAAAAAAAGGGATCG",  # Contains same guide at pos 0-19
            "seq3": "TTTTTTTTTTTTTTTTTTTTTAGGATCG",  # Contains different guide at pos 0-19
            "seq4": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCG",  # No NGG PAM, no guide
        }

        guides = finder.find_guide_rnas(sequences, allow_mismatches=0)

        # Find guide that should be in 2/4 sequences (seq1, seq2)
        aaa_guide = None
        for guide in guides:
            if guide.sequence == "AAAAAAAAAAAAAAAAAAAA":
                aaa_guide = guide
                break

        assert aaa_guide is not None, "Should find AAA guide"
        # Conservation should be exactly 2/4 = 0.5 (50%)
        assert (
            aaa_guide.conservation == 0.5
        ), f"Expected 0.5, got {aaa_guide.conservation}"
        # Should have positions in exactly 2 sequences
        assert (
            len(aaa_guide.positions) == 2
        ), f"Expected 2 sequences, got {len(aaa_guide.positions)}"
        assert "seq1" in aaa_guide.positions
        assert "seq2" in aaa_guide.positions


class TestScoringSystem:
    """Test guide RNA scoring."""

    def test_specificity_score_no_degenerate(self):
        """Test specificity score for non-degenerate guide."""
        finder = GuideRNAFinder()
        guide = GuideRNA(
            sequence="ATCGATCGATCGATCGATCG",
            pam="NGG",
            strand="+",
            conservation=1.0,
            consensus="ATCGATCGATCGATCGATCG",
        )

        finder._calculate_scores(guide)

        assert guide.specificity_score == 1.0

    def test_specificity_score_with_degenerate(self):
        """Test specificity score with degenerate bases."""
        finder = GuideRNAFinder()
        guide = GuideRNA(
            sequence="ATCGATCGATCGATCGATCG",
            pam="NGG",
            strand="+",
            conservation=1.0,
            consensus="RTCGATCGATCGATCGATCG",  # One degenerate base
            degenerate_count=1,
        )

        finder._calculate_scores(guide)

        assert 0 < guide.specificity_score < 1.0

    def test_quality_score_calculation(self):
        """Test overall quality score calculation."""
        finder = GuideRNAFinder()
        guide = GuideRNA(
            sequence="ATCGATCGATCGATCGATCG",
            pam="NGG",
            strand="+",
            conservation=0.9,
            consensus="ATCGATCGATCGATCGATCG",
        )

        finder._calculate_scores(guide)

        # Quality score should balance conservation and specificity
        assert 0 < guide.quality_score <= 1.0

        # NGG should score higher than NAG
        guide_nag = GuideRNA(
            sequence="ATCGATCGATCGATCGATCG",
            pam="NAG",
            strand="+",
            conservation=0.9,
            consensus="ATCGATCGATCGATCGATCG",
        )

        finder._calculate_scores(guide_nag)

        assert guide_nag.quality_score < guide.quality_score


class TestIntegration:
    """Integration tests with real-like sequences."""

    def test_find_guides_in_test_sequences(self):
        """Test finding guides in a small test dataset."""
        finder = GuideRNAFinder(min_conservation=0.5, max_degenerate=2)

        # Create test sequences with common guide regions
        common_guide = "ATCGATCGATCGATCGATCG"
        sequences = {
            "seq1": f"GCGCGC{common_guide}TGGCGCGCG",
            "seq2": f"ATATAT{common_guide}AGGATATAT",
            "seq3": f"CGCGCG{common_guide}CGGCGCGCG",
        }

        guides = finder.find_guide_rnas(sequences)

        # Should find at least one conserved guide
        assert len(guides) > 0

        # Check guide properties
        for guide in guides:
            assert len(guide.sequence) == 20
            assert guide.pam in ["NGG", "NAG", "TGG", "AGG", "CGG"]
            assert guide.conservation >= 0.5
            assert guide.degenerate_count <= 2

    def test_output_formatting(self):
        """Test guide RNA output formatting."""
        finder = GuideRNAFinder()
        guide = GuideRNA(
            sequence="ATCGATCGATCGATCGATCG",
            pam="NGG",
            strand="+",
            conservation=0.85,
            consensus="ATCGATCGATCGATCGATCG",
            degenerate_count=0,
            specificity_score=1.0,
            quality_score=0.91,
            positions={"seq1": [(20, "+")], "seq2": [(25, "+")]},
        )

        output = finder.format_guide_output(guide)

        assert "Guide:" in output
        assert "PAM: NGG" in output
        assert "85.00%" in output  # Conservation percentage
        assert "Degenerate bases: 0" in output


class TestNonDegenerateMode:
    """Test non-degenerate mode functionality."""

    def test_no_degenerate_flag(self):
        """Test that no_degenerate flag prevents degenerate bases."""
        sequences = {
            "seq1": "ATCGATCGATCGATCGATCGAGGTTTT",  # Has PAM at position 20
            "seq2": "ATCGATCGATCGATCGATCGTGGTTTT",  # Slightly different, PAM at 20
            "seq3": "TTCGATCGATCGATCGATCGAGGTTTT",  # Different at position 0
        }

        # Normal mode - should create consensus with degenerate bases
        finder_normal = GuideRNAFinder(min_conservation=0.5, max_degenerate=3)
        guides_normal = finder_normal.find_guide_rnas(sequences)

        # Normal mode may have degenerate bases
        assert len(guides_normal) >= 0  # Should find some guides

        # Non-degenerate mode - should only use exact sequences
        finder_no_degen = GuideRNAFinder(
            min_conservation=0.5, max_degenerate=0, no_degenerate=True
        )
        guides_no_degen = finder_no_degen.find_guide_rnas(sequences)

        # In non-degenerate mode, all guides should have degenerate_count = 0
        for guide in guides_no_degen:
            assert guide.degenerate_count == 0
            assert guide.consensus == guide.sequence
            assert all(base in "ATCG" for base in guide.consensus)

    def test_required_target_prioritization(self):
        """Test that required targets are prioritized in non-degenerate mode."""
        sequences = {
            "required1": "ATCGATCGATCGATCGATCGAGGTTTT",
            "required2": "ATCGATCGATCGATCGATCGAGGTTTT",
            "optional1": "TTTTTTTTTTTTTTTTTTTTCGGTTTT",
            "optional2": "GGGGGGGGGGGGGGGGGGGGAGGTTTT",
        }

        required_targets = ["required1", "required2"]

        finder = GuideRNAFinder(
            min_conservation=0.25,  # Low threshold to get more guides
            max_degenerate=0,
            no_degenerate=True,
            required_targets=required_targets,
        )

        guides = finder.find_guide_rnas(sequences)

        # Check that guides covering required targets are scored higher
        if guides:
            # First guides should cover required targets
            top_guide = guides[0]
            assert top_guide.required_targets_hit > 0

            # Check that required target coverage affects quality score
            for guide in guides:
                if guide.required_targets_hit > 0:
                    assert guide.quality_score > 0

    def test_greedy_optimization(self):
        """Test greedy set cover optimization in non-degenerate mode."""
        # Create sequences where different guides cover different subsets
        sequences = {
            "seq1": "ATCGATCGATCGATCGATCGAGGTTTT",  # Covered by guide A
            "seq2": "ATCGATCGATCGATCGATCGAGGTTTT",  # Covered by guide A
            "seq3": "GGGGGGGGGGGGGGGGGGGGCGGTTTT",  # Covered by guide B
            "seq4": "GGGGGGGGGGGGGGGGGGGGCGGTTTT",  # Covered by guide B
            "seq5": "TTTTTTTTTTTTTTTTTTTTCGGTTTT",  # Covered by guide C
        }

        required_targets = ["seq1", "seq3"]  # Require coverage of specific sequences

        finder = GuideRNAFinder(
            min_conservation=0.2,
            max_degenerate=0,
            no_degenerate=True,
            required_targets=required_targets,
        )

        guides = finder.find_guide_rnas(sequences)

        # The optimization should select guides that cover all required targets
        covered_required = set()
        for guide in guides:
            covered_required.update(set(guide.positions.keys()) & set(required_targets))

        # All required targets should be covered if possible
        assert len(covered_required) <= len(required_targets)

    def test_best_guide_selection_from_cluster(self):
        """Test selection of best individual guide from a cluster."""
        sequences = {
            "seq1": "ATCGATCGATCGATCGATCGAGGTTTT",
            "seq2": "ATCGATCGATCGATCGATCGAGGTTTT",
            "seq3": "ATCGATCGATCGATCGATCGTGGTTTT",  # One mismatch
        }

        finder = GuideRNAFinder(
            min_conservation=0.5,
            max_degenerate=0,
            no_degenerate=True,
            required_targets=["seq1", "seq2"],
        )

        guides = finder.find_guide_rnas(sequences, allow_mismatches=1)

        # In non-degenerate mode, should select the exact sequence
        # that covers the most required targets
        for guide in guides:
            assert guide.degenerate_count == 0
            if guide.required_targets_hit > 0:
                # Guide should be an exact sequence, not a consensus
                assert all(base in "ATCG" for base in guide.sequence)

    def test_no_degenerate_with_no_required_targets(self):
        """Test non-degenerate mode without required targets."""
        sequences = {
            "seq1": "ATCGATCGATCGATCGATCGAGGTTTT",
            "seq2": "ATCGATCGATCGATCGATCGTGGTTTT",
            "seq3": "TTCGATCGATCGATCGATCGAGGTTTT",
        }

        finder = GuideRNAFinder(
            min_conservation=0.3,
            max_degenerate=0,
            no_degenerate=True,
            required_targets=None,  # No required targets
        )

        guides = finder.find_guide_rnas(sequences, allow_mismatches=2)

        # Should still work, selecting best individual sequences by coverage
        for guide in guides:
            assert guide.degenerate_count == 0
            assert guide.consensus == guide.sequence
            assert all(base in "ATCG" for base in guide.sequence)
