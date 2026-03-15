"""Tests for kmer_finder module."""

from ppdesign.kmer_finder import canonical_kmer, count_kmers


class TestCanonicalKmer:
    """Test canonical k-mer computation."""

    def test_canonical_kmer_basic(self):
        """Test basic canonical k-mer functionality."""
        assert canonical_kmer("ATCG") == "ATCG"
        assert canonical_kmer("CGAT") == "ATCG"  # Reverse complement

    def test_canonical_kmer_palindrome(self):
        """Test palindromic sequences."""
        assert canonical_kmer("ATAT") == "ATAT"
        assert canonical_kmer("GCGC") == "GCGC"

    def test_canonical_kmer_reverse(self):
        """Test that reverse is considered."""
        kmer = "AAAG"
        # rev = "GAAA"  # Not used in test
        # rev_comp = "CTTT"  # Not used in test
        assert canonical_kmer(kmer) == "AAAG"  # Original is smallest


class TestCountKmers:
    """Test k-mer counting functionality."""

    def test_count_kmers_simple(self):
        """Test basic k-mer counting."""
        sequence = "ATCGATCG"
        k_values = [3]
        counts = count_kmers(sequence, k_values)

        assert len(counts) == 1
        assert 3 in counts

        # Check specific k-mers
        # In sequence ATCGATCG:
        # Position 0: ATC -> ATC
        # Position 1: TCG -> CGA (reverse complement)
        # Position 2: CGA -> CGA
        # Position 3: GAT -> ATC (reverse complement)
        # Position 4: ATC -> ATC
        # Position 5: TCG -> CGA (reverse complement)
        kmer_counts = counts[3]
        assert (
            kmer_counts[canonical_kmer("ATC")] == 3
        )  # ATC appears at positions 0, 3, 4
        assert (
            kmer_counts[canonical_kmer("TCG")] == 3
        )  # TCG->CGA appears at positions 1, 2, 5

    def test_count_kmers_with_n(self):
        """Test that k-mers with N are skipped."""
        sequence = "ATCNATCG"
        k_values = [3]
        counts = count_kmers(sequence, k_values)

        kmer_counts = counts[3]
        # "TCN" and "CNA" and "NAT" should be skipped
        # Valid k-mers in ATCNATCG:
        # Position 0: ATC -> ATC
        # Position 1: TCN -> skipped (has N)
        # Position 2: CNA -> skipped (has N)
        # Position 3: NAT -> skipped (has N)
        # Position 4: ATC -> ATC
        # Position 5: TCG -> CGA (reverse complement)
        assert "TCN" not in str(kmer_counts)
        assert kmer_counts[canonical_kmer("ATC")] == 2  # Appears at positions 0 and 4

    def test_count_kmers_multiple_k(self):
        """Test counting with multiple k values."""
        sequence = "ATCGATCG"
        k_values = [2, 3, 4]
        counts = count_kmers(sequence, k_values)

        assert len(counts) == 3
        assert all(k in counts for k in k_values)

        # Check counts for different k values
        assert len(counts[2]) > 0
        assert len(counts[3]) > 0
        assert len(counts[4]) > 0
