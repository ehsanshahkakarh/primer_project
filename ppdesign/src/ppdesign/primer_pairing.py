"""Pair primers efficiently using position-based indexing."""

from typing import List, Dict
from .primer_types import Primer, PrimerPair


class PrimerPairMatcher:
    """Pair forward and reverse primers within amplicon size constraints."""

    def __init__(
        self,
        amplicon_min_size: int = 100,
        amplicon_max_size: int = 2000,
        tm_difference_max: float = 5.0,
    ):
        """
        Initialize primer pair matcher.

        Args:
            amplicon_min_size: Minimum amplicon size in bp
            amplicon_max_size: Maximum amplicon size in bp (supports long-read)
            tm_difference_max: Maximum Tm difference between primers (°C)

        Raises:
            ValueError: If constraints are invalid
        """
        if amplicon_min_size >= amplicon_max_size:
            raise ValueError("amplicon_min_size must be < amplicon_max_size")

        self.amplicon_min_size = amplicon_min_size
        self.amplicon_max_size = amplicon_max_size
        self.tm_difference_max = tm_difference_max

    def pair_primers(
        self, forward_primers: List[Primer], reverse_primers: List[Primer]
    ) -> List[PrimerPair]:
        """
        Pair primers efficiently using position indexing (O(N) performance).

        Algorithm:
        1. Index reverse primers by position for O(1) lookup
        2. For each forward primer, only check reverse primers in valid amplicon range
        3. Create and validate pairs meeting all constraints

        Args:
            forward_primers: List of forward primers
            reverse_primers: List of reverse primers

        Returns:
            List of valid primer pairs
        """
        # Index reverse primers by position for fast lookup
        reverse_index = self._index_by_position(reverse_primers)

        pairs = []
        for fwd in forward_primers:
            # Calculate valid position range for reverse primers
            min_pos = fwd.position + self.amplicon_min_size
            max_pos = fwd.position + self.amplicon_max_size

            # Only scan reverse primers within amplicon size window
            for rev_pos in range(min_pos, max_pos + 1):
                if rev_pos not in reverse_index:
                    continue

                # Check all reverse primers at this position
                for rev in reverse_index[rev_pos]:
                    pair = self._create_pair(fwd, rev)
                    if self._is_valid_pair(pair):
                        pairs.append(pair)

        return pairs

    def _index_by_position(self, primers: List[Primer]) -> Dict[int, List[Primer]]:
        """
        Index primers by position for O(1) lookup.

        Args:
            primers: List of primers to index

        Returns:
            Dictionary mapping position to list of primers at that position
        """
        index: Dict[int, List[Primer]] = {}
        for primer in primers:
            index.setdefault(primer.position, []).append(primer)
        return index

    def _create_pair(self, forward: Primer, reverse: Primer) -> PrimerPair:
        """
        Create primer pair with calculated metrics.

        Args:
            forward: Forward primer
            reverse: Reverse primer

        Returns:
            PrimerPair object with calculated amplicon size, Tm difference, etc.
        """
        amplicon_size = reverse.position - forward.position
        tm_diff = abs(forward.tm - reverse.tm)
        cross_dimer = self._calculate_cross_dimer(forward.sequence, reverse.sequence)

        return PrimerPair(
            forward=forward,
            reverse=reverse,
            amplicon_size=amplicon_size,
            tm_difference=tm_diff,
            cross_dimer_score=cross_dimer,
            quality_score=0.0,  # Will be set by scorer
        )

    def _is_valid_pair(self, pair: PrimerPair) -> bool:
        """
        Check if primer pair meets basic constraints.

        Requirements:
        1. Amplicon size within range
        2. Tm difference within threshold
        3. Forward and reverse primers share at least one target sequence

        Args:
            pair: Primer pair to validate

        Returns:
            True if valid, False otherwise
        """
        # Check amplicon size
        if not (
            self.amplicon_min_size <= pair.amplicon_size <= self.amplicon_max_size
        ):
            return False

        # Check Tm difference
        if pair.tm_difference > self.tm_difference_max:
            return False

        # Check target sequence overlap
        fwd_targets = set(pair.forward.target_sequences)
        rev_targets = set(pair.reverse.target_sequences)
        common_targets = fwd_targets & rev_targets

        if not common_targets:
            return False  # No shared targets = invalid pair

        return True

    def _calculate_cross_dimer(self, seq1: str, seq2: str) -> int:
        """
        Calculate cross-dimer score (simple complementarity check).

        Checks for complementary base pairs between forward and reverse primers.

        Args:
            seq1: First primer sequence
            seq2: Second primer sequence

        Returns:
            Number of complementary base pairs when aligned in reverse
        """
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

        # Check forward-reverse complementarity
        # Align seq1 with reverse of seq2
        score = 0
        for b1, b2 in zip(seq1, seq2[::-1]):
            if complement.get(b1) == b2:
                score += 1

        return score
