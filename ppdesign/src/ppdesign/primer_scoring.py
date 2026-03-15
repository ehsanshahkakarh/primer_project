"""Score and rank primer pairs by quality."""

from typing import List
from .primer_types import PrimerPair


class PrimerPairScorer:
    """Score and rank primer pairs based on multiple quality metrics."""

    def __init__(self, tm_difference_max: float = 5.0):
        """
        Initialize primer pair scorer.

        Args:
            tm_difference_max: Maximum Tm difference for normalization
        """
        self.tm_difference_max = tm_difference_max

    def score_pairs(self, pairs: List[PrimerPair]) -> List[PrimerPair]:
        """
        Score all pairs and sort by quality (descending).

        Scoring components:
        - Conservation (0-40 points): Average conservation of both primers
        - Tm matching (0-30 points): Closer Tm difference is better
        - GC content (0-20 points): Closer to 50% is optimal
        - Cross-dimer (0-10 points): Lower dimer score is better

        Args:
            pairs: List of primer pairs to score

        Returns:
            Sorted list with quality_score set for each pair
        """
        # Calculate scores
        for pair in pairs:
            pair.quality_score = self._calculate_score(pair)

        # Sort by quality score (descending)
        return sorted(pairs, key=lambda p: p.quality_score, reverse=True)

    def _calculate_score(self, pair: PrimerPair) -> float:
        """
        Calculate composite quality score (0-100 scale).

        Scoring breakdown:
        - 35% weight: Conservation
        - 25% weight: Tm matching
        - 15% weight: GC content
        - 15% weight: Hairpin stability (primer3 ΔG)
        - 10% weight: Cross-dimer

        Args:
            pair: Primer pair to score

        Returns:
            Quality score from 0-100
        """
        score = 0.0

        # Conservation score (0-35 points)
        score += pair.avg_conservation * 35.0

        # Tm matching score (0-25 points)
        tm_penalty = min(pair.tm_difference / self.tm_difference_max, 1.0)
        score += (1.0 - tm_penalty) * 25.0

        # GC content score (0-15 points)
        gc_distance = abs(pair.avg_gc - 50.0) / 50.0
        score += (1.0 - gc_distance) * 15.0

        # Hairpin stability score (0-15 points)
        # ΔG > -2 kcal/mol = no concern (full points)
        # ΔG < -9 kcal/mol = very stable hairpin (zero points)
        score += self._hairpin_score(pair.forward.hairpin_dg, pair.reverse.hairpin_dg)

        # Cross-dimer score (0-10 points)
        dimer_penalty = min(pair.cross_dimer_score / 10.0, 1.0)
        score += (1.0 - dimer_penalty) * 10.0

        return score

    @staticmethod
    def _hairpin_score(fwd_dg, rev_dg) -> float:
        """Score hairpin stability (0-15 points). Needs primer3 ΔG values."""
        if fwd_dg is None or rev_dg is None:
            return 15.0  # No data → assume no hairpin concern

        def _single(dg):
            # dG > -2: no hairpin concern → 1.0
            # dG < -9: strong hairpin → 0.0
            if dg >= -2.0:
                return 1.0
            if dg <= -9.0:
                return 0.0
            return (dg + 9.0) / 7.0

        avg = (_single(fwd_dg) + _single(rev_dg)) / 2.0
        return avg * 15.0

    def get_top_pairs(
        self, pairs: List[PrimerPair], n: int = 10
    ) -> List[PrimerPair]:
        """
        Get top N scoring primer pairs.

        Args:
            pairs: List of scored primer pairs
            n: Number of pairs to return

        Returns:
            Top N pairs by quality score
        """
        scored_pairs = self.score_pairs(pairs)
        return scored_pairs[:n]

    def filter_by_score(
        self, pairs: List[PrimerPair], min_score: float
    ) -> List[PrimerPair]:
        """
        Filter pairs by minimum quality score.

        Args:
            pairs: List of primer pairs
            min_score: Minimum quality score threshold

        Returns:
            Pairs with quality_score >= min_score
        """
        scored_pairs = self.score_pairs(pairs)
        return [pair for pair in scored_pairs if pair.quality_score >= min_score]
