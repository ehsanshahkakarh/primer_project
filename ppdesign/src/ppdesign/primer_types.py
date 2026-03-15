"""Data structures for primer pair design."""

from dataclasses import dataclass
from typing import List, Tuple, Optional


@dataclass
class Primer:
    """
    Single primer (forward or reverse).

    Coordinate system:
    - 0-based inclusive start positions
    - Forward primers: sequence as-is, position = start
    - Reverse primers: stored as reverse-complement, position = binding start on forward strand

    Tm calculation: BioPython MeltingTemp.Tm_NN (Na=50mM, dnac1=250nM)
    Hairpin dG: primer3-py if available, else None
    """

    sequence: str  # 18-25bp, ACGT only
    position: int  # Start position (0-based)
    strand: str  # '+' or '-'
    tm: float  # Melting temperature (°C)
    gc_content: float  # Percentage 0-100
    hairpin_dg: Optional[float]  # ΔG kcal/mol (negative = stable), None if unavailable
    self_dimer_score: int  # Simple complementarity count
    conservation: float  # Fraction 0.0-1.0
    target_sequences: List[str]  # Sequence IDs where this primer is found


@dataclass
class PrimerPair:
    """Forward + reverse primer pair for PCR amplification."""

    forward: Primer
    reverse: Primer
    amplicon_size: int  # reverse.position - forward.position
    tm_difference: float  # abs(forward.tm - reverse.tm)
    cross_dimer_score: int  # Complementarity between primers
    quality_score: float  # Composite score (0-100)

    @property
    def amplicon_range(self) -> Tuple[int, int]:
        """
        Return (start, end) of amplicon on forward strand.

        Returns:
            Tuple of (forward_start, reverse_end)
        """
        return (self.forward.position, self.reverse.position)

    @property
    def avg_tm(self) -> float:
        """Return average melting temperature of the pair."""
        return (self.forward.tm + self.reverse.tm) / 2.0

    @property
    def avg_gc(self) -> float:
        """Return average GC content of the pair."""
        return (self.forward.gc_content + self.reverse.gc_content) / 2.0

    @property
    def avg_conservation(self) -> float:
        """Return average conservation of the pair."""
        return (self.forward.conservation + self.reverse.conservation) / 2.0
