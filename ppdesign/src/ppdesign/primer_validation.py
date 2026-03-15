"""Validate primers against PCR best practices."""

import re
from typing import Tuple, List
from .primer_types import Primer, PrimerPair


class PrimerValidator:
    """Validate primers and primer pairs against PCR best practices."""

    def __init__(
        self,
        max_poly_x: int = 4,
        require_3_prime_gc_clamp: bool = True,
        max_cross_dimer_score: int = 8,
    ):
        """
        Initialize primer validator.

        Args:
            max_poly_x: Maximum allowed poly-X run length
            require_3_prime_gc_clamp: Require G or C in last 5 bases
            max_cross_dimer_score: Maximum allowed cross-dimer score
        """
        self.max_poly_x = max_poly_x
        self.require_3_prime_gc_clamp = require_3_prime_gc_clamp
        self.max_cross_dimer_score = max_cross_dimer_score

    def validate_primer(self, primer: Primer) -> Tuple[bool, List[str]]:
        """
        Validate single primer against PCR best practices.

        Checks:
        - Sequence composition (ACGT only)
        - Poly-X runs (no long repeats)
        - 3' GC clamp (stability)
        - 3' end stability (avoid AT-rich ends)

        Args:
            primer: Primer to validate

        Returns:
            Tuple of (is_valid, list_of_issues)
        """
        issues = []

        # Check sequence composition
        if not re.match(r"^[ACGT]+$", primer.sequence):
            issues.append("Non-ACGT characters detected")

        # Check for poly-X runs
        if self._has_poly_x_run(primer.sequence):
            issues.append(f"Poly-X run longer than {self.max_poly_x} bases")

        # Check 3' GC clamp
        if self.require_3_prime_gc_clamp and not self._has_gc_clamp(primer.sequence):
            issues.append("Missing 3' GC clamp (no G/C in last 5 bases)")

        # Check 3' stability
        if self._has_3_prime_instability(primer.sequence):
            issues.append("3' end instability (AT-rich last 3 bases)")

        return (len(issues) == 0, issues)

    def validate_pair(self, pair: PrimerPair) -> Tuple[bool, List[str]]:
        """
        Validate primer pair.

        Checks:
        - Individual primer validation
        - Cross-dimer formation

        Args:
            pair: Primer pair to validate

        Returns:
            Tuple of (is_valid, list_of_issues)
        """
        issues = []

        # Validate forward primer
        fwd_valid, fwd_issues = self.validate_primer(pair.forward)
        issues.extend([f"Forward: {issue}" for issue in fwd_issues])

        # Validate reverse primer
        rev_valid, rev_issues = self.validate_primer(pair.reverse)
        issues.extend([f"Reverse: {issue}" for issue in rev_issues])

        # Check cross-dimer score
        if pair.cross_dimer_score > self.max_cross_dimer_score:
            issues.append(
                f"Cross-dimer score {pair.cross_dimer_score} exceeds maximum {self.max_cross_dimer_score}"
            )

        return (len(issues) == 0, issues)

    def _has_poly_x_run(self, sequence: str) -> bool:
        """
        Check for runs of identical bases.

        Poly-X runs can cause non-specific binding and secondary structure issues.

        Args:
            sequence: Sequence to check

        Returns:
            True if poly-X run detected, False otherwise
        """
        for base in "ACGT":
            repeat = base * (self.max_poly_x + 1)
            if repeat in sequence:
                return True
        return False

    def _has_gc_clamp(self, sequence: str) -> bool:
        """
        Check for G or C in last 5 bases (3' GC clamp).

        GC clamp improves primer stability and extension efficiency.

        Args:
            sequence: Sequence to check

        Returns:
            True if GC clamp present, False otherwise
        """
        last_5 = sequence[-5:]
        return "G" in last_5 or "C" in last_5

    def _has_3_prime_instability(self, sequence: str) -> bool:
        """
        Check for AT-rich 3' end (instability warning).

        AT-rich 3' ends can reduce primer binding specificity.

        Args:
            sequence: Sequence to check

        Returns:
            True if 3' end is AT-rich, False otherwise
        """
        last_3 = sequence[-3:]
        at_count = last_3.count("A") + last_3.count("T")
        return at_count >= 3
