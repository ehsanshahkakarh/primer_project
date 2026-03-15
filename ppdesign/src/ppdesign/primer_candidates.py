"""Generate primer candidates from conserved regions."""

import re
from typing import List, Optional
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

from .primer_types import Primer
from .conserved_finder import ConservedRegion


class PrimerCandidateGenerator:
    """Generate forward and reverse primer candidates from conserved regions."""

    def __init__(
        self,
        primer_min_length: int = 18,
        primer_max_length: int = 25,
        tm_min: float = 55.0,
        tm_max: float = 65.0,
        gc_min: float = 40.0,
        gc_max: float = 60.0,
        conservation_threshold: float = 0.8,
        max_degenerate_positions: int = 2,
    ):
        """
        Initialize primer candidate generator.

        Args:
            primer_min_length: Minimum primer length in bp
            primer_max_length: Maximum primer length in bp
            tm_min: Minimum melting temperature (°C)
            tm_max: Maximum melting temperature (°C)
            gc_min: Minimum GC content (%)
            gc_max: Maximum GC content (%)
            conservation_threshold: Minimum conservation fraction
            max_degenerate_positions: Maximum IUPAC degenerate positions allowed (0-5)

        Raises:
            ValueError: If constraints are invalid
        """
        if primer_min_length >= primer_max_length:
            raise ValueError("primer_min_length must be < primer_max_length")
        if not (0 <= gc_min <= 100 and 0 <= gc_max <= 100 and gc_min < gc_max):
            raise ValueError("Invalid GC range")
        if tm_min >= tm_max:
            raise ValueError("tm_min must be < tm_max")
        if not (0 <= max_degenerate_positions <= 5):
            raise ValueError("max_degenerate_positions must be between 0 and 5")

        self.primer_min_length = primer_min_length
        self.primer_max_length = primer_max_length
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.conservation_threshold = conservation_threshold
        self.max_degenerate_positions = max_degenerate_positions

    def generate_forward_primers(
        self, conserved_regions: List[ConservedRegion]
    ) -> List[Primer]:
        """
        Extract forward primers from conserved regions.

        Primers are extracted from individual sequences (not consensus) to handle
        IUPAC codes gracefully and generate more primer candidates.

        Args:
            conserved_regions: List of conserved regions to process

        Returns:
            List of valid forward primers
        """
        primers = []

        for region in conserved_regions:
            # Skip if conservation too low
            if region.conservation < self.conservation_threshold:
                continue

            # Use first target sequence's position as reference
            reference_seq_id = sorted(region.positions.keys())[0]
            reference_position = region.positions[reference_seq_id]

            # Extract candidates from EACH sequence, not consensus
            candidate_seqs = set()
            for seq in region.sequences:
                for length in range(self.primer_min_length, self.primer_max_length + 1):
                    if length <= len(seq):
                        candidate_seqs.add(seq[:length])

            # Validate each unique candidate
            for candidate_seq in candidate_seqs:
                primer = self._create_primer_if_valid(
                    sequence=candidate_seq,
                    position=reference_position,
                    strand="+",
                    conservation=region.conservation,
                    target_ids=list(region.positions.keys()),
                )

                if primer is not None:
                    primers.append(primer)

        return primers

    def generate_reverse_primers(
        self, conserved_regions: List[ConservedRegion]
    ) -> List[Primer]:
        """
        Extract reverse primers from conserved regions (stored as reverse-complement).

        Primers are extracted from individual sequences (not consensus) to handle
        IUPAC codes gracefully and generate more primer candidates.

        Args:
            conserved_regions: List of conserved regions to process

        Returns:
            List of valid reverse primers (sequences are reverse-complemented)
        """
        primers = []

        for region in conserved_regions:
            # Skip if conservation too low
            if region.conservation < self.conservation_threshold:
                continue

            # Use first target sequence's position as reference
            reference_seq_id = sorted(region.positions.keys())[0]
            reference_position = region.positions[reference_seq_id]

            # Extract candidates from EACH sequence, not consensus
            candidate_seqs = set()
            for seq in region.sequences:
                for length in range(self.primer_min_length, self.primer_max_length + 1):
                    if length <= len(seq):
                        # Extract from end of region and reverse complement
                        end_seq = seq[-length:]
                        rc_seq = str(Seq(end_seq).reverse_complement())
                        candidate_seqs.add(rc_seq)

            # Position is the binding site end on forward strand
            # Use length of first sequence as reference
            region_length = len(region.sequences[0]) if region.sequences else 0
            primer_position = reference_position + region_length

            # Validate each unique candidate
            for rc_seq in candidate_seqs:
                primer = self._create_primer_if_valid(
                    sequence=rc_seq,
                    position=primer_position,
                    strand="-",
                    conservation=region.conservation,
                    target_ids=list(region.positions.keys()),
                )

                if primer is not None:
                    primers.append(primer)

        return primers

    def _create_primer_if_valid(
        self,
        sequence: str,
        position: int,
        strand: str,
        conservation: float,
        target_ids: List[str],
    ) -> Optional[Primer]:
        """
        Create primer if it passes all validation checks.

        Returns:
            Primer object if valid, None otherwise
        """
        # Validate sequence composition
        if not self._is_valid_sequence(sequence):
            return None

        # Calculate properties
        tm = self._calculate_tm(sequence)
        gc = self._calculate_gc(sequence)

        # Check thermodynamic constraints
        if not (self.tm_min <= tm <= self.tm_max):
            return None
        if not (self.gc_min <= gc <= self.gc_max):
            return None

        # Create primer
        return Primer(
            sequence=sequence,
            position=position,
            strand=strand,
            tm=tm,
            gc_content=gc,
            hairpin_dg=self._check_hairpin(sequence),
            self_dimer_score=self._simple_self_dimer(sequence),
            conservation=conservation,
            target_sequences=target_ids,
        )

    def _is_valid_sequence(self, seq: str) -> bool:
        """
        Check if sequence meets composition requirements.

        Validates that sequence contains only valid IUPAC codes and has at most
        max_degenerate_positions non-ACGT bases.

        Args:
            seq: Sequence to check

        Returns:
            True if valid, False otherwise
        """
        # Check all bases are valid IUPAC codes
        if not re.match(r"^[ACGTRYMKSWBDHVN]+$", seq):
            return False

        # Count degenerate positions (not ACGT)
        degenerate_count = sum(1 for base in seq if base not in "ACGT")

        # Check against threshold
        return degenerate_count <= self.max_degenerate_positions

    def _calculate_tm(self, sequence: str) -> float:
        """
        Calculate melting temperature using BioPython (standard PCR conditions).

        Uses nearest-neighbor thermodynamics:
        - Na+ = 50 mM
        - Primer concentration = 250 nM

        Args:
            sequence: Primer sequence

        Returns:
            Melting temperature in °C
        """
        return mt.Tm_NN(sequence, Na=50, dnac1=250)

    def _calculate_gc(self, sequence: str) -> float:
        """
        Calculate GC percentage.

        Args:
            sequence: Primer sequence

        Returns:
            GC percentage (0-100)
        """
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100.0

    def _check_hairpin(self, sequence: str) -> Optional[float]:
        """
        Check hairpin formation using primer3-py if available.

        Args:
            sequence: Primer sequence

        Returns:
            Hairpin ΔG in kcal/mol if primer3 available, None otherwise
        """
        try:
            import primer3

            result = primer3.calc_hairpin(sequence)
            return result.dg / 1000.0  # Convert from cal/mol to kcal/mol
        except ImportError:
            return None

    def _simple_self_dimer(self, sequence: str) -> int:
        """
        Calculate simple self-complementarity score.

        Counts complementary base pairs when sequence is aligned with its reverse.

        Args:
            sequence: Primer sequence

        Returns:
            Number of complementary base pairs
        """
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        rev_seq = sequence[::-1]
        score = sum(
            1 for b1, b2 in zip(sequence, rev_seq) if complement.get(b1) == b2
        )
        return score
