#!/usr/bin/env python3
"""
CRISPR-Cas9 guide RNA finder for conserved regions.
Identifies guide RNAs with SpCas9 PAM sites (NGG/NAG) across multiple sequences.
"""

import logging
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass, field
from collections import defaultdict, Counter
import re
import numpy as np
from Bio.Seq import Seq


logger = logging.getLogger(__name__)


@dataclass
class GuideRNA:
    """Represents a CRISPR guide RNA candidate."""

    sequence: str  # 20bp guide sequence
    pam: str  # PAM sequence (NGG or NAG)
    strand: str  # '+' or '-'
    conservation: float  # Fraction of input sequences containing this guide (sequence-level coverage)
    positions: Dict[str, List[Tuple[int, str]]] = field(
        default_factory=dict
    )  # seq_id -> [(position, strand)]
    consensus: Optional[str] = None  # Consensus with degenerate bases if needed
    degenerate_count: int = 0  # Number of degenerate bases
    specificity_score: float = 0.0  # Score based on non-degenerate content
    quality_score: float = 0.0  # Overall quality score
    required_targets_hit: int = 0  # Number of required targets hit by this guide
    required_targets_total: int = 0  # Total number of required targets


class GuideRNAFinder:
    """Find conserved CRISPR guide RNAs across multiple sequences."""

    # IUPAC degenerate base codes
    IUPAC_CODES = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
    }

    def __init__(
        self,
        min_conservation: float = 0.8,
        max_degenerate: int = 3,
        pam_patterns: Optional[List[str]] = None,
        required_targets: Optional[List[str]] = None,
        no_degenerate: bool = False,
        min_guides: int = 1,
        max_guides: Optional[int] = None,
        min_additional_coverage: float = 0.05,
        perfect_coverage: bool = False,
        min_coverage_per_target: int = 3,
        max_total_grnas: int = 20,
        min_gc: Optional[float] = None,
        max_gc: Optional[float] = None,
        cluster_threshold: Optional[float] = None,
    ):
        """
        Initialize guide RNA finder.

        Args:
            min_conservation: Minimum fraction of input sequences that must contain the guide
                            (e.g., 0.8 = guide must be present in at least 80% of sequences)
            max_degenerate: Maximum number of degenerate bases allowed in consensus
            pam_patterns: List of PAM patterns to search for (default: ['NGG', 'NAG'])
            required_targets: List of sequence IDs that must be targeted by guide RNAs
            no_degenerate: If True, disable degenerate bases and use exact sequences only
            min_guides: Minimum number of guides to return (useful for getting multiple guides)
            max_guides: Maximum number of guides to return (None = no limit)
            min_additional_coverage: Minimum additional coverage threshold for selecting extra guides
            perfect_coverage: If True, enable perfect coverage mode with set cover optimization
            min_coverage_per_target: Minimum number of perfect gRNA matches required per target sequence
            max_total_grnas: Maximum total number of gRNAs to select in perfect coverage mode
            min_gc: Minimum GC content percentage for guide RNAs (0-100)
            max_gc: Maximum GC content percentage for guide RNAs (0-100)
            cluster_threshold: Similarity threshold for clustering sequences (0.0-1.0)
        """
        self.min_conservation = min_conservation
        self.max_degenerate = max_degenerate
        self.no_degenerate = no_degenerate
        self.min_guides = min_guides
        self.max_guides = max_guides
        self.min_additional_coverage = min_additional_coverage
        self.perfect_coverage = perfect_coverage
        self.min_coverage_per_target = min_coverage_per_target
        self.max_total_grnas = max_total_grnas
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.cluster_threshold = cluster_threshold
        self.pam_patterns = pam_patterns or ["NGG", "NAG"]
        self.guide_length = 20  # Standard SpCas9 guide length
        self.required_targets = set(required_targets) if required_targets else set()

    def _calculate_gc_content(self, seq: str) -> float:
        """Calculate GC content percentage of a sequence."""
        if not seq:
            return 0.0
        gc_count = seq.upper().count("G") + seq.upper().count("C")
        return (gc_count / len(seq)) * 100

    def _passes_gc_filter(self, seq: str) -> bool:
        """Check if sequence passes GC content filter."""
        if self.min_gc is None and self.max_gc is None:
            return True

        gc_content = self._calculate_gc_content(seq)

        if self.min_gc is not None and gc_content < self.min_gc:
            return False
        if self.max_gc is not None and gc_content > self.max_gc:
            return False

        return True

    def _cluster_sequences(self, sequences: Dict[str, str]) -> Dict[str, str]:
        """Cluster similar sequences based on similarity threshold."""
        if not self.cluster_threshold:
            return sequences

        from Bio import pairwise2

        clustered = {}
        processed = set()

        seq_ids = list(sequences.keys())
        for i, seq_id in enumerate(seq_ids):
            if seq_id in processed:
                continue

            # This sequence is a cluster representative
            clustered[seq_id] = sequences[seq_id]
            processed.add(seq_id)

            # Find similar sequences
            for j in range(i + 1, len(seq_ids)):
                other_id = seq_ids[j]
                if other_id in processed:
                    continue

                # Calculate similarity
                alignments = pairwise2.align.globalxx(
                    sequences[seq_id], sequences[other_id], one_alignment_only=True
                )
                if alignments:
                    score = alignments[0].score
                    max_len = max(len(sequences[seq_id]), len(sequences[other_id]))
                    similarity = score / max_len

                    if similarity >= self.cluster_threshold:
                        processed.add(other_id)

        logger.info(
            f"Clustered {len(sequences)} sequences to {len(clustered)} representatives"
        )
        return clustered

    def find_guide_rnas(
        self, sequences: Dict[str, str], allow_mismatches: int = 2
    ) -> List[GuideRNA]:
        """
        Find conserved guide RNAs across sequences.

        Args:
            sequences: Dictionary of sequence_id -> sequence
            allow_mismatches: Number of mismatches allowed when clustering guides

        Returns:
            List of GuideRNA objects sorted by quality score
        """
        if len(sequences) < 2:
            logger.warning("Need at least 2 sequences to find conserved guide RNAs")
            return []

        # Apply sequence clustering if threshold is specified
        if self.cluster_threshold is not None:
            sequences = self._cluster_sequences(sequences)
            logger.info(f"Clustered to {len(sequences)} representative sequences")

        # Perfect coverage mode uses different algorithm
        if self.perfect_coverage:
            return self._find_perfect_coverage_guides(sequences)

        # Step 1: Find all PAM sites and extract guides
        all_guides = self._extract_all_guides(sequences)
        logger.info(f"Found {len(all_guides)} total guide RNA candidates")

        # Step 2: Cluster similar guides
        guide_clusters = self._cluster_guides(all_guides, allow_mismatches)
        logger.info(f"Formed {len(guide_clusters)} guide clusters")

        # Step 3: Generate consensus and filter by conservation
        conserved_guides = []
        for cluster in guide_clusters:
            guide = self._process_guide_cluster(cluster, sequences)
            if guide and guide.conservation >= self.min_conservation:
                conserved_guides.append(guide)

        logger.info(f"Found {len(conserved_guides)} conserved guide RNAs")

        # Step 4: Calculate scores and sort
        for guide in conserved_guides:
            self._calculate_scores(guide)

        conserved_guides.sort(key=lambda x: x.quality_score, reverse=True)

        # Step 5: Apply greedy optimization in non-degenerate mode with required targets
        if self.no_degenerate and self.required_targets:
            conserved_guides = self._optimize_guide_selection(
                conserved_guides, sequences
            )
            logger.info(
                f"Optimized to {len(conserved_guides)} guides for maximum coverage"
            )

        return conserved_guides

    def find_pam_sites(self, sequence: str, pam_pattern: str) -> List[Tuple[int, str]]:
        """
        Find PAM sites in a sequence.

        Args:
            sequence: DNA sequence to search
            pam_pattern: PAM pattern (e.g., 'NGG', 'NAG')

        Returns:
            List of (position, strand) tuples
        """
        pam_sites = []

        # Convert PAM pattern to regex
        pattern = pam_pattern.replace("N", "[ATCG]")
        regex = re.compile(pattern)

        # Search forward strand
        for match in regex.finditer(sequence):
            pam_sites.append((match.start(), "+"))

        # Search reverse strand
        rev_comp = str(Seq(sequence).reverse_complement())
        for match in regex.finditer(rev_comp):
            # Convert position to forward strand coordinates
            pos = len(sequence) - match.end()
            pam_sites.append((pos, "-"))

        return pam_sites

    def extract_guide_sequence(
        self, sequence: str, pam_pos: int, strand: str, pam_length: int = 3
    ) -> Optional[str]:
        """
        Extract guide RNA sequence upstream of PAM.

        Args:
            sequence: DNA sequence
            pam_pos: PAM start position
            strand: '+' or '-'
            pam_length: Length of PAM sequence

        Returns:
            20bp guide sequence or None if too close to edge
        """
        if strand == "+":
            # Guide is upstream of PAM on forward strand
            guide_start = pam_pos - self.guide_length
            if guide_start < 0:
                return None
            guide_seq = sequence[guide_start:pam_pos]

        else:  # strand == '-'
            # Guide is downstream of PAM position on forward strand (upstream on reverse)
            guide_end = pam_pos + pam_length + self.guide_length
            if guide_end > len(sequence):
                return None
            guide_seq = sequence[pam_pos + pam_length : guide_end]
            # Reverse complement to get actual guide sequence
            guide_seq = str(Seq(guide_seq).reverse_complement())

        if len(guide_seq) != self.guide_length:
            return None

        return guide_seq.upper()

    def _extract_all_guides(self, sequences: Dict[str, str]) -> List[Dict]:
        """Extract all possible guide RNAs from all sequences."""
        all_guides = []

        for seq_id, sequence in sequences.items():
            sequence = sequence.upper()

            for pam_pattern in self.pam_patterns:
                pam_sites = self.find_pam_sites(sequence, pam_pattern)

                for pam_pos, strand in pam_sites:
                    guide_seq = self.extract_guide_sequence(
                        sequence, pam_pos, strand, len(pam_pattern)
                    )

                    if (
                        guide_seq
                        and "N" not in guide_seq
                        and self._passes_gc_filter(guide_seq)
                    ):  # Skip guides with N or wrong GC
                        all_guides.append(
                            {
                                "sequence": guide_seq,
                                "pam": pam_pattern,
                                "strand": strand,
                                "seq_id": seq_id,
                                "position": pam_pos,
                            }
                        )

        return all_guides

    def _cluster_guides(
        self, guides: List[Dict], max_mismatches: int
    ) -> List[List[Dict]]:
        """Cluster similar guide sequences allowing mismatches using vectorized Hamming distance."""
        if not guides:
            return []

        base_map = {"A": 0, "C": 1, "G": 2, "T": 3}

        encoded_guides = []
        for guide in guides:
            seq = guide["sequence"]
            arr = np.fromiter((base_map[base] for base in seq), dtype=np.uint8)
            encoded_guides.append(arr)

        clusters: List[List[Dict]] = []
        center_matrix = None

        for idx, guide in enumerate(guides):
            arr = encoded_guides[idx]
            assigned = False

            if center_matrix is not None and center_matrix.size:
                diffs = np.count_nonzero(center_matrix != arr, axis=1)
                match = np.where(diffs <= max_mismatches)[0]
                if match.size:
                    clusters[int(match[0])].append(guide)
                    assigned = True

            if not assigned:
                clusters.append([guide])
                if center_matrix is None or not center_matrix.size:
                    center_matrix = arr[np.newaxis, :]
                else:
                    center_matrix = np.vstack([center_matrix, arr])

        return clusters

    def _count_mismatches(self, seq1: str, seq2: str) -> int:
        """Count mismatches between two sequences."""
        if len(seq1) != len(seq2):
            return float("inf")
        return sum(1 for a, b in zip(seq1, seq2) if a != b)

    def _process_guide_cluster(
        self, cluster: List[Dict], all_sequences: Dict[str, str]
    ) -> Optional[GuideRNA]:
        """
        Process a cluster of similar guides to create a consensus GuideRNA.

        The conservation metric represents the fraction of input sequences that contain
        at least one guide RNA from this cluster. This measures sequence-level coverage
        rather than individual guide occurrence frequency.

        Args:
            cluster: List of guide dictionaries from the same cluster
            all_sequences: Dictionary of all input sequences {seq_id: sequence}

        Returns:
            GuideRNA object with calculated conservation metric, or None if invalid

        Conservation Calculation:
            conservation = number_of_sequences_with_guide / total_input_sequences

            Examples:
            - If 3 out of 4 sequences contain the guide: conservation = 0.75 (75%)
            - If all 5 sequences contain the guide: conservation = 1.0 (100%)
            - If only 1 out of 6 sequences contains the guide: conservation = 0.167 (16.7%)
        """
        if not cluster:
            return None

        # Track which sequences contain this guide cluster
        # Key insight: seq_coverage tracks unique sequences containing ANY guide from this cluster
        seq_coverage = defaultdict(list)
        for guide in cluster:
            seq_coverage[guide["seq_id"]].append((guide["position"], guide["strand"]))

        # Conservation = fraction of input sequences containing this guide
        # This counts sequence-level presence, not individual guide occurrences
        conservation = len(seq_coverage) / len(all_sequences)

        # In non-degenerate mode, select best individual sequence from cluster
        if self.no_degenerate:
            return self._select_best_guide_from_cluster(
                cluster, seq_coverage, conservation, all_sequences
            )

        # Generate consensus sequence
        consensus = self._generate_consensus([g["sequence"] for g in cluster])
        degenerate_count = sum(1 for base in consensus if base not in "ATCG")

        if degenerate_count > self.max_degenerate:
            return None

        # Use most common PAM
        pam_counts = Counter(g["pam"] for g in cluster)
        most_common_pam = pam_counts.most_common(1)[0][0]

        # Determine predominant strand
        strand_counts = Counter(g["strand"] for g in cluster)
        predominant_strand = strand_counts.most_common(1)[0][0]

        # Calculate required target coverage
        required_hit = 0
        if self.required_targets:
            required_hit = len(self.required_targets & set(seq_coverage.keys()))

        return GuideRNA(
            sequence=cluster[0]["sequence"],  # Representative sequence
            pam=most_common_pam,
            strand=predominant_strand,
            conservation=conservation,
            positions=dict(seq_coverage),
            consensus=consensus,
            degenerate_count=degenerate_count,
            required_targets_hit=required_hit,
            required_targets_total=len(self.required_targets),
        )

    def _select_best_guide_from_cluster(
        self,
        cluster: List[Dict],
        seq_coverage: Dict[str, List[Tuple[int, str]]],
        conservation: float,
        all_sequences: Dict[str, str],
    ) -> Optional[GuideRNA]:
        """
        Select the best individual guide from a cluster for non-degenerate mode.

        Prioritizes guides that:
        1. Cover the most required targets
        2. Have the highest individual coverage
        3. Appear most frequently in the cluster
        """
        best_guide = None
        best_score = -1

        # Count occurrences of each unique sequence in the cluster
        sequence_counts = Counter(g["sequence"] for g in cluster)

        for seq, count in sequence_counts.items():
            # Get all guides with this exact sequence
            guides_with_seq = [g for g in cluster if g["sequence"] == seq]

            # Track coverage for this specific sequence
            specific_coverage = defaultdict(list)
            for guide in guides_with_seq:
                specific_coverage[guide["seq_id"]].append(
                    (guide["position"], guide["strand"])
                )

            # Calculate metrics for this specific sequence
            seq_conservation = len(specific_coverage) / len(all_sequences)

            # Calculate required target coverage for this sequence
            required_hit = 0
            if self.required_targets:
                required_hit = len(
                    self.required_targets & set(specific_coverage.keys())
                )

            # Score: prioritize required targets, then coverage, then frequency
            score = (
                required_hit * 1000  # Heavily prioritize required targets
                + seq_conservation * 100  # Then conservation
                + count  # Then frequency in cluster
            )

            if score > best_score:
                best_score = score

                # Get most common PAM and strand for this sequence
                pam_counts = Counter(g["pam"] for g in guides_with_seq)
                most_common_pam = pam_counts.most_common(1)[0][0]

                strand_counts = Counter(g["strand"] for g in guides_with_seq)
                predominant_strand = strand_counts.most_common(1)[0][0]

                best_guide = GuideRNA(
                    sequence=seq,
                    pam=most_common_pam,
                    strand=predominant_strand,
                    conservation=seq_conservation,
                    positions=dict(specific_coverage),
                    consensus=seq,  # In non-degenerate mode, consensus = sequence
                    degenerate_count=0,
                    required_targets_hit=required_hit,
                    required_targets_total=len(self.required_targets),
                )

        return best_guide

    def _generate_consensus(self, sequences: List[str]) -> str:
        """Generate consensus sequence with IUPAC codes for degenerate positions."""
        if not sequences:
            return ""

        if len(sequences) == 1:
            return sequences[0]

        consensus = []
        for i in range(len(sequences[0])):
            bases = Counter(seq[i] for seq in sequences if i < len(seq))

            if len(bases) == 1:
                consensus.append(list(bases.keys())[0])
            else:
                # Find appropriate IUPAC code
                base_set = set(bases.keys())
                iupac_code = self._get_iupac_code(base_set)
                consensus.append(iupac_code)

        return "".join(consensus)

    def _get_iupac_code(self, bases: Set[str]) -> str:
        """Get IUPAC code for a set of bases."""
        if len(bases) == 1:
            return list(bases)[0]

        for code, base_list in self.IUPAC_CODES.items():
            if set(base_list) == bases:
                return code

        return "N"  # Any base

    def _optimize_guide_selection(
        self, guides: List[GuideRNA], sequences: Dict[str, str]
    ) -> List[GuideRNA]:
        """
        Optimize guide selection using greedy set cover algorithm.

        Prioritizes:
        1. 100% coverage of required targets
        2. Maximum coverage of additional sequences
        3. Minimum number of guides
        """
        optimized = []
        uncovered_required = self.required_targets.copy()
        uncovered_sequences = set(sequences.keys())

        # Phase 1: Cover all required targets
        while uncovered_required and guides:
            best_guide = None
            best_new_coverage = 0

            for guide in guides:
                # Count how many uncovered required targets this guide hits
                guide_targets = set(guide.positions.keys())
                new_required = uncovered_required & guide_targets

                if len(new_required) > best_new_coverage:
                    best_new_coverage = len(new_required)
                    best_guide = guide

            if best_guide and best_new_coverage > 0:
                optimized.append(best_guide)
                guides.remove(best_guide)

                # Update coverage
                guide_targets = set(best_guide.positions.keys())
                uncovered_required -= guide_targets
                uncovered_sequences -= guide_targets

                logger.debug(
                    f"Selected guide covering {best_new_coverage} required targets, "
                    f"{len(uncovered_required)} required remaining"
                )
            else:
                # No more guides can cover required targets
                if uncovered_required:
                    logger.warning(
                        f"Could not achieve 100% coverage of required targets. "
                        f"Missing: {uncovered_required}"
                    )
                break

        # Phase 2: Maximize coverage of remaining sequences
        # Continue until we've covered enough sequences or run out of useful guides
        # Use configurable threshold

        while (
            uncovered_sequences
            and guides
            and (self.max_guides is None or len(optimized) < self.max_guides)
        ):
            best_guide = None
            best_new_coverage = 0

            for guide in guides:
                # Count how many uncovered sequences this guide hits
                guide_targets = set(guide.positions.keys())
                new_sequences = uncovered_sequences & guide_targets

                if len(new_sequences) > best_new_coverage:
                    best_new_coverage = len(new_sequences)
                    best_guide = guide

            # Add guide if it meets threshold OR we need more guides to reach min_guides
            coverage_gain = (
                best_new_coverage / len(sequences) if best_new_coverage > 0 else 0
            )
            if best_guide and (
                coverage_gain >= self.min_additional_coverage
                or len(optimized) < self.min_guides
            ):
                optimized.append(best_guide)
                guides.remove(best_guide)

                # Update coverage
                guide_targets = set(best_guide.positions.keys())
                uncovered_sequences -= guide_targets

                logger.debug(
                    f"Selected guide covering {best_new_coverage} additional sequences, "
                    f"{len(uncovered_sequences)} uncovered remaining"
                )
            else:
                # No more guides provide sufficient additional coverage
                # But check if we need more guides to meet minimum
                if len(optimized) < self.min_guides and guides:
                    # Select best remaining guide by quality score
                    best_guide = max(guides, key=lambda g: g.quality_score)
                    optimized.append(best_guide)
                    guides.remove(best_guide)
                    logger.debug(
                        f"Added guide to meet minimum count (now {len(optimized)} guides)"
                    )
                else:
                    break

        # Log final coverage statistics
        total_covered = len(sequences) - len(uncovered_sequences)
        coverage_pct = (total_covered / len(sequences)) * 100
        logger.info(
            f"Optimization complete: {len(optimized)} guides selected, "
            f"covering {total_covered}/{len(sequences)} sequences ({coverage_pct:.1f}%)"
        )

        if not uncovered_required:
            logger.info("✓ All required targets covered")

        # Check if we found fewer guides than requested
        if len(optimized) < self.min_guides:
            logger.warning(
                f"Only {len(optimized)} guides found (requested {self.min_guides}). "
                f"Consider using degenerate mode (remove --no-degenerate) for more guides."
            )

        return optimized

    def _calculate_scores(self, guide: GuideRNA):
        """Calculate specificity and quality scores for a guide RNA."""
        # Specificity score based on non-degenerate content
        if guide.consensus:
            non_degen = sum(1 for base in guide.consensus if base in "ATCG")
            guide.specificity_score = non_degen / len(guide.consensus)
        else:
            guide.specificity_score = 1.0

        # Position-weighted specificity (5' end more important)
        if guide.consensus:
            position_weights = [
                1.5 if i < 10 else 1.0 for i in range(len(guide.consensus))
            ]
            weighted_score = sum(
                w if guide.consensus[i] in "ATCG" else 0
                for i, w in enumerate(position_weights)
            ) / sum(position_weights)
            guide.specificity_score = (guide.specificity_score + weighted_score) / 2

        # Calculate quality score with heavy weight on required targets
        if guide.required_targets_total > 0:
            # If required targets are specified, prioritize them heavily (60% weight)
            required_coverage = (
                guide.required_targets_hit / guide.required_targets_total
            )
            guide.quality_score = (
                required_coverage * 0.6  # Required targets get 60% weight
                + guide.conservation * 0.25  # General conservation 25%
                + guide.specificity_score * 0.15  # Specificity 15%
            )
        else:
            # Standard scoring when no required targets
            guide.quality_score = (
                guide.conservation * 0.6  # Weight conservation higher
                + guide.specificity_score * 0.4
            )

        # Penalize secondary PAM (NAG)
        if guide.pam == "NAG":
            guide.quality_score *= 0.9

    def format_guide_output(self, guide: GuideRNA) -> str:
        """Format guide RNA for output."""
        target_count = len(guide.positions)

        output = (
            f"Guide: {guide.consensus or guide.sequence}\n"
            f"PAM: {guide.pam} ({guide.strand} strand)\n"
            f"Conservation: {guide.conservation:.2%} ({target_count} sequences)\n"
        )

        # Add required target information if applicable
        if guide.required_targets_total > 0:
            output += f"Required targets: {guide.required_targets_hit}/{guide.required_targets_total}\n"

        output += (
            f"Degenerate bases: {guide.degenerate_count}\n"
            f"Specificity score: {guide.specificity_score:.3f}\n"
            f"Quality score: {guide.quality_score:.3f}\n"
        )

        return output

    # Perfect Coverage Mode Implementation

    def _find_perfect_coverage_guides(
        self, sequences: Dict[str, str]
    ) -> List[GuideRNA]:
        """
        Find guides using perfect coverage mode with greedy set cover algorithm.

        This mode ensures each target sequence has at least N perfect gRNA matches.
        """
        logger.info(
            f"Running perfect coverage mode: min {self.min_coverage_per_target} guides per target"
        )

        # Step 1: Find ALL possible gRNAs for each sequence
        all_guides_per_seq = self._find_all_grnas_per_sequence(sequences)

        total_guides = sum(len(guides) for guides in all_guides_per_seq.values())
        logger.info(f"Found {total_guides} total gRNA candidates across all sequences")

        # Step 2: Create unique guide set with exact sequence matching
        unique_guides = self._create_unique_guide_set(all_guides_per_seq, sequences)
        logger.info(f"Created {len(unique_guides)} unique gRNA candidates")

        # Step 3: Apply greedy set cover algorithm
        selected_guides = self._greedy_set_cover_perfect(unique_guides, sequences)

        # Step 4: Calculate scores and sort
        for guide in selected_guides:
            self._calculate_scores(guide)

        selected_guides.sort(key=lambda x: x.quality_score, reverse=True)

        # Step 5: Generate coverage statistics
        coverage_stats = self._calculate_perfect_coverage_stats(
            selected_guides, sequences
        )
        coverage_stats['min_coverage_per_target'] = self.min_coverage_per_target
        coverage_stats['max_total_grnas'] = self.max_total_grnas
        logger.info(f"Perfect coverage achieved: {coverage_stats}")

        return selected_guides

    def _find_all_grnas_per_sequence(
        self, sequences: Dict[str, str]
    ) -> Dict[str, List[Dict]]:
        """Find ALL possible gRNAs in each sequence (not just conserved ones)."""
        all_guides_per_seq = {}

        for seq_id, sequence in sequences.items():
            sequence = sequence.upper()
            seq_guides = []

            for pam_pattern in self.pam_patterns:
                pam_sites = self.find_pam_sites(sequence, pam_pattern)

                for pam_pos, strand in pam_sites:
                    guide_seq = self.extract_guide_sequence(
                        sequence, pam_pos, strand, len(pam_pattern)
                    )

                    if guide_seq and "N" not in guide_seq:
                        seq_guides.append(
                            {
                                "sequence": guide_seq,
                                "pam": pam_pattern,
                                "strand": strand,
                                "seq_id": seq_id,
                                "position": pam_pos,
                            }
                        )

            all_guides_per_seq[seq_id] = seq_guides
            logger.debug(f"Found {len(seq_guides)} gRNAs in sequence {seq_id}")

        return all_guides_per_seq

    def _create_unique_guide_set(
        self, all_guides_per_seq: Dict[str, List[Dict]], sequences: Dict[str, str]
    ) -> List[GuideRNA]:
        """Create unique guide set from all possible guides."""
        # Group by exact sequence match
        guide_groups = defaultdict(list)

        for seq_id, guides in all_guides_per_seq.items():
            for guide in guides:
                key = (guide["sequence"], guide["pam"], guide["strand"])
                guide_groups[key].append(guide)

        # Convert to GuideRNA objects
        unique_guides = []
        for (seq, pam, strand), guide_list in guide_groups.items():
            # Track positions for this guide across sequences
            positions = defaultdict(list)
            for guide in guide_list:
                positions[guide["seq_id"]].append((guide["position"], guide["strand"]))

            # Calculate conservation (fraction of sequences containing this guide)
            conservation = len(positions) / len(sequences)

            # Calculate required target coverage
            required_hit = 0
            if self.required_targets:
                required_hit = len(self.required_targets & set(positions.keys()))

            guide_rna = GuideRNA(
                sequence=seq,
                pam=pam,
                strand=strand,
                conservation=conservation,
                positions=dict(positions),
                consensus=seq,  # Perfect coverage uses exact sequences
                degenerate_count=0,
                required_targets_hit=required_hit,
                required_targets_total=len(self.required_targets),
            )

            unique_guides.append(guide_rna)

        return unique_guides

    def _greedy_set_cover_perfect(
        self, guides: List[GuideRNA], sequences: Dict[str, str]
    ) -> List[GuideRNA]:
        """
        Greedy set cover algorithm ensuring each sequence has at least N perfect gRNA matches.
        """
        selected = []
        target_coverage = {seq_id: 0 for seq_id in sequences.keys()}
        remaining_guides = guides.copy()

        logger.info(f"Starting greedy set cover for {len(sequences)} targets")
        logger.info(
            f"Target: {self.min_coverage_per_target} guides per target, max {self.max_total_grnas} total"
        )

        iteration = 0
        while remaining_guides and len(selected) < self.max_total_grnas:
            iteration += 1

            # Find sequences that still need more coverage
            undercover_targets = {
                seq_id: count
                for seq_id, count in target_coverage.items()
                if count < self.min_coverage_per_target
            }

            if not undercover_targets:
                logger.info("All targets have sufficient coverage")
                break

            # Score guides by how much they help undercover targets
            best_guide = None
            best_score = -1
            best_new_coverage = 0

            for guide in remaining_guides:
                # Calculate benefit for undercover targets
                undercover_benefit = 0
                total_new_coverage = 0

                for seq_id in guide.positions.keys():
                    if seq_id in undercover_targets:
                        # Higher weight for targets that need more coverage
                        need = self.min_coverage_per_target - target_coverage[seq_id]
                        undercover_benefit += min(1, need)  # Cap at 1 per target

                    if target_coverage[seq_id] < self.min_coverage_per_target:
                        total_new_coverage += 1

                # Score: prioritize undercover targets, then quality
                score = (
                    undercover_benefit * 1000  # Heavily weight undercover targets
                    + total_new_coverage * 100  # Then new coverage
                    + guide.conservation * 10  # Then conservation
                )

                if score > best_score:
                    best_score = score
                    best_guide = guide
                    best_new_coverage = total_new_coverage

            if best_guide:
                selected.append(best_guide)
                remaining_guides.remove(best_guide)

                # Update coverage counts
                for seq_id in best_guide.positions.keys():
                    target_coverage[seq_id] += 1

                undercover_count = len(undercover_targets)
                logger.debug(
                    f"Iteration {iteration}: Selected guide covering {best_new_coverage} targets, "
                    f"{undercover_count} targets still need coverage"
                )
            else:
                logger.warning("No more beneficial guides found")
                break

        # Final statistics
        final_undercover = sum(
            1
            for count in target_coverage.values()
            if count < self.min_coverage_per_target
        )
        coverage_achieved = len(sequences) - final_undercover

        logger.info(f"Selected {len(selected)} guides")
        logger.info(
            f"Coverage achieved: {coverage_achieved}/{len(sequences)} targets "
            f"({coverage_achieved/len(sequences):.1%})"
        )

        if final_undercover > 0:
            undercover_seqs = [
                seq_id
                for seq_id, count in target_coverage.items()
                if count < self.min_coverage_per_target
            ]
            logger.warning(
                f"{final_undercover} targets with insufficient coverage: {undercover_seqs[:5]}..."
            )

        return selected

    def _calculate_perfect_coverage_stats(
        self, guides: List[GuideRNA], sequences: Dict[str, str]
    ) -> Dict:
        """Calculate detailed coverage statistics for perfect coverage mode."""
        # Create coverage matrix
        coverage_matrix = self._generate_coverage_matrix(guides, sequences)

        # Calculate per-target statistics
        target_stats = {}
        for seq_id in sequences.keys():
            covering_guides = [
                i for i, guide in enumerate(guides) if seq_id in guide.positions
            ]
            target_stats[seq_id] = {
                "guide_count": len(covering_guides),
                "guide_indices": covering_guides,
                "sufficient_coverage": len(covering_guides)
                >= self.min_coverage_per_target,
            }

        # Overall statistics
        total_targets = len(sequences)
        covered_targets = sum(
            1 for stats in target_stats.values() if stats["sufficient_coverage"]
        )
        avg_coverage = (
            sum(stats["guide_count"] for stats in target_stats.values()) / total_targets
        )

        stats = {
            "total_guides": len(guides),
            "total_targets": total_targets,
            "covered_targets": covered_targets,
            "coverage_percentage": (
                covered_targets / total_targets if total_targets > 0 else 0
            ),
            "average_guides_per_target": avg_coverage,
            "min_guides_per_target": (
                min(stats["guide_count"] for stats in target_stats.values())
                if target_stats
                else 0
            ),
            "max_guides_per_target": (
                max(stats["guide_count"] for stats in target_stats.values())
                if target_stats
                else 0
            ),
            "target_stats": target_stats,
            "coverage_matrix": coverage_matrix,
        }

        return stats

    def _generate_coverage_matrix(
        self, guides: List[GuideRNA], sequences: Dict[str, str]
    ) -> List[List[int]]:
        """Generate binary coverage matrix: guides x targets."""
        matrix = []
        seq_ids = list(sequences.keys())

        for guide in guides:
            row = []
            for seq_id in seq_ids:
                # 1 if guide covers this target, 0 otherwise
                row.append(1 if seq_id in guide.positions else 0)
            matrix.append(row)

        return matrix
