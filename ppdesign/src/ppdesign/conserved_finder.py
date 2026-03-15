"""Utilities for finding conserved regions in sequences."""

import logging
import math
import shutil
from typing import List, Dict, Optional
from collections import defaultdict, Counter
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from dataclasses import dataclass
import subprocess
import tempfile
from pathlib import Path
import hashlib

logger = logging.getLogger(__name__)


@dataclass
class ConservedRegion:
    """Represents a conserved region found across sequences."""

    start: int
    end: int
    sequences: List[str]
    conservation: float
    consensus: str
    positions: Dict[str, int]  # sequence_id -> position in that sequence


class KmerBasedFinder:
    """Find conserved regions using k-mer approach."""

    def __init__(self, kmer_size: int = 20, min_conservation: float = 0.8):
        """
        Initialize k-mer based conserved region finder.

        Args:
            kmer_size: Size of k-mers to use for initial search
            min_conservation: Minimum fraction of sequences that must contain the region
        """
        self.kmer_size = kmer_size
        self.min_conservation = min_conservation
        self._alignment_cache = {}

    def find_conserved_regions(
        self,
        sequences: Dict[str, str],
        min_length: int = 15,
        max_length: int = 30,
        allow_mismatches: int = 2,
        align_regions: bool = False,
    ) -> List[ConservedRegion]:
        """
        Find conserved regions across sequences using k-mer seeding and extension.

        Args:
            sequences: Dictionary of sequence_id -> sequence
            min_length: Minimum length of conserved region
            max_length: Maximum length of conserved region
            allow_mismatches: Number of mismatches allowed in conserved region
            align_regions: If True, use MAFFT to refine regions with nucleotide-level alignment

        Returns:
            List of conserved regions
        """
        if len(sequences) < 2:
            logger.warning("Need at least 2 sequences to find conserved regions")
            return []

        # Step 1: Build k-mer index
        kmer_index = self._build_kmer_index(sequences)

        # Step 2: Find shared k-mers
        shared_kmers = self._find_shared_kmers(kmer_index, len(sequences))

        # Step 3: Extend k-mers to find longer conserved regions
        conserved_regions = []
        processed_positions = set()

        for kmer, positions in shared_kmers.items():
            # Check if we've already processed this region
            skip = False
            for seq_id, pos in positions.items():
                if (seq_id, pos) in processed_positions:
                    skip = True
                    break

            if skip:
                continue

            # Try to extend this k-mer
            region = self._extend_kmer(
                kmer, positions, sequences, min_length, max_length, allow_mismatches
            )

            if region:
                # Optionally refine with nucleotide-level alignment
                if align_regions:
                    region = self._refine_with_alignment(
                        region, sequences, min_length, max_length
                    )

                if region:  # May be None after alignment refinement
                    conserved_regions.append(region)
                    # Mark positions as processed
                    for seq_id, pos in region.positions.items():
                        for i in range(region.start, region.end):
                            processed_positions.add((seq_id, pos + i - region.start))

        # Sort by conservation and length
        conserved_regions.sort(
            key=lambda r: (r.conservation, r.end - r.start), reverse=True
        )

        return conserved_regions

    def _build_kmer_index(
        self, sequences: Dict[str, str]
    ) -> Dict[str, Dict[str, List[int]]]:
        """Build index of k-mer positions in each sequence."""
        kmer_index = defaultdict(lambda: defaultdict(list))

        for seq_id, sequence in sequences.items():
            seq_upper = sequence.upper()
            for i in range(len(seq_upper) - self.kmer_size + 1):
                kmer = seq_upper[i : i + self.kmer_size]
                if "N" not in kmer:  # Skip k-mers with ambiguous bases
                    kmer_index[kmer][seq_id].append(i)

        return kmer_index

    def _find_shared_kmers(
        self, kmer_index: Dict[str, Dict[str, List[int]]], total_sequences: int
    ) -> Dict[str, Dict[str, int]]:
        """Find k-mers present in sufficient fraction of sequences."""
        shared_kmers = {}
        min_sequences = max(1, math.ceil(total_sequences * self.min_conservation))

        for kmer, seq_positions in kmer_index.items():
            if len(seq_positions) >= min_sequences:
                # Take first position in each sequence
                positions = {
                    seq_id: pos_list[0] for seq_id, pos_list in seq_positions.items()
                }
                shared_kmers[kmer] = positions

        logger.info(
            f"Found {len(shared_kmers)} shared k-mers in at least {min_sequences} sequences"
        )
        return shared_kmers

    def _extend_kmer(
        self,
        kmer: str,
        positions: Dict[str, int],
        sequences: Dict[str, str],
        min_length: int,
        max_length: int,
        allow_mismatches: int,
    ) -> Optional[ConservedRegion]:
        """Extend k-mers and enforce mismatch limits."""
        left_extension = 0
        right_extension = self.kmer_size

        # Grow to the left
        while True:
            bases = []
            for seq_id, pos in positions.items():
                new_pos = pos - left_extension - 1
                if new_pos < 0:
                    bases = []
                    break
                bases.append(sequences[seq_id][new_pos].upper())
            if not bases:
                break
            base_counts = Counter(bases)
            if base_counts.most_common(1)[0][1] >= len(positions) * self.min_conservation:
                left_extension += 1
            else:
                break

        # Grow to the right
        while right_extension < max_length:
            bases = []
            for seq_id, pos in positions.items():
                new_pos = pos + right_extension
                if new_pos >= len(sequences[seq_id]):
                    bases = []
                    break
                bases.append(sequences[seq_id][new_pos].upper())
            if not bases:
                break
            base_counts = Counter(bases)
            if base_counts.most_common(1)[0][1] >= len(positions) * self.min_conservation:
                right_extension += 1
            else:
                break

        total_length = left_extension + right_extension
        if total_length < min_length:
            return None

        extended = {}
        final_positions: Dict[str, int] = {}
        for seq_id, pos in positions.items():
            start = pos - left_extension
            end = pos + right_extension
            extended[seq_id] = sequences[seq_id][start:end].upper()
            final_positions[seq_id] = start

        required_sequences = max(1, math.ceil(self.min_conservation * len(sequences)))

        def compute_mismatches(window_dict: Dict[str, str]) -> Dict[str, int]:
            mismatch_counts = {sid: 0 for sid in window_dict}
            for idx in range(total_length):
                col_bases = Counter(window[idx] for window in window_dict.values())
                consensus_base = col_bases.most_common(1)[0][0]
                for sid, window in window_dict.items():
                    if window[idx] != consensus_base:
                        mismatch_counts[sid] += 1
            return mismatch_counts

        filtered = dict(extended)
        if allow_mismatches >= 0:
            while True:
                mismatch_counts = compute_mismatches(filtered)
                to_remove = [sid for sid, count in mismatch_counts.items() if count > allow_mismatches]
                if not to_remove:
                    break
                for sid in to_remove:
                    filtered.pop(sid, None)
                    final_positions.pop(sid, None)
                if len(filtered) < required_sequences:
                    return None

        if not filtered:
            return None

        final_sequences = list(filtered.values())
        consensus = self._create_consensus(final_sequences)
        conservation = len(filtered) / len(sequences)

        return ConservedRegion(
            start=0,
            end=total_length,
            sequences=final_sequences,
            conservation=conservation,
            consensus=consensus,
            positions={sid: final_positions[sid] for sid in filtered},
        )

    def _create_consensus(self, sequences: List[str]) -> str:
        """Create consensus sequence with IUPAC ambiguity codes."""
        consensus = []
        seq_length = len(sequences[0])

        iupac_codes = {
            frozenset(["A"]): "A",
            frozenset(["C"]): "C",
            frozenset(["G"]): "G",
            frozenset(["T"]): "T",
            frozenset(["A", "G"]): "R",
            frozenset(["C", "T"]): "Y",
            frozenset(["G", "C"]): "S",
            frozenset(["A", "T"]): "W",
            frozenset(["G", "T"]): "K",
            frozenset(["A", "C"]): "M",
            frozenset(["C", "G", "T"]): "B",
            frozenset(["A", "G", "T"]): "D",
            frozenset(["A", "C", "T"]): "H",
            frozenset(["A", "C", "G"]): "V",
            frozenset(["A", "C", "G", "T"]): "N",
        }

        for i in range(seq_length):
            bases = set(seq[i] for seq in sequences if seq[i] != "N")
            consensus.append(iupac_codes.get(frozenset(bases), "N"))

        return "".join(consensus)

    def _refine_with_alignment(
        self,
        region: ConservedRegion,
        sequences: Dict[str, str],
        min_length: int,
        max_length: int,
    ) -> Optional[ConservedRegion]:
        """
        Refine a conserved region using MAFFT alignment.

        Args:
            region: The conserved region to refine
            sequences: Full sequences dict
            min_length: Minimum region length
            max_length: Maximum region length

        Returns:
            Refined ConservedRegion or None if refinement fails
        """
        # Check if MAFFT is available
        if shutil.which("mafft") is None:
            logger.warning(
                "MAFFT not found, skipping alignment refinement. "
                "Install MAFFT for highest quality primer design."
            )
            return region

        # Skip alignment for perfect matches (no IUPAC codes)
        if all(base in "ACGT" for base in region.consensus):
            return region

        try:
            # Create alignment
            alignment, conserved_cols = self._align_region_sequences(region.sequences)

            if not conserved_cols:
                logger.debug(
                    f"No conserved columns after alignment for region at {region.start}"
                )
                return None

            # Extract ungapped, conserved segments
            refined_sequences = self._extract_conserved_segments(
                alignment, conserved_cols, min_length, max_length
            )

            if not refined_sequences or len(refined_sequences[0]) < min_length:
                return None

            # Update region with refined sequences
            consensus = self._create_consensus(refined_sequences)
            conservation = len(refined_sequences) / len(sequences)

            return ConservedRegion(
                start=region.start,
                end=region.start + len(refined_sequences[0]),
                sequences=refined_sequences,
                conservation=conservation,
                consensus=consensus,
                positions=region.positions,
            )

        except Exception as e:
            logger.debug(f"Alignment refinement failed: {e}, using original region")
            return region

    def _align_region_sequences(self, sequences: List[str]) -> tuple:
        """
        Align sequences using MAFFT and identify conserved columns.

        Args:
            sequences: List of sequences to align

        Returns:
            Tuple of (alignment object, list of conserved column indices)
        """
        # Write sequences to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp:
            for i, seq in enumerate(sequences):
                record = SeqRecord(Seq(seq), id=f"seq{i}", description="")
                SeqIO.write([record], tmp, "fasta")
            tmp_path = tmp.name

        aligned_path = tmp_path + ".aln"

        try:
            # Run MAFFT
            cmd = ["mafft", "--auto", "--quiet", tmp_path]
            with open(aligned_path, "w") as out:
                result = subprocess.run(
                    cmd, stdout=out, stderr=subprocess.PIPE, text=True
                )
                if result.returncode != 0:
                    raise RuntimeError(f"MAFFT failed: {result.stderr}")

            # Parse alignment
            alignment = AlignIO.read(aligned_path, "fasta")

            # Calculate column-wise conservation
            conserved_cols = self._identify_conserved_columns(
                alignment, self.min_conservation
            )

            return alignment, conserved_cols

        finally:
            Path(tmp_path).unlink(missing_ok=True)
            Path(aligned_path).unlink(missing_ok=True)

    def _identify_conserved_columns(self, alignment, threshold: float) -> List[int]:
        """
        Identify columns in alignment that meet conservation threshold.

        Args:
            alignment: Bio.Align alignment object
            threshold: Minimum conservation fraction

        Returns:
            List of column indices that are conserved
        """
        conserved_indices = []
        alignment_length = alignment.get_alignment_length()
        num_sequences = len(alignment)

        for col_idx in range(alignment_length):
            col = alignment[:, col_idx]

            # Count non-gap bases
            bases = [base for base in col if base != "-"]

            if not bases:
                continue

            # Calculate most common base frequency
            base_counts = Counter(bases)
            most_common_count = base_counts.most_common(1)[0][1]

            # Check if conservation meets threshold
            conservation = most_common_count / num_sequences

            if conservation >= threshold:
                conserved_indices.append(col_idx)

        return conserved_indices

    def _extract_conserved_segments(
        self,
        alignment,
        conserved_cols: List[int],
        min_length: int,
        max_length: int,
    ) -> List[str]:
        """
        Extract ungapped conserved segments from alignment.

        Args:
            alignment: Bio.Align alignment object
            conserved_cols: List of conserved column indices
            min_length: Minimum segment length
            max_length: Maximum segment length

        Returns:
            List of ungapped sequences from conserved region
        """
        if not conserved_cols:
            return []

        # Find longest contiguous stretch of conserved columns
        best_start = 0
        best_end = 0
        current_start = conserved_cols[0]
        current_end = conserved_cols[0]

        for i in range(1, len(conserved_cols)):
            if conserved_cols[i] == conserved_cols[i - 1] + 1:
                current_end = conserved_cols[i]
            else:
                if current_end - current_start > best_end - best_start:
                    best_start = current_start
                    best_end = current_end
                current_start = conserved_cols[i]
                current_end = conserved_cols[i]

        # Check final stretch
        if current_end - current_start > best_end - best_start:
            best_start = current_start
            best_end = current_end

        # Extract sequences without gaps
        segment_length = best_end - best_start + 1

        if segment_length < min_length or segment_length > max_length:
            return []

        ungapped_seqs = []
        for record in alignment:
            segment = str(record.seq[best_start : best_end + 1])
            ungapped = segment.replace("-", "")

            # Only include if not too many gaps
            if len(ungapped) >= min_length * 0.8:  # Allow up to 20% gaps
                ungapped_seqs.append(ungapped)

        return ungapped_seqs


class MSABasedFinder:
    """Find conserved regions using multiple sequence alignment."""

    def __init__(self, min_conservation: float = 0.8):
        """
        Initialize MSA-based conserved region finder.

        Args:
            min_conservation: Minimum fraction of sequences that must be conserved
        """
        self.min_conservation = min_conservation
        self._alignment_cache = {}

    def find_conserved_regions(
        self,
        sequences: Dict[str, str],
        min_length: int = 15,
        max_length: int = 30,
        window_size: int = None,
        step_size: int = None,
    ) -> List[ConservedRegion]:
        """
        Find conserved regions using multiple sequence alignment.

        Args:
            sequences: Dictionary of sequence_id -> sequence
            min_length: Minimum length of conserved region
            max_length: Maximum length of conserved region
            window_size: Size of sliding window (default: max_length)
            step_size: Step size for sliding window (default: 1)

        Returns:
            List of conserved regions
        """
        if len(sequences) < 2:
            logger.warning("Need at least 2 sequences to find conserved regions")
            return []

        # Set defaults
        if window_size is None:
            window_size = max_length
        if step_size is None:
            step_size = 1

        # Create temporary FASTA file
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp:
            for seq_id, seq in sequences.items():
                record = SeqRecord(Seq(seq), id=seq_id, description="")
                SeqIO.write([record], tmp, "fasta")
            tmp_path = tmp.name

        aligned_path = tmp_path + ".aln"
        cache_key = hashlib.sha1(
            "||".join(f"{sid}:{seq}" for sid, seq in sorted(sequences.items())).encode()
        ).hexdigest() + f"|{min_length}|{max_length}|{window_size}|{step_size}"

        try:
            if cache_key in self._alignment_cache:
                alignment = self._alignment_cache[cache_key]
            else:
                if shutil.which("mafft") is None:
                    raise RuntimeError('MAFFT executable not found on PATH.')

                cmd = ["mafft", "--auto", "--quiet", tmp_path]

                with open(aligned_path, "w") as out:
                    result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
                    if result.returncode != 0:
                        logger.error(f"MAFFT failed: {result.stderr.decode()}")
                        return []

                alignment = AlignIO.read(aligned_path, "fasta")
                self._alignment_cache[cache_key] = alignment

            conserved_regions = self._find_conserved_windows(
                alignment, min_length, max_length, window_size, step_size
            )

            return conserved_regions

        finally:
            Path(tmp_path).unlink(missing_ok=True)
            Path(aligned_path).unlink(missing_ok=True)

    def _find_conserved_windows(
        self,
        alignment,
        min_length: int,
        max_length: int,
        window_size: int,
        step_size: int,
    ) -> List[ConservedRegion]:
        """Find conserved regions in alignment using sliding window."""
        conserved_regions = []
        alignment_length = alignment.get_alignment_length()
        num_sequences = len(alignment)
        min_conserved = int(num_sequences * self.min_conservation)

        # Slide window across alignment
        for start in range(0, alignment_length - window_size + 1, step_size):
            end = min(start + window_size, alignment_length)

            # Get window sequences (without gaps)
            window_seqs = {}
            positions = {}

            for record in alignment:
                # Get subsequence and remove gaps
                subseq = str(record.seq[start:end]).replace("-", "")

                if len(subseq) >= min_length:
                    window_seqs[record.id] = subseq
                    # Calculate actual position in original sequence
                    gaps_before = str(record.seq[:start]).count("-")
                    positions[record.id] = start - gaps_before

            # Check if enough sequences have content in this window
            if len(window_seqs) >= min_conserved:
                # Check conservation
                if self._is_conserved(list(window_seqs.values())):
                    consensus = self._create_consensus(list(window_seqs.values()))

                    region = ConservedRegion(
                        start=0,
                        end=len(consensus),
                        sequences=list(window_seqs.values()),
                        conservation=len(window_seqs) / num_sequences,
                        consensus=consensus,
                        positions=positions,
                    )

                    conserved_regions.append(region)

        # Remove redundant regions
        conserved_regions = self._remove_redundant_regions(conserved_regions)

        return conserved_regions

    def _is_conserved(self, sequences: List[str]) -> bool:
        """Check if sequences are sufficiently conserved."""
        if not sequences:
            return False

        # Use first sequence as reference
        ref_seq = sequences[0]
        conserved_count = 1

        for seq in sequences[1:]:
            # Simple similarity check (could be improved)
            if self._sequence_similarity(ref_seq, seq) >= 0.8:
                conserved_count += 1

        return conserved_count >= len(sequences) * self.min_conservation

    def _sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate simple sequence similarity."""
        if len(seq1) != len(seq2):
            # Handle different lengths
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]

        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1) if seq1 else 0

    def _create_consensus(self, sequences: List[str]) -> str:
        """Create consensus sequence with IUPAC ambiguity codes."""
        if not sequences:
            return ""

        # Find minimum length
        min_len = min(len(seq) for seq in sequences)

        consensus = []
        iupac_codes = {
            frozenset(["A"]): "A",
            frozenset(["C"]): "C",
            frozenset(["G"]): "G",
            frozenset(["T"]): "T",
            frozenset(["A", "G"]): "R",
            frozenset(["C", "T"]): "Y",
            frozenset(["G", "C"]): "S",
            frozenset(["A", "T"]): "W",
            frozenset(["G", "T"]): "K",
            frozenset(["A", "C"]): "M",
            frozenset(["C", "G", "T"]): "B",
            frozenset(["A", "G", "T"]): "D",
            frozenset(["A", "C", "T"]): "H",
            frozenset(["A", "C", "G"]): "V",
            frozenset(["A", "C", "G", "T"]): "N",
        }

        for i in range(min_len):
            bases = set(
                seq[i].upper() for seq in sequences if i < len(seq) and seq[i] != "N"
            )
            consensus.append(iupac_codes.get(frozenset(bases), "N"))

        return "".join(consensus)

    def _remove_redundant_regions(
        self, regions: List[ConservedRegion]
    ) -> List[ConservedRegion]:
        """Remove overlapping or redundant regions."""
        if not regions:
            return []

        # Sort by conservation and length
        regions.sort(key=lambda r: (r.conservation, r.end - r.start), reverse=True)

        kept_regions = []

        for region in regions:
            # Check if this region overlaps with any kept region
            is_redundant = False

            for kept in kept_regions:
                # Check overlap in consensus
                if (
                    region.consensus in kept.consensus
                    or kept.consensus in region.consensus
                ):
                    is_redundant = True
                    break

            if not is_redundant:
                kept_regions.append(region)

        return kept_regions


class MSAConservedRegionScanner:
    """Find conserved primer sites by scanning a MAFFT multiple sequence alignment.

    Unlike MSABasedFinder which uses pairwise similarity windows, this class
    computes per-column conservation and finds all contiguous conserved stretches
    -- the natural approach for primer site discovery from homologous sequences.
    """

    MAFFT_COMMANDS = {
        "linsi": ["mafft", "--localpair", "--maxiterate", "1000", "--quiet"],
        "qinsi": ["mafft-qinsi", "--quiet"],
        "auto": ["mafft", "--auto", "--quiet"],
    }

    def __init__(self, min_conservation: float = 0.8, mafft_strategy: str = "linsi"):
        if mafft_strategy not in self.MAFFT_COMMANDS:
            raise ValueError(
                f"Unknown MAFFT strategy '{mafft_strategy}', "
                f"choose from: {list(self.MAFFT_COMMANDS)}"
            )
        self.min_conservation = min_conservation
        self.mafft_strategy = mafft_strategy

    def find_conserved_regions(
        self,
        sequences: Dict[str, str],
        min_length: int = 18,
        max_length: int = 25,
    ) -> List[ConservedRegion]:
        """Find conserved regions suitable as primer binding sites.

        Args:
            sequences: Dict of sequence_id -> sequence
            min_length: Minimum region length (bp)
            max_length: Maximum region length (bp)

        Returns:
            List of primer-length ConservedRegion objects
        """
        if len(sequences) < 2:
            logger.warning("Need at least 2 sequences to find conserved regions")
            return []

        executable = self.MAFFT_COMMANDS[self.mafft_strategy][0]
        if shutil.which(executable) is None:
            raise RuntimeError(f"'{executable}' not found on PATH.")

        alignment = self._run_mafft(sequences)
        col_cons = self._column_conservation(alignment)
        runs = self._find_conserved_runs(col_cons, min_length)

        logger.info(
            f"Found {len(runs)} conserved runs (>= {min_length} columns) "
            f"in {alignment.get_alignment_length()}-column alignment"
        )

        regions = []
        seq_ids = [record.id for record in alignment]
        step = max(1, min_length // 2)

        for run_start, run_end in runs:
            for win_start in range(run_start, run_end - min_length + 1, step):
                win_end = min(win_start + max_length, run_end)
                if win_end - win_start < min_length:
                    continue
                region = self._extract_region(
                    alignment, seq_ids, len(sequences), win_start, win_end
                )
                if region is not None:
                    regions.append(region)

        regions = self._remove_redundant(regions)
        logger.info(f"Produced {len(regions)} primer-site regions after deduplication")
        return regions

    def _run_mafft(self, sequences: Dict[str, str]):
        """Run MAFFT and return parsed alignment."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp:
            for seq_id, seq in sequences.items():
                record = SeqRecord(Seq(seq), id=seq_id, description="")
                SeqIO.write([record], tmp, "fasta")
            tmp_path = tmp.name

        aligned_path = tmp_path + ".aln"
        try:
            cmd = list(self.MAFFT_COMMANDS[self.mafft_strategy]) + [tmp_path]
            with open(aligned_path, "w") as out:
                result = subprocess.run(
                    cmd, stdout=out, stderr=subprocess.PIPE, text=True
                )
                if result.returncode != 0:
                    raise RuntimeError(f"MAFFT failed: {result.stderr}")
            return AlignIO.read(aligned_path, "fasta")
        finally:
            Path(tmp_path).unlink(missing_ok=True)
            Path(aligned_path).unlink(missing_ok=True)

    def _column_conservation(self, alignment) -> np.ndarray:
        """Compute per-column conservation (most-common-base count / total sequences)."""
        aln_len = alignment.get_alignment_length()
        num_seqs = len(alignment)
        cons = np.zeros(aln_len)

        for col in range(aln_len):
            bases = [
                alignment[i, col].upper()
                for i in range(num_seqs)
                if alignment[i, col] != "-"
            ]
            if not bases:
                continue
            cons[col] = Counter(bases).most_common(1)[0][1] / num_seqs

        return cons

    def _find_conserved_runs(
        self, conservation: np.ndarray, min_length: int
    ) -> List[tuple]:
        """Find contiguous stretches of columns meeting the conservation threshold."""
        above = conservation >= self.min_conservation
        runs = []
        i = 0
        while i < len(above):
            if above[i]:
                j = i
                while j < len(above) and above[j]:
                    j += 1
                if j - i >= min_length:
                    runs.append((i, j))
                i = j
            else:
                i += 1
        return runs

    def _extract_region(
        self, alignment, seq_ids, total_sequences, col_start, col_end
    ) -> Optional[ConservedRegion]:
        """Extract ungapped sequences and original positions for a column range."""
        region_seqs = []
        positions = {}
        window_len = col_end - col_start

        for i, record in enumerate(alignment):
            subseq = str(record.seq[col_start:col_end])
            ungapped = subseq.replace("-", "").upper()

            if len(ungapped) < window_len * 0.8:
                continue

            region_seqs.append(ungapped)
            gaps_before = str(record.seq[:col_start]).count("-")
            positions[seq_ids[i]] = col_start - gaps_before

        min_required = max(1, math.ceil(self.min_conservation * total_sequences))
        if len(region_seqs) < min_required:
            return None

        consensus = self._create_consensus(region_seqs)
        return ConservedRegion(
            start=0,
            end=len(consensus),
            sequences=region_seqs,
            conservation=len(region_seqs) / total_sequences,
            consensus=consensus,
            positions=positions,
        )

    def _create_consensus(self, sequences: List[str]) -> str:
        """Create consensus with IUPAC ambiguity codes."""
        if not sequences:
            return ""

        min_len = min(len(seq) for seq in sequences)
        iupac_codes = {
            frozenset(["A"]): "A",
            frozenset(["C"]): "C",
            frozenset(["G"]): "G",
            frozenset(["T"]): "T",
            frozenset(["A", "G"]): "R",
            frozenset(["C", "T"]): "Y",
            frozenset(["G", "C"]): "S",
            frozenset(["A", "T"]): "W",
            frozenset(["G", "T"]): "K",
            frozenset(["A", "C"]): "M",
            frozenset(["C", "G", "T"]): "B",
            frozenset(["A", "G", "T"]): "D",
            frozenset(["A", "C", "T"]): "H",
            frozenset(["A", "C", "G"]): "V",
            frozenset(["A", "C", "G", "T"]): "N",
        }

        consensus = []
        for i in range(min_len):
            bases = set(
                seq[i].upper() for seq in sequences if i < len(seq) and seq[i] != "N"
            )
            consensus.append(iupac_codes.get(frozenset(bases), "N"))

        return "".join(consensus)

    def _remove_redundant(
        self, regions: List[ConservedRegion]
    ) -> List[ConservedRegion]:
        """Remove regions with overlapping positions in the reference sequence."""
        if not regions:
            return []

        regions.sort(key=lambda r: r.conservation, reverse=True)
        kept = []

        for region in regions:
            if not region.positions:
                continue
            ref_id = sorted(region.positions.keys())[0]
            ref_pos = region.positions[ref_id]
            region_len = len(region.sequences[0]) if region.sequences else 0

            is_redundant = False
            for k in kept:
                if ref_id in k.positions:
                    if abs(k.positions[ref_id] - ref_pos) < region_len // 2:
                        is_redundant = True
                        break

            if not is_redundant:
                kept.append(region)

        return kept


class Minimap2BasedFinder:
    """Find conserved regions using minimap2 alignment for contigs/genomes."""

    def __init__(self, min_conservation: float = 0.8, identity_threshold: float = 0.9):
        """
        Initialize minimap2-based conserved region finder.

        Args:
            min_conservation: Minimum fraction of sequences that must be conserved
            identity_threshold: Minimum identity threshold for alignments
        """
        self.min_conservation = min_conservation
        self._alignment_cache = {}
        self.identity_threshold = identity_threshold
        self.logger = logging.getLogger(__name__)

    def find_conserved_regions(
        self,
        sequences: Dict[str, str],
        min_length: int = 15,
        max_length: int = 30,
        preset: str = "asm20",
    ) -> List[ConservedRegion]:
        """
        Find conserved regions using minimap2 all-vs-all alignment.

        Args:
            sequences: Dictionary of sequence_id -> sequence
            min_length: Minimum length of conserved region
            max_length: Maximum length of conserved region
            preset: Minimap2 preset (asm20 for <20% divergence, asm10 for <10%, asm5 for <5%)

        Returns:
            List of conserved regions
        """
        if len(sequences) < 2:
            return []

        with tempfile.TemporaryDirectory() as tmpdir:
            # Write sequences to temporary file
            seq_file = Path(tmpdir) / "sequences.fna"
            with open(seq_file, "w") as f:
                for seq_id, seq in sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")

            if shutil.which("minimap2") is None:
                raise RuntimeError('minimap2 executable not found on PATH.')

            cache_key = hashlib.sha1(
                "||".join(f"{sid}:{seq}" for sid, seq in sorted(sequences.items())).encode()
            ).hexdigest() + f"|{preset}|{min_length}|{max_length}"

            if cache_key in self._alignment_cache:
                return self._alignment_cache[cache_key]

            paf_file = Path(tmpdir) / "alignments.paf"
            cmd = [
                "minimap2",
                f"-x{preset}",
                "-c",
                "--cs",
                str(seq_file),
                str(seq_file),
            ]

            with open(paf_file, "w") as f:
                result = subprocess.run(
                    cmd, stdout=f, stderr=subprocess.PIPE, text=True
                )
                if result.returncode != 0:
                    self.logger.error(f"Minimap2 failed: {result.stderr}")
                    return []

            alignment_blocks = self._parse_paf_alignments(paf_file, sequences)
            conserved_regions = self._find_conserved_blocks(
                alignment_blocks, sequences, min_length, max_length
            )

            self._alignment_cache[cache_key] = conserved_regions

            return conserved_regions

    def _parse_paf_alignments(
        self, paf_file: Path, sequences: Dict[str, str]
    ) -> List[Dict]:
        """Parse PAF file to extract high-identity alignment blocks."""
        alignments = []

        with open(paf_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 12:
                    continue

                query_name = parts[0]
                # query_len = int(parts[1])  # Not used
                query_start = int(parts[2])
                query_end = int(parts[3])
                strand = parts[4]
                target_name = parts[5]
                # target_len = int(parts[6])  # Not used
                target_start = int(parts[7])
                target_end = int(parts[8])
                num_matches = int(parts[9])
                alignment_len = int(parts[10])

                # Skip self-alignments
                if query_name == target_name:
                    continue

                # Calculate identity
                identity = num_matches / alignment_len if alignment_len > 0 else 0

                # Only keep high-identity alignments
                if identity >= self.identity_threshold:
                    # Extract CS tag if available
                    cs_tag = None
                    for field in parts[12:]:
                        if field.startswith("cs:Z:"):
                            cs_tag = field[5:]
                            break

                    alignments.append(
                        {
                            "query": query_name,
                            "target": target_name,
                            "query_start": query_start,
                            "query_end": query_end,
                            "target_start": target_start,
                            "target_end": target_end,
                            "strand": strand,
                            "identity": identity,
                            "cs_tag": cs_tag,
                        }
                    )

        return alignments

    def _find_conserved_blocks(
        self,
        alignments: List[Dict],
        sequences: Dict[str, str],
        min_length: int,
        max_length: int,
    ) -> List[ConservedRegion]:
        """Find blocks conserved across multiple sequences from pairwise alignments."""
        conserved_regions = []

        # Group alignments by reference sequence
        seq_list = list(sequences.keys())
        if not seq_list:
            return []

        # Use first sequence as reference
        ref_seq_id = seq_list[0]
        ref_seq = sequences[ref_seq_id]
        ref_len = len(ref_seq)

        # Find regions where multiple sequences align to reference
        coverage = np.zeros(ref_len, dtype=int)
        position_map = defaultdict(
            lambda: defaultdict(list)
        )  # ref_pos -> {seq_id: [target_positions]}

        for aln in alignments:
            if aln["query"] == ref_seq_id:
                # Reference is query
                for i in range(aln["query_start"], aln["query_end"]):
                    coverage[i] += 1
                    target_pos = aln["target_start"] + (i - aln["query_start"])
                    if aln["strand"] == "-":
                        target_pos = aln["target_end"] - (i - aln["query_start"]) - 1
                    position_map[i][aln["target"]].append(target_pos)

            elif aln["target"] == ref_seq_id:
                # Reference is target
                for i in range(aln["target_start"], aln["target_end"]):
                    coverage[i] += 1
                    query_pos = aln["query_start"] + (i - aln["target_start"])
                    position_map[i][aln["query"]].append(query_pos)

        # Find regions with high coverage
        min_sequences = max(1, math.ceil(self.min_conservation * len(sequences)))

        # Scan for conserved regions
        i = 0
        while i < ref_len - min_length:
            if (
                coverage[i] >= min_sequences - 1
            ):  # -1 because reference itself isn't counted
                # Extend region while coverage is high
                j = i
                while (
                    j < min(i + max_length, ref_len)
                    and coverage[j] >= min_sequences - 1
                ):
                    j += 1

                region_length = j - i
                if min_length <= region_length <= max_length:
                    # Extract sequences for this region
                    region_seqs = [ref_seq[i:j]]
                    positions = {ref_seq_id: i}

                    # Get corresponding sequences from other genomes
                    for seq_id in sequences:
                        if seq_id == ref_seq_id:
                            continue

                        # Find most common position mapping for this region
                        pos_counter = Counter()
                        for pos in range(i, j):
                            if seq_id in position_map[pos]:
                                for target_pos in position_map[pos][seq_id]:
                                    pos_counter[target_pos - (pos - i)] += 1

                        if pos_counter:
                            # Use most common starting position
                            start_pos = pos_counter.most_common(1)[0][0]
                            if 0 <= start_pos < len(sequences[seq_id]) - region_length:
                                region_seqs.append(
                                    sequences[seq_id][
                                        start_pos : start_pos + region_length
                                    ]
                                )
                                positions[seq_id] = start_pos

                    # Check if enough sequences have this region
                    if len(region_seqs) >= min_sequences:
                        consensus = self._create_consensus(region_seqs)
                        conservation = len(region_seqs) / len(sequences)

                        conserved_regions.append(
                            ConservedRegion(
                                start=i,
                                end=j,
                                sequences=region_seqs,
                                conservation=conservation,
                                consensus=consensus,
                                positions=positions,
                            )
                        )

                i = j
            else:
                i += 1

        return conserved_regions

    def _create_consensus(self, sequences: List[str]) -> str:
        """Create consensus sequence with IUPAC ambiguity codes."""
        if not sequences:
            return ""

        consensus = []
        seq_length = len(sequences[0])

        iupac_codes = {
            frozenset(["A"]): "A",
            frozenset(["C"]): "C",
            frozenset(["G"]): "G",
            frozenset(["T"]): "T",
            frozenset(["A", "G"]): "R",
            frozenset(["C", "T"]): "Y",
            frozenset(["G", "C"]): "S",
            frozenset(["A", "T"]): "W",
            frozenset(["G", "T"]): "K",
            frozenset(["A", "C"]): "M",
            frozenset(["C", "G", "T"]): "B",
            frozenset(["A", "G", "T"]): "D",
            frozenset(["A", "C", "T"]): "H",
            frozenset(["A", "C", "G"]): "V",
            frozenset(["A", "C", "G", "T"]): "N",
        }

        for i in range(seq_length):
            bases = set()
            for seq in sequences:
                if i < len(seq):
                    bases.add(seq[i].upper())

            consensus.append(iupac_codes.get(frozenset(bases), "N"))

        return "".join(consensus)
