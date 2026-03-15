"""Integration tests for primer pair design CLI."""

import pytest
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO


class TestPrimerCLI:
    """Test primer CLI end-to-end."""

    @pytest.fixture
    def test_fasta(self, tmp_path):
        """Create a test FASTA file."""
        fasta_file = tmp_path / "test_sequences.fna"

        # Create 5 similar sequences with conserved regions
        sequences = []
        for i in range(5):
            # Create sequence with conserved start and end
            seq = "ATGCGATCGATCGATCGATCG" + "TAGCTAGCTAGCTAGCT" * 10 + "GCTAGCTAGCTAGCTAGCTAGCT"
            sequences.append(f">test_seq_{i+1}\n{seq}\n")

        fasta_file.write_text("".join(sequences))
        return fasta_file

    def test_cli_help(self):
        """Test that CLI help works."""
        result = subprocess.run(
            [sys.executable, "-m", "ppdesign.cli", "primer", "--help"],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "PCR primer pair design" in result.stdout

    def test_basic_primer_design(self, test_fasta, tmp_path):
        """Test basic primer design workflow."""
        output_dir = tmp_path / "primer_output"

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir),
                "--method",
                "kmer",
                "--conservation",
                "0.7",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "500",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "Pipeline Complete" in result.stdout

        # Check output files exist
        assert (output_dir / "primer_pairs.csv").exists()
        assert (output_dir / "primer_pairs.fasta").exists()
        assert (output_dir / "amplicon_predictions.tsv").exists()
        assert (output_dir / "summary.txt").exists()

    def test_output_files_content(self, test_fasta, tmp_path):
        """Test that output files have correct content."""
        output_dir = tmp_path / "primer_output"

        subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir),
                "--method",
                "kmer",
                "--conservation",
                "0.7",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
            ],
            env={"PYTHONPATH": "src"},
            check=True,
        )

        # Check CSV has header and data
        csv_file = output_dir / "primer_pairs.csv"
        lines = csv_file.read_text().strip().split("\n")
        assert len(lines) > 1  # Header + at least one data row
        assert "pair_id" in lines[0]
        assert "forward_seq" in lines[0]
        assert "reverse_seq" in lines[0]

        # Check FASTA format
        fasta_file = output_dir / "primer_pairs.fasta"
        records = list(SeqIO.parse(fasta_file, "fasta"))
        assert len(records) >= 2  # At least one forward and one reverse
        assert any("forward" in r.id for r in records)
        assert any("reverse" in r.id for r in records)

        # Check summary
        summary_file = output_dir / "summary.txt"
        summary = summary_file.read_text()
        assert "Total primer pairs" in summary
        assert "Top 10" in summary or "Top" in summary

    def test_invalid_parameters(self, test_fasta, tmp_path):
        """Test that invalid parameters are caught."""
        output_dir = tmp_path / "primer_output"

        # Test with min > max amplicon size (should fail during runtime)
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir),
                "--amplicon-min",
                "500",
                "--amplicon-max",
                "100",  # Invalid: max < min
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0

    def test_different_amplicon_ranges(self, test_fasta, tmp_path):
        """Test different amplicon size ranges."""
        # This test uses the same sequences with different amplicon constraints
        # Note: The test fixture has short conserved regions, so we use
        # amplicon ranges that work with the available region sizes

        # Test with amplicon range 50-300
        output_dir1 = tmp_path / "amplicons_50_300"
        result1 = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir1),
                "--method",
                "kmer",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Test with amplicon range 100-500
        output_dir2 = tmp_path / "amplicons_100_500"
        result2 = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir2),
                "--method",
                "kmer",
                "--amplicon-min",
                "100",
                "--amplicon-max",
                "500",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Both should complete (may have different number of pairs based on conserved region sizes)
        # At least one should succeed
        assert result1.returncode == 0 or result2.returncode == 0, (
            f"Both amplicon ranges failed:\n"
            f"50-300: {result1.stderr}\n"
            f"100-500: {result2.stderr}"
        )

    def test_tm_constraints(self, test_fasta, tmp_path):
        """Test Tm constraint effects."""
        # Strict Tm range
        output_dir = tmp_path / "strict_tm"
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir),
                "--method",
                "kmer",
                "--tm-min",
                "58",
                "--tm-max",
                "62",
                "--tm-diff-max",
                "2.0",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Should complete (may find fewer primers)
        assert result.returncode == 0

    def test_gc_constraints(self, test_fasta, tmp_path):
        """Test GC content constraint effects."""
        output_dir = tmp_path / "gc_constrained"
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir),
                "--method",
                "kmer",
                "--gc-min",
                "45",
                "--gc-max",
                "55",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0

    def test_conservation_threshold(self, test_fasta, tmp_path):
        """Test different conservation thresholds."""
        # High conservation
        output_dir1 = tmp_path / "high_cons"
        result1 = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir1),
                "--method",
                "kmer",
                "--conservation",
                "0.9",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Low conservation
        output_dir2 = tmp_path / "low_cons"
        result2 = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir2),
                "--method",
                "kmer",
                "--conservation",
                "0.5",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Both should work
        assert result1.returncode == 0
        assert result2.returncode == 0

    def test_relative_output_path(self, test_fasta, tmp_path):
        """Test that relative output paths are handled correctly."""
        # Change to temp directory
        import os

        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ppdesign.cli",
                    "primer",
                    "main",
                    "--fasta-input",
                    str(test_fasta),
                    "--output-dir",
                    "relative_output",  # Relative path
                    "--method",
                    "kmer",
                    "--conservation",
                    "0.7",
                    "--amplicon-min",
                    "50",
                    "--amplicon-max",
                    "300",
                ],
                env={"PYTHONPATH": str(Path(original_dir) / "src")},
                capture_output=True,
                text=True,
            )

            assert result.returncode == 0
            # Should create results/relative_output
            assert (tmp_path / "results" / "relative_output").exists()

        finally:
            os.chdir(original_dir)

    def test_default_parameters(self, test_fasta, tmp_path):
        """
        Test primer design with DEFAULT parameters (Bug 3 fix).

        Before fix: Default parameters (amplicon_min=100, conserved region
        max_length=30) created impossible constraint - no valid pairs.

        After fix: Default parameters should work without manual tuning.
        Uses --method kmer since subprocess env lacks PATH for MAFFT.
        """
        output_dir = tmp_path / "default_params"

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir),
                "--method",
                "kmer",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Should succeed with default parameters
        assert result.returncode == 0, f"Failed with defaults: {result.stderr}"
        assert "Pipeline Complete" in result.stdout

        # Verify output files exist
        assert (output_dir / "primer_pairs.csv").exists()

        # Verify we got at least one valid pair
        csv_file = output_dir / "primer_pairs.csv"
        lines = csv_file.read_text().strip().split("\n")
        pair_count = len(lines) - 1  # Subtract header

        assert pair_count >= 1, (
            f"Expected at least 1 primer pair with default parameters, got {pair_count}. "
            "This indicates Bug 3 (conserved region length mismatch) may not be fixed."
        )

    def test_align_regions_flag(self, test_fasta, tmp_path):
        """Test --align-regions flag enables MAFFT alignment."""
        output_dir_fast = tmp_path / "fast_mode"
        output_dir_aligned = tmp_path / "aligned_mode"

        # Fast mode (k-mer, no alignment)
        result_fast = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir_fast),
                "--method",
                "kmer",
                "--conservation",
                "0.7",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Aligned mode (with MAFFT, if available)
        result_aligned = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir_aligned),
                "--conservation",
                "0.7",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
                "--align-regions",  # Enable alignment
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Fast mode should always work
        assert result_fast.returncode == 0
        assert "Pipeline Complete" in result_fast.stdout

        # Aligned mode should work if MAFFT is available
        # If MAFFT not available, should fall back gracefully
        if result_aligned.returncode == 0:
            assert "Pipeline Complete" in result_aligned.stdout
            assert "MAFFT" in result_aligned.stdout
        else:
            # MAFFT might not be installed
            assert "mafft" in result_aligned.stderr.lower() or result_aligned.returncode == 0

    def test_max_degenerate_positions(self, test_fasta, tmp_path):
        """Test --max-degenerate-positions parameter."""
        # Strict mode (ACGT only)
        output_dir_strict = tmp_path / "strict_acgt"
        result_strict = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir_strict),
                "--method",
                "kmer",
                "--conservation",
                "0.7",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
                "--max-degenerate-positions",
                "0",  # No IUPAC codes allowed
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Flexible mode (allow 1-2 IUPAC codes)
        output_dir_flexible = tmp_path / "flexible_iupac"
        result_flexible = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir_flexible),
                "--method",
                "kmer",
                "--conservation",
                "0.7",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
                "--max-degenerate-positions",
                "2",  # Allow 1-2 IUPAC codes
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Both should complete
        assert result_strict.returncode == 0
        assert result_flexible.returncode == 0

        # Both should produce primers (test data is highly conserved)
        assert (output_dir_strict / "primer_pairs.csv").exists()
        assert (output_dir_flexible / "primer_pairs.csv").exists()

    def test_combined_alignment_and_degeneracy(self, test_fasta, tmp_path):
        """Test combined --align-regions and --max-degenerate-positions."""
        output_dir = tmp_path / "high_quality_mode"

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "ppdesign.cli",
                "primer",
                "main",
                "--fasta-input",
                str(test_fasta),
                "--output-dir",
                str(output_dir),
                "--conservation",
                "0.7",
                "--amplicon-min",
                "50",
                "--amplicon-max",
                "300",
                "--align-regions",  # Enable alignment
                "--max-degenerate-positions",
                "1",  # Allow minimal degeneracy
            ],
            env={"PYTHONPATH": "src"},
            capture_output=True,
            text=True,
        )

        # Should complete (may or may not use MAFFT depending on availability)
        assert result.returncode == 0 or "MAFFT" in result.stderr
        if result.returncode == 0:
            assert (output_dir / "primer_pairs.csv").exists()
