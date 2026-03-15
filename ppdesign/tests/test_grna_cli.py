#!/usr/bin/env python3
"""Tests for CRISPR guide RNA CLI interface."""

import pytest
from pathlib import Path
import tempfile
import shutil
from typer.testing import CliRunner
from ppdesign.probedesign_grna import app
import pandas as pd
from Bio import SeqIO


class TestGRNACLI:
    """Test guide RNA CLI interface."""

    @pytest.fixture
    def runner(self):
        """Create CLI test runner."""
        return CliRunner()

    @pytest.fixture
    def test_sequences(self):
        """Create test sequence directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            seq_dir = Path(tmpdir) / "sequences"
            seq_dir.mkdir()

            # Create test sequences with conserved regions
            test_fasta = seq_dir / "test.fna"
            with open(test_fasta, "w") as f:
                # Common guide region with NGG PAM
                common = "ATCGATCGATCGATCGATCGTGG"
                f.write(f">seq1\nGCGCGC{common}GCGCGC\n")
                f.write(f">seq2\nATATAT{common}ATATAT\n")
                f.write(f">seq3\nTTTTTT{common}TTTTTT\n")

            yield seq_dir

    def test_basic_cli_run(self, runner, test_sequences):
        """Test basic CLI execution."""
        with tempfile.TemporaryDirectory():
            result = runner.invoke(
                app,
                [
                    "--fasta-dir",
                    str(test_sequences),
                    "--output-dir",
                    "test_output",
                    "--conservation",
                    "0.5",
                    "--max-degenerate",
                    "2",
                ],
            )

            # Check command executed successfully
            assert result.exit_code == 0
            assert "Loading sequences" in result.output
            assert "Found" in result.output
            assert "guide rna" in result.output.lower()

    def test_output_files_created(self, runner, test_sequences):
        """Test that all expected output files are created."""
        with tempfile.TemporaryDirectory():
            # Run CLI
            result = runner.invoke(
                app,
                [
                    "--fasta-dir",
                    str(test_sequences),
                    "--output-dir",
                    "test_output",
                    "--conservation",
                    "0.5",
                ],
            )

            assert result.exit_code == 0

            # Check output files exist
            output_dir = Path("results") / "test_output"
            assert output_dir.exists()
            assert (output_dir / "guide_rnas.csv").exists()
            assert (output_dir / "guide_rnas.fasta").exists()
            assert (output_dir / "summary.txt").exists()
            assert (output_dir / "target_mapping.tsv").exists()

            # Clean up
            if output_dir.exists():
                shutil.rmtree(output_dir)

    def test_csv_output_format(self, runner, test_sequences):
        """Test CSV output format and content."""
        with tempfile.TemporaryDirectory():
            result = runner.invoke(
                app,
                [
                    "--fasta-dir",
                    str(test_sequences),
                    "--output-dir",
                    "test_csv",
                    "--conservation",
                    "0.5",
                ],
            )

            assert result.exit_code == 0

            # Read CSV output
            csv_file = Path("results") / "test_csv" / "guide_rnas.csv"
            assert csv_file.exists()

            df = pd.read_csv(csv_file)

            # Check required columns
            required_cols = [
                "Guide_ID",
                "Sequence",
                "PAM",
                "Strand",
                "Conservation",
                "Target_Count",
                "Degenerate_Bases",
                "Specificity_Score",
                "Quality_Score",
            ]
            for col in required_cols:
                assert col in df.columns

            # Check data types and values
            assert len(df) > 0  # Should find at least one guide
            assert df["Guide_ID"].str.startswith("gRNA_").all()
            assert df["PAM"].isin(["NGG", "NAG", "TGG", "AGG", "CGG"]).all()
            assert df["Strand"].isin(["+", "-"]).all()

            # Clean up
            shutil.rmtree(Path("results") / "test_csv")

    def test_fasta_output_format(self, runner, test_sequences):
        """Test FASTA output format."""
        with tempfile.TemporaryDirectory():
            result = runner.invoke(
                app,
                [
                    "--fasta-dir",
                    str(test_sequences),
                    "--output-dir",
                    "test_fasta",
                    "--conservation",
                    "0.5",
                ],
            )

            assert result.exit_code == 0

            # Read FASTA output
            fasta_file = Path("results") / "test_fasta" / "guide_rnas.fasta"
            assert fasta_file.exists()

            records = list(SeqIO.parse(fasta_file, "fasta"))
            assert len(records) > 0

            # Check FASTA format
            for record in records:
                assert record.id.startswith("gRNA_")
                assert len(str(record.seq)) == 20  # Guide should be 20bp
                assert "PAM=" in record.description
                assert "conservation=" in record.description

            # Clean up
            shutil.rmtree(Path("results") / "test_fasta")

    def test_custom_pam_patterns(self, runner, test_sequences):
        """Test custom PAM pattern specification."""
        with tempfile.TemporaryDirectory():
            result = runner.invoke(
                app,
                [
                    "--fasta-dir",
                    str(test_sequences),
                    "--output-dir",
                    "test_custom_pam",
                    "--pam-type",
                    "TGG,AGG",
                    "--conservation",
                    "0.5",
                ],
            )

            # Should complete successfully
            assert result.exit_code == 0
            assert "Using PAM patterns: TGG, AGG" in result.output

            # Clean up
            output_dir = Path("results") / "test_custom_pam"
            if output_dir.exists():
                shutil.rmtree(output_dir)

    def test_conservation_threshold(self, runner, test_sequences):
        """Test conservation threshold filtering."""
        # High conservation threshold
        result_high = runner.invoke(
            app,
            [
                "--fasta-dir",
                str(test_sequences),
                "--output-dir",
                "test_high_cons",
                "--conservation",
                "1.0",  # 100% conservation required
            ],
        )

        # With 100% conservation, might find fewer or no guides
        assert "conserved guide rnas" in result_high.output.lower()

        # Clean up
        output_dir = Path("results") / "test_high_cons"
        if output_dir.exists():
            shutil.rmtree(output_dir)

    def test_top_n_filtering(self, runner, test_sequences):
        """Test limiting output to top N guides."""
        with tempfile.TemporaryDirectory():
            result = runner.invoke(
                app,
                [
                    "--fasta-dir",
                    str(test_sequences),
                    "--output-dir",
                    "test_top_n",
                    "--conservation",
                    "0.3",
                    "--top-n",
                    "5",
                ],
            )

            assert result.exit_code == 0

            # Check CSV has at most 5 guides
            csv_file = Path("results") / "test_top_n" / "guide_rnas.csv"
            if csv_file.exists():
                df = pd.read_csv(csv_file)
                assert len(df) <= 5

            # Clean up
            shutil.rmtree(Path("results") / "test_top_n")

    def test_verbose_mode(self, runner, test_sequences):
        """Test verbose output mode."""
        result = runner.invoke(
            app,
            [
                "--fasta-dir",
                str(test_sequences),
                "--output-dir",
                "test_verbose",
                "--conservation",
                "0.5",
                "--verbose",
            ],
        )

        assert result.exit_code == 0
        # Verbose mode should show more detailed logging
        assert len(result.output) > 0

        # Clean up
        output_dir = Path("results") / "test_verbose"
        if output_dir.exists():
            shutil.rmtree(output_dir)

    def test_error_no_sequences(self, runner):
        """Test error handling when no sequences found."""
        with tempfile.TemporaryDirectory() as tmpdir:
            empty_dir = Path(tmpdir) / "empty"
            empty_dir.mkdir()

            result = runner.invoke(
                app, ["--fasta-dir", str(empty_dir), "--output-dir", "test_error"]
            )

            assert result.exit_code == 1
            assert "No FASTA files found" in result.output

    def test_error_single_sequence(self, runner):
        """Test error when only one sequence provided."""
        with tempfile.TemporaryDirectory() as tmpdir:
            seq_dir = Path(tmpdir) / "single"
            seq_dir.mkdir()

            # Create single sequence
            test_fasta = seq_dir / "single.fna"
            with open(test_fasta, "w") as f:
                f.write(">seq1\nATCGATCGATCGTGGATCG\n")

            result = runner.invoke(
                app, ["--fasta-dir", str(seq_dir), "--output-dir", "test_single"]
            )

            assert result.exit_code == 1
            assert "Need at least 2 sequences" in result.output

    def test_summary_statistics(self, runner, test_sequences):
        """Test summary statistics generation."""
        with tempfile.TemporaryDirectory():
            result = runner.invoke(
                app,
                [
                    "--fasta-dir",
                    str(test_sequences),
                    "--output-dir",
                    "test_summary",
                    "--conservation",
                    "0.5",
                ],
            )

            assert result.exit_code == 0

            # Check summary file
            summary_file = Path("results") / "test_summary" / "summary.txt"
            assert summary_file.exists()

            with open(summary_file, "r") as f:
                content = f.read()
                assert "CRISPR Guide RNA Design Summary" in content
                assert "Input Parameters" in content
                assert "Results" in content
                assert "Conservation Threshold" in content

            # Clean up
            shutil.rmtree(Path("results") / "test_summary")
