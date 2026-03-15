"""Integration tests for MSA-first primer pipeline on real Burkholderia 16S data."""

import csv
import os
import shutil
import subprocess
import sys

import pytest
from pathlib import Path


BURKHOLDERIA_FASTA = Path(__file__).parent / "data" / "burkholderia_16s.fna"
MAFFT_AVAILABLE = shutil.which("mafft") is not None
def _qinsi_functional():
    """Check if mafft-qinsi actually works (needs MXSCARNA extension)."""
    if shutil.which("mafft-qinsi") is None:
        return False
    import tempfile
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(">s1\nACGT\n>s2\nACGT\n")
            tmp = f.name
        result = subprocess.run(
            ["mafft-qinsi", "--quiet", tmp],
            capture_output=True, text=True, timeout=30,
        )
        Path(tmp).unlink(missing_ok=True)
        return result.returncode == 0
    except Exception:
        return False

QINSI_AVAILABLE = _qinsi_functional()


@pytest.mark.skipif(not MAFFT_AVAILABLE, reason="MAFFT not installed")
@pytest.mark.skipif(
    not BURKHOLDERIA_FASTA.exists(), reason="Burkholderia 16S test data not found"
)
class TestPrimerMSAPipeline:
    """Test MSA-based primer design on real Burkholderia 16S sequences."""

    def _run_primer_cli(self, tmp_path, extra_args=None):
        """Run the primer CLI and return (result, output_dir)."""
        output_dir = tmp_path / "output"
        cmd = [
            sys.executable,
            "-m",
            "ppdesign.cli",
            "primer",
            "main",
            "--fasta-input",
            str(BURKHOLDERIA_FASTA),
            "--output-dir",
            str(output_dir),
        ]
        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(
            cmd,
            env={"PYTHONPATH": "src", "PATH": os.environ.get("PATH", "")},
            capture_output=True,
            text=True,
        )
        return result, output_dir

    def _read_primer_pairs(self, output_dir):
        """Read primer_pairs.csv and return list of dicts."""
        csv_path = output_dir / "primer_pairs.csv"
        if not csv_path.exists():
            return []
        with open(csv_path) as f:
            return list(csv.DictReader(f))

    def test_short_amplicon_300bp(self, tmp_path):
        """300bp amplicon from 16S -- should find primers in V3-V4 region."""
        result, output_dir = self._run_primer_cli(
            tmp_path,
            [
                "--amplicon-min", "200",
                "--amplicon-max", "400",
                "--conservation", "0.8",
            ],
        )

        assert result.returncode == 0, (
            f"Pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )
        assert "Pipeline Complete" in result.stdout

        pairs = self._read_primer_pairs(output_dir)
        assert len(pairs) >= 1, "Expected at least 1 primer pair for 200-400bp amplicon"

        for pair in pairs:
            amplicon = int(pair["amplicon_size"])
            assert 200 <= amplicon <= 400, f"Amplicon {amplicon}bp outside 200-400 range"

    def test_long_amplicon_1500bp(self, tmp_path):
        """~1.5kb amplicon -- near-full-length 16S primers."""
        result, output_dir = self._run_primer_cli(
            tmp_path,
            [
                "--amplicon-min", "1300",
                "--amplicon-max", "1600",
                "--conservation", "0.8",
            ],
        )

        assert result.returncode == 0, (
            f"Pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )

        pairs = self._read_primer_pairs(output_dir)
        assert len(pairs) >= 1, "Expected at least 1 primer pair for 1300-1600bp amplicon"

        for pair in pairs:
            amplicon = int(pair["amplicon_size"])
            assert 1300 <= amplicon <= 1600, (
                f"Amplicon {amplicon}bp outside 1300-1600 range"
            )

    @pytest.mark.skipif(not QINSI_AVAILABLE, reason="mafft-qinsi not installed")
    def test_structure_aware(self, tmp_path):
        """Q-INS-i mode for rRNA."""
        result, output_dir = self._run_primer_cli(
            tmp_path,
            [
                "--structure-aware",
                "--amplicon-min", "200",
                "--amplicon-max", "400",
                "--conservation", "0.8",
            ],
        )

        assert result.returncode == 0, (
            f"Q-INS-i pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )
        assert "structure-aware" in result.stdout

    def test_kmer_fallback(self, tmp_path):
        """Legacy k-mer mode runs without errors.

        K-mer mode may not find amplicon-length conserved regions in diverse
        16S data (that's exactly why MSA mode exists), so we just verify the
        mode activates and doesn't crash.
        """
        result, output_dir = self._run_primer_cli(
            tmp_path,
            [
                "--method", "kmer",
                "--amplicon-min", "200",
                "--amplicon-max", "400",
                "--conservation", "0.6",
            ],
        )

        assert "k-mer" in result.stdout
        # K-mer mode may legitimately find 0 regions for diverse 16S data;
        # exit 1 with "No conserved regions found" is acceptable
        assert result.returncode == 0 or "No conserved regions found" in result.stdout
