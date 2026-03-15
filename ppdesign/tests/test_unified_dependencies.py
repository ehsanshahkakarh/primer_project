import pytest

from ppdesign.probedesign_unified import UnifiedProbeDesign


def test_unified_pipeline_dependency_check(monkeypatch):
    def fake_which(tool):
        return None

    monkeypatch.setattr("ppdesign.probedesign_unified.shutil.which", fake_which)

    with pytest.raises(RuntimeError) as excinfo:
        UnifiedProbeDesign(output_dir="tmp_unified", threads=1)

    assert "Missing required tools" in str(excinfo.value)
