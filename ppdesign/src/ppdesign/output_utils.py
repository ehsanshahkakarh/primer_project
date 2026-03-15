"""Output directory utilities for consistent path handling across pipelines."""

from pathlib import Path


def resolve_output_dir(output_dir: str | Path) -> Path:
    """
    Resolve output directory consistently across all pipelines.

    Rules:
    - Absolute paths: use as-is (preserves user intent)
    - Relative paths: prepend 'results/' (standard location)

    Args:
        output_dir: Output directory path (string or Path object)

    Returns:
        Resolved Path object

    Examples:
        >>> resolve_output_dir("/absolute/path")
        Path("/absolute/path")

        >>> resolve_output_dir("relative/path")
        Path("results/relative/path")

        >>> resolve_output_dir(Path("/absolute/path"))
        Path("/absolute/path")
    """
    output_path = Path(output_dir)
    return output_path if output_path.is_absolute() else Path("results") / output_path
