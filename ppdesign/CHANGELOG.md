# Changelog

## [0.11.0] - 2026-02-02

### Added
- **Optional MAFFT alignment refinement** for conserved regions via `--align-regions` flag
  - Enables nucleotide-level alignment of regions found by k-mer search
  - Calculates column-wise conservation metrics
  - Extracts ungapped, conserved segments for optimal primer design
  - Gracefully falls back if MAFFT not installed
  - Default: OFF (fast k-mer mode), enable for highest quality
- **Flexible degeneracy control** via `--max-degenerate-positions` parameter
  - Controls how many IUPAC degenerate codes allowed per primer
  - Range: 0 (strict ACGT only) to 5 (flexible)
  - Default: 2 (allows minimal variation)
  - Validates all IUPAC codes (RYMKSWBDHVN)
- **Individual sequence extraction** for primer candidates
  - Primers now extracted from each sequence in region, not just consensus
  - Handles regions with IUPAC codes gracefully
  - Generates more primer candidates per region
  - Maintains "near-perfect match" requirement

### Changed
- **BREAKING**: Primer extraction strategy changed from consensus-only to individual sequences
  - Before: Extracted primers from consensus sequence (failed if IUPAC codes present)
  - After: Extracts from all individual sequences in region
  - Impact: More primers generated, especially for variable regions
  - Migration: No code changes needed, output may have more primers
- Primer sequence validation now accepts IUPAC codes up to `max_degenerate_positions` threshold
- Conserved region finding now accepts `align_regions` parameter (default: False)

### Documentation
- Updated `docs/primer_pipeline.md` with alignment mode usage guide
- Added performance comparison (k-mer vs MAFFT)
- Added "when to use alignment" decision guide
- Added degeneracy control recommendations
- Added MAFFT installation instructions
- Added example commands for fast vs high-quality modes

### Testing
- Added `tests/test_conserved_finder.py` with 7 new tests for alignment functionality
- Added `tests/test_primer_candidates.py` with 11 new tests for IUPAC handling
- Added 3 new integration tests in `test_primer_cli.py` for new CLI flags
- All tests pass with new features

## [0.10.1] - 2026-02-02

### Fixed
- **CRITICAL**: Primer position calculation now uses actual genomic coordinates instead of relative positions
  - All forward primers were incorrectly reporting position=0
  - All reverse primers were reporting position=region_length
  - Fix: Use actual genomic positions from `region.positions` dict
  - Impact: Amplicon sizes and position-based pairing now work correctly
- **CRITICAL**: Primer pairing now validates that forward and reverse primers share at least one target sequence
  - Before: Invalid pairs could be created (forward from seq1 + reverse from seq2)
  - Fix: Added target sequence overlap check in `_is_valid_pair()`
  - Impact: All output primer pairs are now guaranteed to amplify valid targets
- **CRITICAL**: Conserved region length requirements now match amplicon size constraints
  - Before: KmerBasedFinder used max_length=30bp but amplicon_min=100bp (impossible constraint)
  - Fix: Automatically calculate min/max region lengths from amplicon parameters
  - Formula: min_region_length = amplicon_min + primer_max_length
  - Impact: Default parameters now work without manual tuning
- **HIGH**: Nucleotide pipeline dependency check now correctly looks for `mafft` instead of `muscle`
  - Consistent with actual MSA implementation
- **MEDIUM**: Nucleotide pipeline now creates output directory if it doesn't exist
  - Ensures consistent behavior with other pipelines

### Changed
- Conserved region length parameters in primer pipeline now automatically calculated from amplicon size requirements
- Added parameter validation display showing min conserved region needed based on amplicon/primer constraints

### Testing
- Added 3 new unit tests for critical bug fixes (position diversity, target overlap, parameter validation)
- Added integration test for default parameter usage
- Updated E2E verification script to test default parameters
- All 31 unit tests and 10 integration tests now pass

## [0.10.0] - 2026-02-02

### Added
- **PCR primer pair design module** (`ppdesign primer`) for forward/reverse primer pairing
  - Position-based pairing optimization (O(N) vs O(N²) performance)
  - Amplicon size control (100bp-2kb, supports long-read applications)
  - Primer Tm matching validation (within 5°C by default)
  - Cross-dimer checking between primer pairs
  - Quality scoring system (0-100 scale) based on conservation, Tm, GC, and dimer metrics
- New primer design modules (all <300 lines, KISS-compliant):
  - `primer_types.py`: Data structures with 0-based coordinate system
  - `primer_candidates.py`: Forward/reverse primer generation from conserved regions
  - `primer_pairing.py`: Efficient position-indexed pairing algorithm
  - `primer_validation.py`: PCR best practices validation (poly-X, GC clamp, 3' stability)
  - `primer_scoring.py`: Composite quality scoring and ranking
  - `primer_output.py`: CSV/FASTA/TSV/summary output writers
- Shared `output_utils.py` helper for consistent output directory handling
- New documentation: primer pipeline usage and parameters
- Test dataset for primer design validation (`resources/test_data/primer/`)

### Fixed
- **CLI documentation**: All commands now use correct unified syntax (`ppdesign <subcommand>`)
  - Updated README.md and all docs files
  - Fixed examples in unified_pipeline.md, nucleotide_pipeline.md, post_processing.md, results.md, overview.md
- **Database documentation**: Clarified NeLLi database requires JGI infrastructure (set PPDESIGN_DUCKDB_PATH)
- **Output directory handling**: Consistent behavior across all pipelines
  - Absolute paths: preserved as-is
  - Relative paths: automatically prefixed with `results/`
  - All pipelines now use shared `resolve_output_dir()` function
- **Dependency checks**: Added pre-flight validation to nucleotide pipeline
  - Checks for minimap2 availability (minimap2 method)
  - Checks for muscle availability (MSA method)
- Removed non-existent `pixi run jupyter-tailscale` reference from README

### Changed
- Unified CLI now includes `ppdesign primer` subcommand
- All pipelines (unified, nucleotide, grna) now use shared output directory helper
- Updated architecture to show 3 distinct modes:
  1. Gene-based probe design (unified)
  2. Direct nucleotide probe design (nucleotide)
  3. PCR primer pair design (primer) ← NEW
  4. CRISPR guide RNA design (grna)

## [0.9.0] - 2025-10-05

### Added
- Published consolidated release notes highlighting the PPDesign probe, primer,
  and gRNA design workflows.
- Introduced a release badge in the project README so GitHub always displays the
  latest tagged version at a glance.

### Changed
- Bumped the project version metadata (`pixi.toml`, `ppdesign.__init__`) to 0.9.0
  in preparation for tagging the release.
- Updated developer instructions to use the Pixi `dev` feature, ensuring test and
  lint commands run with a fully provisioned toolchain.
- Archived superseded codon-alignment and CRISPR helper scripts under
  `legacy/modules/` to keep the active pipeline surface lean.

### Fixed
- Clarified how to invoke the automated regression suite so new releases are
  validated consistently before publishing.

---

See `ppdesign/README.md` for detailed usage guides and
`docs/` for pipeline-specific documentation.
