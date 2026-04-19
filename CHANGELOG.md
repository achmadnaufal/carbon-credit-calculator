# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Monte Carlo uncertainty propagation module (`src/uncertainty.py`)
  that samples DBH, allometric residual, carbon fraction, and
  root-to-shoot ratio to produce a confidence interval on total tCO2e.
- `UncertaintyConfig` dataclass with configurable iterations, CVs,
  confidence level, and random seed.
- `propagate_tree` and `propagate_population` helpers returning an
  `UncertaintyResult` with mean, median, std, CI, and raw samples.
- Chave et al. (2014) pan-tropical allometric model via
  `estimate_agb_chave2014(dbh, height, wood_density)`.
- `tree_co2e` single-tree convenience function with method selector.
- `CarbonCreditCalculator.compute_tree_co2e` adds per-tree `agb_kg`,
  `carbon_kg`, `co2e_tonnes` columns.
- `CarbonCreditCalculator.compute_credits` returns project totals,
  per-species breakdown, optional Monte Carlo uncertainty block, and
  VCS-style buffer pool deduction (`buffer_pool_pct` config key).
- Wood density (`wd`) and four additional species (Eucalyptus
  urophylla, Swietenia macrophylla, Gliricidia sepium, Paraserianthes
  falcataria) in `SPECIES_PARAMS`.
- Case-insensitive species lookup with generic-tropical fallback.
- Google-style docstrings (Args, Returns, Raises, Example) on all
  public functions.
- 70-case pytest suite covering allometrics, conversions, Monte
  Carlo propagation, validation, and end-to-end pipeline.
- `demo/sample_data.csv` — 20-row realistic tropical/agroforestry
  inventory for quickstart.
- Public package API in `src/__init__.py`.

### Changed
- `CarbonCreditCalculator.analyze` now returns carbon-specific totals
  (delegates to `compute_credits`) instead of a generic numeric summary.
- `load_data` raises `FileNotFoundError` for missing files and
  `ValueError` for unsupported extensions or empty files.
- `validate` now enforces required columns (`species`, `dbh_cm`),
  positive DBH/height, non-NaN species, and raises `ValueError` with
  specific messages.
- Example script (`examples/basic_usage.py`) updated to use the demo
  inventory and show Monte Carlo output.
- README expanded with methodology, data format, scope disclaimer,
  and quickstart loading `demo/sample_data.csv`.

### Fixed
- `preprocess` no longer relies on in-place mutation order; returns a
  defensive copy so callers' DataFrames are never mutated.

## [0.1.0] - 2024

### Added
- Initial implementation: CSV/Excel loader, generic numeric summary,
  basic species power-law allometry, sample data generator.
