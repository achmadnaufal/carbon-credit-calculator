# Carbon Credit Calculator

A Python tool for estimating carbon sequestration and tradable carbon
credits from tree-planting projects. It converts field-measured tree
dimensions into above-ground biomass (AGB), root biomass, carbon
stock, and CO2-equivalent (tCO2e) using peer-reviewed allometric
equations, with optional Monte Carlo uncertainty propagation.

## What it does

1. **Loads** a tree inventory from CSV or Excel (`species`, `dbh_cm`,
   optional `height_m`, `age_years`, `plot_id`, `date_measured`).
2. **Validates** the data — rejects empty frames, missing columns,
   non-positive DBH/height, NaNs, and missing species.
3. **Computes per-tree** AGB, total biomass (shoots + roots via
   species-specific root-to-shoot ratios), carbon mass, and tCO2e.
4. **Aggregates** to project-level totals with a per-species
   breakdown.
5. **Propagates uncertainty** through all coefficients and returns a
   confidence interval on total tCO2e (Monte Carlo).
6. **Deducts** an optional VCS-style buffer pool to report net
   tradable credits.

## Features

- Species-specific power-law allometry (`AGB = a * DBH^b`) for common
  agroforestry and tropical plantation species.
- Chave et al. (2014) pan-tropical allometric model using wood density
  and tree height (optional, set `allometry_method="chave2014"`).
- Root-to-shoot ratios per species (Mokany et al. 2006).
- IPCC (2006) default carbon fraction of 0.47 (configurable).
- Stoichiometric 44/12 C to CO2e conversion.
- Monte Carlo uncertainty propagation on DBH measurement, allometric
  residual, carbon fraction, and root-to-shoot ratio. Configurable
  sample size, confidence level, and random seed.
- VCS-style buffer pool deduction on gross credits.
- Google-style docstrings on all public functions.

## Install

```bash
pip install -r requirements.txt
```

Runs on Python 3.9+.

## Quickstart

```python
from src import CarbonCreditCalculator, UncertaintyConfig

calc = CarbonCreditCalculator({"buffer_pool_pct": 0.20})
cfg = UncertaintyConfig(n_iterations=5000, random_seed=42)

result = calc.run(
    "demo/sample_data.csv",
    uncertainty=True,
    uncertainty_config=cfg,
)

print(f"Gross tCO2e:          {result['total_tco2e']:.4f}")
print(f"Buffer pool (20%):    {result['buffer_tco2e']:.4f}")
print(f"Net tradable credits: {result['net_credits_tco2e']:.4f}")

u = result["uncertainty"]
print(f"90% CI on tCO2e: [{u['ci_lower_tco2e']:.4f}, {u['ci_upper_tco2e']:.4f}]")
```

Or run the bundled example:

```bash
python examples/basic_usage.py
```

## Data format

Required columns: `species`, `dbh_cm`. Recommended additional
columns: `tree_id`, `height_m`, `age_years`, `plot_id`,
`date_measured`. See [`demo/sample_data.csv`](demo/sample_data.csv)
for a 20-row example.

## Methodology

- **AGB (default):** species-specific power law
  `AGB_kg = a * DBH_cm^b`, coefficients compiled from published
  equations for each species (see `src/allometrics.py`).
- **AGB (Chave 2014):** `AGB = 0.0673 * (rho * D^2 * H)^0.976` where
  `rho` is wood density (g/cm^3), `D` is DBH (cm), `H` is height (m).
- **Below-ground biomass:** AGB x species-specific root-to-shoot
  ratio (default 0.26, IPCC 2006 / Mokany et al. 2006).
- **Carbon:** total biomass x 0.47 (IPCC 2006 Table 4.3).
- **CO2e:** carbon x 44/12 (stoichiometric ratio).

## Monte Carlo uncertainty

Default coefficients of variation (all configurable via
`UncertaintyConfig`):

| Source                  | Default | Reference                                 |
|-------------------------|---------|-------------------------------------------|
| DBH measurement         | CV 5%   | Field protocols                           |
| Allometric residual     | CV 20%  | Chave et al. 2014 (pan-tropical ~36%)     |
| Carbon fraction         | SD 0.02 | IPCC 2006 (0.47 +/- 0.02)                 |
| Root-to-shoot ratio     | CV 15%  | Mokany et al. 2006                        |

Per-iteration biomass is drawn from a lognormal error distribution to
keep samples strictly positive. Per-tree errors are treated as
independent when aggregating across a plot.

## Scope and honesty

This tool is a **first-pass estimation** suitable for feasibility
studies, educational use, and plot-level scoping. It does **not**:

- Account for baseline (without-project) carbon stocks.
- Model tree mortality, harvest leakage, or soil carbon dynamics.
- Implement a specific carbon registry's methodology end-to-end
  (VCS VM0047, Gold Standard AR, etc.). The buffer pool deduction is
  a convenience placeholder, not a registry-compliant calculation.
- Replace field verification or a certified validator/verifier body.

Always verify outputs against your registry's required methodology
before claiming or trading credits.

## Tests

```bash
pytest
```

70+ unit and integration tests cover allometric calculations,
conversions, Monte Carlo behaviour, DataFrame validation, and the
end-to-end pipeline.

## Project structure

```
carbon-credit-calculator/
  src/
    allometrics.py      # AGB equations, carbon/CO2e conversions
    uncertainty.py      # Monte Carlo propagation
    main.py             # CarbonCreditCalculator pipeline
    data_generator.py   # Synthetic data for stress tests
  demo/
    sample_data.csv     # 20-tree example inventory
  examples/
    basic_usage.py      # End-to-end quickstart
  tests/                # pytest suite
  requirements.txt
  CHANGELOG.md
  README.md
```

## License

MIT.
