"""Monte Carlo uncertainty propagation for carbon stock estimates.

Allometric equations, carbon fractions, and root-to-shoot ratios all
carry meaningful uncertainty. Reporting a single point estimate of
tCO2e without a confidence interval can lead to over-crediting. This
module propagates these uncertainties through the biomass -> carbon ->
CO2e pipeline using Monte Carlo sampling.

Default coefficient-of-variation values are conservative mid-range
estimates drawn from the IPCC 2019 Refinement and Chave et al. (2014).
They can be overridden per call via :class:`UncertaintyConfig`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Tuple

import numpy as np

from .allometrics import (
    CARBON_FRACTION,
    CO2_TO_C,
    get_species_params,
)


@dataclass(frozen=True)
class UncertaintyConfig:
    """Configuration for Monte Carlo uncertainty propagation.

    Attributes:
        n_iterations: Number of Monte Carlo samples. Higher = narrower
            confidence intervals on the CI itself. 10,000 is typical.
        dbh_cv: Coefficient of variation for DBH measurement error.
        allometric_cv: CV for allometric model residual error. Chave
            2014 reports ~0.357 (36%) for pan-tropical moist forests.
        carbon_fraction_sd: Standard deviation (absolute) of the carbon
            fraction. IPCC reports 0.47 +/- 0.02.
        rs_ratio_cv: CV for root-to-shoot ratio (Mokany et al. 2006).
        confidence_level: Two-sided confidence level (e.g. 0.90 -> 5th
            and 95th percentiles).
        random_seed: Seed for reproducibility. None = non-deterministic.
    """

    n_iterations: int = 10_000
    dbh_cv: float = 0.05
    allometric_cv: float = 0.20
    carbon_fraction_sd: float = 0.02
    rs_ratio_cv: float = 0.15
    confidence_level: float = 0.90
    random_seed: int | None = 42

    def __post_init__(self) -> None:
        if self.n_iterations < 100:
            raise ValueError(
                f"n_iterations must be >= 100, got {self.n_iterations}"
            )
        if not (0 < self.confidence_level < 1):
            raise ValueError(
                f"confidence_level must be in (0, 1), got {self.confidence_level}"
            )
        for name, val in [
            ("dbh_cv", self.dbh_cv),
            ("allometric_cv", self.allometric_cv),
            ("carbon_fraction_sd", self.carbon_fraction_sd),
            ("rs_ratio_cv", self.rs_ratio_cv),
        ]:
            if val < 0 or not math.isfinite(val):
                raise ValueError(f"{name} must be >= 0 and finite, got {val}")


@dataclass(frozen=True)
class UncertaintyResult:
    """Result of Monte Carlo propagation for a single tree or population.

    Attributes:
        mean_tco2e: Mean of the CO2e distribution (tonnes).
        median_tco2e: Median of the CO2e distribution (tonnes).
        std_tco2e: Standard deviation of the CO2e distribution.
        ci_lower: Lower bound of the confidence interval.
        ci_upper: Upper bound of the confidence interval.
        confidence_level: Confidence level used (e.g. 0.90).
        n_iterations: Number of Monte Carlo samples drawn.
        samples: Raw sample array (may be downsampled if huge).
    """

    mean_tco2e: float
    median_tco2e: float
    std_tco2e: float
    ci_lower: float
    ci_upper: float
    confidence_level: float
    n_iterations: int
    samples: np.ndarray = field(repr=False)

    def as_dict(self) -> dict:
        """Return a plain-dict view suitable for JSON or DataFrame rows."""
        return {
            "mean_tco2e": round(float(self.mean_tco2e), 6),
            "median_tco2e": round(float(self.median_tco2e), 6),
            "std_tco2e": round(float(self.std_tco2e), 6),
            "ci_lower_tco2e": round(float(self.ci_lower), 6),
            "ci_upper_tco2e": round(float(self.ci_upper), 6),
            "confidence_level": self.confidence_level,
            "n_iterations": self.n_iterations,
        }


def _lognormal_params(mean: float, cv: float) -> Tuple[float, float]:
    """Compute lognormal (mu, sigma) that yield the given mean and CV.

    The lognormal family is preferred for strictly-positive quantities
    like biomass because it cannot produce negative samples.
    """
    variance = (cv * mean) ** 2
    sigma2 = math.log(1.0 + variance / (mean ** 2)) if mean > 0 else 0.0
    sigma = math.sqrt(sigma2)
    mu = math.log(mean) - 0.5 * sigma2 if mean > 0 else 0.0
    return mu, sigma


def propagate_tree(
    species: str,
    dbh_cm: float,
    config: UncertaintyConfig | None = None,
) -> UncertaintyResult:
    """Run Monte Carlo uncertainty propagation for a single tree.

    Samples measurement error on DBH, allometric model residual error,
    carbon fraction, and root-to-shoot ratio, then computes the
    resulting tCO2e distribution.

    Args:
        species: Scientific name (unknown -> generic tropical).
        dbh_cm: Measured DBH in centimetres.
        config: Uncertainty configuration. Default values are used if
            ``None``.

    Returns:
        :class:`UncertaintyResult` with mean, median, std, and CI.

    Raises:
        ValueError: If ``dbh_cm`` is not a positive finite number.

    Example:
        >>> r = propagate_tree("Acacia mangium", 20.0)
        >>> r.ci_lower < r.mean_tco2e < r.ci_upper
        True
    """
    if dbh_cm is None or not math.isfinite(float(dbh_cm)) or float(dbh_cm) <= 0:
        raise ValueError(f"dbh_cm must be > 0 and finite, got {dbh_cm!r}")
    cfg = config or UncertaintyConfig()
    rng = np.random.default_rng(cfg.random_seed)
    params = get_species_params(
        species if isinstance(species, str) else "generic_tropical"
    )

    # DBH measurement error — normal, clipped at a small positive value.
    dbh_samples = rng.normal(dbh_cm, max(dbh_cm * cfg.dbh_cv, 1e-9), cfg.n_iterations)
    dbh_samples = np.clip(dbh_samples, 1e-3, None)

    # Point AGB, then multiplicative lognormal model error.
    agb_point = params["a"] * (dbh_samples ** params["b"])
    if cfg.allometric_cv > 0:
        mu, sigma = _lognormal_params(1.0, cfg.allometric_cv)
        model_err = rng.lognormal(mu, sigma, cfg.n_iterations)
    else:
        model_err = np.ones(cfg.n_iterations)
    agb_samples = agb_point * model_err

    # Root-to-shoot ratio — truncated normal at 0.
    if cfg.rs_ratio_cv > 0:
        rs_samples = rng.normal(
            params["rs"], params["rs"] * cfg.rs_ratio_cv, cfg.n_iterations
        )
        rs_samples = np.clip(rs_samples, 0.0, None)
    else:
        rs_samples = np.full(cfg.n_iterations, params["rs"])

    # Carbon fraction — normal, clipped to (0, 1].
    if cfg.carbon_fraction_sd > 0:
        cf_samples = rng.normal(
            CARBON_FRACTION, cfg.carbon_fraction_sd, cfg.n_iterations
        )
        cf_samples = np.clip(cf_samples, 1e-6, 1.0)
    else:
        cf_samples = np.full(cfg.n_iterations, CARBON_FRACTION)

    total_biomass_kg = agb_samples * (1.0 + rs_samples)
    carbon_kg = total_biomass_kg * cf_samples
    co2e_t = carbon_kg * CO2_TO_C / 1000.0

    return _summarize(co2e_t, cfg)


def propagate_population(
    trees: Iterable[Tuple[str, float]],
    config: UncertaintyConfig | None = None,
) -> UncertaintyResult:
    """Run Monte Carlo propagation across many trees and sum per-iteration.

    Per-tree sources of error are treated as independent, which is the
    standard assumption when trees are sampled from a population. For
    systematic errors (e.g. one biased DBH tape), run each tree
    separately and add the confidence intervals conservatively.

    Args:
        trees: Iterable of ``(species, dbh_cm)`` pairs.
        config: Uncertainty configuration. ``None`` uses defaults.

    Returns:
        :class:`UncertaintyResult` for the total project CO2e.

    Raises:
        ValueError: If ``trees`` is empty or contains invalid entries.

    Example:
        >>> r = propagate_population([("Acacia mangium", 20.0), ("Tectona grandis", 25.0)])
        >>> r.mean_tco2e > 0
        True
    """
    trees = list(trees)
    if not trees:
        raise ValueError("trees must contain at least one entry")
    cfg = config or UncertaintyConfig()
    total_samples = np.zeros(cfg.n_iterations)
    for idx, entry in enumerate(trees):
        if not isinstance(entry, tuple) or len(entry) != 2:
            raise ValueError(
                f"trees[{idx}] must be a (species, dbh_cm) tuple, got {entry!r}"
            )
        species, dbh_cm = entry
        # Give each tree a different RNG stream for independence.
        seed = None if cfg.random_seed is None else cfg.random_seed + idx
        per_tree_cfg = UncertaintyConfig(
            n_iterations=cfg.n_iterations,
            dbh_cv=cfg.dbh_cv,
            allometric_cv=cfg.allometric_cv,
            carbon_fraction_sd=cfg.carbon_fraction_sd,
            rs_ratio_cv=cfg.rs_ratio_cv,
            confidence_level=cfg.confidence_level,
            random_seed=seed,
        )
        tree_result = propagate_tree(species, dbh_cm, per_tree_cfg)
        total_samples += tree_result.samples
    return _summarize(total_samples, cfg)


def _summarize(samples: np.ndarray, cfg: UncertaintyConfig) -> UncertaintyResult:
    """Summarize Monte Carlo samples into an :class:`UncertaintyResult`."""
    alpha = 1.0 - cfg.confidence_level
    lower_pct = 100.0 * (alpha / 2.0)
    upper_pct = 100.0 * (1.0 - alpha / 2.0)
    ci_lower, ci_upper = np.percentile(samples, [lower_pct, upper_pct])
    return UncertaintyResult(
        mean_tco2e=float(samples.mean()),
        median_tco2e=float(np.median(samples)),
        std_tco2e=float(samples.std(ddof=1)),
        ci_lower=float(ci_lower),
        ci_upper=float(ci_upper),
        confidence_level=cfg.confidence_level,
        n_iterations=cfg.n_iterations,
        samples=samples,
    )
