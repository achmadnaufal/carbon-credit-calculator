"""Allometric equations for tropical tree biomass estimation.

Implements species-specific power-law allometry of the form
``AGB = a * DBH^b`` for above-ground biomass, plus carbon and
CO2-equivalent conversions.

References:
    IPCC (2006) Guidelines for National Greenhouse Gas Inventories,
        Volume 4: Agriculture, Forestry and Other Land Use.
    Chave, J. et al. (2014) Improved allometric models to estimate the
        aboveground biomass of tropical trees. Global Change Biology 20,
        3177-3190. (Pan-tropical moist forest model, used for the
        ``chave2014`` method.)
    Mokany, K., Raison, R. J., Prokushkin, A. S. (2006) Critical analysis
        of root:shoot ratios in terrestrial biomes.
"""
from __future__ import annotations

import math
from typing import Optional

# IPCC default carbon fraction of dry biomass (Table 4.3, IPCC 2006).
CARBON_FRACTION: float = 0.47

# Molecular mass ratio CO2 : C = 44 / 12.
CO2_TO_C: float = 44.0 / 12.0

# Species-specific parameters for AGB = a * DBH^b (kg, cm).
# ``rs`` is the root-to-shoot biomass ratio. ``wd`` is wood density
# (g/cm^3, oven-dry mass / green volume) used by the Chave 2014 model.
SPECIES_PARAMS: dict = {
    "Acacia mangium":          {"a": 0.1277, "b": 2.3944, "rs": 0.26, "wd": 0.52},
    "Tectona grandis":         {"a": 0.1532, "b": 2.3210, "rs": 0.25, "wd": 0.66},
    "Eucalyptus pellita":      {"a": 0.1000, "b": 2.4500, "rs": 0.24, "wd": 0.78},
    "Eucalyptus urophylla":    {"a": 0.0998, "b": 2.4450, "rs": 0.24, "wd": 0.72},
    "Swietenia macrophylla":   {"a": 0.1340, "b": 2.3700, "rs": 0.26, "wd": 0.55},
    "Gliricidia sepium":       {"a": 0.1200, "b": 2.3500, "rs": 0.27, "wd": 0.65},
    "Paraserianthes falcataria": {"a": 0.1100, "b": 2.3600, "rs": 0.26, "wd": 0.33},
    "Rhizophora apiculata":    {"a": 0.1083, "b": 2.5530, "rs": 0.30, "wd": 0.87},
    "Avicennia marina":        {"a": 0.0940, "b": 2.6100, "rs": 0.32, "wd": 0.67},
    "generic_tropical":        {"a": 0.1184, "b": 2.3800, "rs": 0.26, "wd": 0.60},
}


def _require_positive(name: str, value: float) -> None:
    """Raise ValueError if ``value`` is not a positive finite number.

    Args:
        name: Parameter name used in the error message.
        value: Value to validate.

    Raises:
        ValueError: If ``value`` is ``None``, NaN, non-numeric, or not > 0.
    """
    if value is None:
        raise ValueError(f"{name} must not be None")
    try:
        v = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric, got {value!r}") from exc
    if math.isnan(v) or math.isinf(v):
        raise ValueError(f"{name} must be a finite number, got {value!r}")
    if v <= 0:
        raise ValueError(f"{name} must be > 0, got {v}")


def get_species_params(species: str) -> dict:
    """Return allometric parameters for a species, with fallback.

    Unknown species fall back to ``generic_tropical`` parameters. Species
    names are matched case-insensitively against the known keys.

    Args:
        species: Scientific name (e.g. ``"Tectona grandis"``).

    Returns:
        Dict with keys ``a``, ``b``, ``rs``, ``wd``.

    Raises:
        ValueError: If ``species`` is empty or not a string.

    Example:
        >>> get_species_params("Tectona grandis")["b"]
        2.321
    """
    if not isinstance(species, str) or not species.strip():
        raise ValueError("species must be a non-empty string")
    if species in SPECIES_PARAMS:
        return SPECIES_PARAMS[species]
    lowered = species.strip().lower()
    for key, params in SPECIES_PARAMS.items():
        if key.lower() == lowered:
            return params
    return SPECIES_PARAMS["generic_tropical"]


def estimate_agb(species: str, dbh_cm: float) -> float:
    """Estimate above-ground biomass using species-specific allometry.

    Uses the power-law form ``AGB = a * DBH^b`` with species-specific
    coefficients. Falls back to generic tropical parameters when the
    species is not in :data:`SPECIES_PARAMS`.

    Args:
        species: Scientific name of the tree.
        dbh_cm: Diameter at breast height (1.3 m) in centimetres.

    Returns:
        Above-ground biomass in kilograms (dry weight).

    Raises:
        ValueError: If ``dbh_cm`` is not a positive finite number.

    Example:
        >>> round(estimate_agb("Acacia mangium", 20.0), 2)
        166.49
    """
    _require_positive("dbh_cm", dbh_cm)
    params = get_species_params(species if isinstance(species, str) else "generic_tropical")
    return float(params["a"] * (float(dbh_cm) ** params["b"]))


def estimate_agb_chave2014(
    dbh_cm: float,
    height_m: float,
    wood_density: float,
) -> float:
    """Estimate AGB using the Chave et al. (2014) pan-tropical model.

    Applies the form ``AGB = 0.0673 * (rho * D^2 * H)^0.976`` where
    ``rho`` is wood density (g/cm^3), ``D`` is DBH (cm) and ``H`` is
    total tree height (m). Returns biomass in kilograms.

    Args:
        dbh_cm: Diameter at breast height in centimetres.
        height_m: Total tree height in metres.
        wood_density: Wood density in g/cm^3 (typical range 0.2-1.0).

    Returns:
        Above-ground biomass in kilograms (dry weight).

    Raises:
        ValueError: If any input is not a positive finite number, or if
            ``wood_density`` is implausibly large (> 1.5 g/cm^3).

    Example:
        >>> round(estimate_agb_chave2014(20.0, 15.0, 0.6), 2)
        199.05
    """
    _require_positive("dbh_cm", dbh_cm)
    _require_positive("height_m", height_m)
    _require_positive("wood_density", wood_density)
    if wood_density > 1.5:
        raise ValueError(
            f"wood_density must be <= 1.5 g/cm^3, got {wood_density}"
        )
    return float(
        0.0673
        * ((float(wood_density) * (float(dbh_cm) ** 2) * float(height_m)) ** 0.976)
    )


def biomass_to_carbon(biomass_kg: float, carbon_fraction: float = CARBON_FRACTION) -> float:
    """Convert dry biomass to elemental carbon mass.

    Args:
        biomass_kg: Dry biomass in kilograms.
        carbon_fraction: Fraction of biomass that is carbon (default 0.47).

    Returns:
        Carbon mass in kilograms.

    Raises:
        ValueError: If ``biomass_kg`` is negative or not finite, or if
            ``carbon_fraction`` is outside (0, 1].

    Example:
        >>> round(biomass_to_carbon(100.0), 2)
        47.0
    """
    if biomass_kg is None or not math.isfinite(float(biomass_kg)):
        raise ValueError(f"biomass_kg must be finite, got {biomass_kg!r}")
    if float(biomass_kg) < 0:
        raise ValueError(f"biomass_kg must be >= 0, got {biomass_kg}")
    if not (0 < float(carbon_fraction) <= 1):
        raise ValueError(
            f"carbon_fraction must be in (0, 1], got {carbon_fraction}"
        )
    return float(biomass_kg) * float(carbon_fraction)


def carbon_to_co2e(carbon_kg: float) -> float:
    """Convert carbon mass to CO2-equivalent using the 44/12 stoichiometric ratio.

    Args:
        carbon_kg: Elemental carbon mass in kilograms.

    Returns:
        CO2-equivalent mass in kilograms.

    Raises:
        ValueError: If ``carbon_kg`` is negative or not finite.

    Example:
        >>> round(carbon_to_co2e(12.0), 4)
        44.0
    """
    if carbon_kg is None or not math.isfinite(float(carbon_kg)):
        raise ValueError(f"carbon_kg must be finite, got {carbon_kg!r}")
    if float(carbon_kg) < 0:
        raise ValueError(f"carbon_kg must be >= 0, got {carbon_kg}")
    return float(carbon_kg) * CO2_TO_C


def biomass_to_co2e(
    biomass_kg: float,
    rs_ratio: float = 0.26,
    carbon_fraction: float = CARBON_FRACTION,
) -> float:
    """Convert AGB to total CO2-equivalent (tonnes) including roots.

    Applies the root-to-shoot ratio to account for below-ground biomass,
    then converts to carbon and finally to CO2e in tonnes.

    Args:
        biomass_kg: Above-ground biomass in kilograms.
        rs_ratio: Root-to-shoot biomass ratio (default 0.26).
        carbon_fraction: Carbon fraction of dry biomass (default 0.47).

    Returns:
        Total sequestered CO2e in tonnes (tCO2e).

    Raises:
        ValueError: If ``biomass_kg`` is negative, ``rs_ratio`` is
            negative, or inputs are not finite.

    Example:
        >>> round(biomass_to_co2e(1000.0), 4)
        2.1714
    """
    if biomass_kg is None or not math.isfinite(float(biomass_kg)):
        raise ValueError(f"biomass_kg must be finite, got {biomass_kg!r}")
    if float(biomass_kg) < 0:
        raise ValueError(f"biomass_kg must be >= 0, got {biomass_kg}")
    if rs_ratio is None or not math.isfinite(float(rs_ratio)) or float(rs_ratio) < 0:
        raise ValueError(f"rs_ratio must be >= 0 and finite, got {rs_ratio!r}")
    total_biomass = float(biomass_kg) * (1.0 + float(rs_ratio))
    carbon_kg = biomass_to_carbon(total_biomass, carbon_fraction)
    return carbon_to_co2e(carbon_kg) / 1000.0


def tree_co2e(
    species: str,
    dbh_cm: float,
    height_m: Optional[float] = None,
    method: str = "power",
) -> float:
    """Compute total CO2e for a single tree, roots + shoots.

    Args:
        species: Scientific name (unknown species -> generic tropical).
        dbh_cm: DBH in centimetres.
        height_m: Height in metres (required for ``method="chave2014"``).
        method: One of ``"power"`` (species-specific power law) or
            ``"chave2014"`` (pan-tropical model using wood density).

    Returns:
        Total CO2e in tonnes (tCO2e) for this tree.

    Raises:
        ValueError: If inputs are invalid or ``method`` is unsupported.

    Example:
        >>> round(tree_co2e("Acacia mangium", 20.0), 4)
        0.3615
    """
    params = get_species_params(species if isinstance(species, str) else "generic_tropical")
    if method == "power":
        agb = estimate_agb(species, dbh_cm)
    elif method == "chave2014":
        if height_m is None:
            raise ValueError("height_m is required for method='chave2014'")
        agb = estimate_agb_chave2014(dbh_cm, height_m, params["wd"])
    else:
        raise ValueError(f"Unknown method {method!r}; use 'power' or 'chave2014'")
    return biomass_to_co2e(agb, rs_ratio=params["rs"])
