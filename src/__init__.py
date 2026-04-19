"""Carbon credit calculator package.

Public API:
    CarbonCreditCalculator - end-to-end credit pipeline.
    estimate_agb, estimate_agb_chave2014 - allometric equations.
    biomass_to_co2e, biomass_to_carbon, carbon_to_co2e - conversions.
    UncertaintyConfig, UncertaintyResult - Monte Carlo configuration.
    propagate_tree, propagate_population - uncertainty propagation.
"""
from .allometrics import (
    CARBON_FRACTION,
    CO2_TO_C,
    SPECIES_PARAMS,
    biomass_to_carbon,
    biomass_to_co2e,
    carbon_to_co2e,
    estimate_agb,
    estimate_agb_chave2014,
    get_species_params,
    tree_co2e,
)
from .main import CarbonCreditCalculator
from .uncertainty import (
    UncertaintyConfig,
    UncertaintyResult,
    propagate_population,
    propagate_tree,
)

__all__ = [
    "CARBON_FRACTION",
    "CO2_TO_C",
    "SPECIES_PARAMS",
    "CarbonCreditCalculator",
    "UncertaintyConfig",
    "UncertaintyResult",
    "biomass_to_carbon",
    "biomass_to_co2e",
    "carbon_to_co2e",
    "estimate_agb",
    "estimate_agb_chave2014",
    "get_species_params",
    "propagate_population",
    "propagate_tree",
    "tree_co2e",
]
