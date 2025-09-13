"""
Allometric equations for biomass estimation.
Based on peer-reviewed equations for tropical species.
"""
import numpy as np

CARBON_FRACTION = 0.47
CO2_TO_C = 44 / 12

SPECIES_PARAMS = {
    "Acacia mangium":      {"a": 0.1277, "b": 2.3944, "rs": 0.26},
    "Tectona grandis":     {"a": 0.1532, "b": 2.3210, "rs": 0.25},
    "Eucalyptus pellita":  {"a": 0.1000, "b": 2.4500, "rs": 0.24},
    "Rhizophora apiculata":{"a": 0.1083, "b": 2.5530, "rs": 0.30},
    "Avicennia marina":    {"a": 0.0940, "b": 2.6100, "rs": 0.32},
    "generic_tropical":    {"a": 0.1184, "b": 2.3800, "rs": 0.26},
}

def estimate_agb(species: str, dbh_cm: float) -> float:
    """Above-ground biomass in kg using allometric equation."""
    p = SPECIES_PARAMS.get(species, SPECIES_PARAMS["generic_tropical"])
    return p["a"] * (dbh_cm ** p["b"])

def biomass_to_co2e(biomass_kg: float, rs_ratio: float = 0.26) -> float:
    """Convert total biomass to tCO2e."""
    total = biomass_kg * (1 + rs_ratio)
    carbon_kg = total * CARBON_FRACTION
    return carbon_kg * CO2_TO_C / 1000
