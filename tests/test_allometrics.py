"""Unit tests for the allometrics module."""
import math

import pytest

from src.allometrics import (
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


class TestConstants:
    def test_carbon_fraction_ipcc_default(self):
        assert CARBON_FRACTION == pytest.approx(0.47)

    def test_co2_c_stoichiometric_ratio(self):
        assert CO2_TO_C == pytest.approx(44.0 / 12.0)


class TestSpeciesParams:
    def test_known_species_returns_exact_params(self):
        p = get_species_params("Tectona grandis")
        assert p["a"] == pytest.approx(0.1532)
        assert p["b"] == pytest.approx(2.3210)

    def test_unknown_species_falls_back_to_generic(self):
        p = get_species_params("Unobtainium giganteum")
        assert p == SPECIES_PARAMS["generic_tropical"]

    def test_case_insensitive_lookup(self):
        p = get_species_params("tectona grandis")
        assert p["b"] == pytest.approx(2.3210)

    def test_empty_species_raises(self):
        with pytest.raises(ValueError):
            get_species_params("")

    def test_none_species_raises(self):
        with pytest.raises(ValueError):
            get_species_params(None)  # type: ignore[arg-type]


class TestEstimateAgb:
    def test_agb_known_species_within_sensible_range(self):
        # 20 cm Acacia mangium: expect ~100-200 kg dry biomass.
        agb = estimate_agb("Acacia mangium", 20.0)
        assert 100 < agb < 200

    def test_agb_scales_monotonically_with_dbh(self):
        assert (
            estimate_agb("Tectona grandis", 10.0)
            < estimate_agb("Tectona grandis", 20.0)
            < estimate_agb("Tectona grandis", 40.0)
        )

    def test_agb_mature_tree_under_5_tonnes(self):
        # 60 cm DBH tropical tree should be < 5 tonnes dry biomass.
        agb = estimate_agb("generic_tropical", 60.0)
        assert agb < 5000

    def test_negative_dbh_raises(self):
        with pytest.raises(ValueError):
            estimate_agb("Tectona grandis", -10.0)

    def test_zero_dbh_raises(self):
        with pytest.raises(ValueError):
            estimate_agb("Tectona grandis", 0.0)

    def test_nan_dbh_raises(self):
        with pytest.raises(ValueError):
            estimate_agb("Tectona grandis", float("nan"))


class TestChave2014:
    def test_chave_with_realistic_inputs(self):
        # 30 cm DBH, 20 m tall, wd=0.6 should be a few hundred kg.
        agb = estimate_agb_chave2014(30.0, 20.0, 0.6)
        assert 200 < agb < 800

    def test_chave_scales_with_wood_density(self):
        light = estimate_agb_chave2014(20.0, 15.0, 0.3)
        heavy = estimate_agb_chave2014(20.0, 15.0, 0.9)
        assert heavy > light

    def test_chave_rejects_negative_height(self):
        with pytest.raises(ValueError):
            estimate_agb_chave2014(20.0, -5.0, 0.6)

    def test_chave_rejects_implausible_wood_density(self):
        with pytest.raises(ValueError):
            estimate_agb_chave2014(20.0, 15.0, 2.0)


class TestConversions:
    def test_biomass_to_carbon_exact(self):
        assert biomass_to_carbon(100.0) == pytest.approx(47.0)

    def test_biomass_to_carbon_custom_fraction(self):
        assert biomass_to_carbon(100.0, 0.50) == pytest.approx(50.0)

    def test_biomass_to_carbon_zero_is_zero(self):
        assert biomass_to_carbon(0.0) == 0.0

    def test_biomass_to_carbon_negative_raises(self):
        with pytest.raises(ValueError):
            biomass_to_carbon(-1.0)

    def test_biomass_to_carbon_invalid_fraction_raises(self):
        with pytest.raises(ValueError):
            biomass_to_carbon(100.0, 1.5)
        with pytest.raises(ValueError):
            biomass_to_carbon(100.0, 0.0)

    def test_carbon_to_co2e_exact_44_12(self):
        # 12 kg C -> 44 kg CO2e exactly.
        assert carbon_to_co2e(12.0) == pytest.approx(44.0)

    def test_carbon_to_co2e_negative_raises(self):
        with pytest.raises(ValueError):
            carbon_to_co2e(-1.0)

    def test_biomass_to_co2e_end_to_end(self):
        # Hand calculation: 1000 kg AGB * 1.26 * 0.47 * 44/12 / 1000
        expected = 1000.0 * 1.26 * 0.47 * (44.0 / 12.0) / 1000.0
        assert biomass_to_co2e(1000.0, 0.26) == pytest.approx(expected)

    def test_biomass_to_co2e_negative_rs_raises(self):
        with pytest.raises(ValueError):
            biomass_to_co2e(100.0, rs_ratio=-0.1)


class TestTreeCo2e:
    def test_power_method(self):
        tco2e = tree_co2e("Acacia mangium", 20.0)
        assert 0.2 < tco2e < 0.5

    def test_chave_method_requires_height(self):
        with pytest.raises(ValueError):
            tree_co2e("Acacia mangium", 20.0, method="chave2014")

    def test_chave_method_runs_with_height(self):
        tco2e = tree_co2e("Acacia mangium", 20.0, height_m=15.0, method="chave2014")
        assert tco2e > 0

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError):
            tree_co2e("Acacia mangium", 20.0, method="bogus")
