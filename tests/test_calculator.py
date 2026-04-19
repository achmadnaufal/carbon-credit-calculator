"""Tests for the CarbonCreditCalculator end-to-end pipeline."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.main import CarbonCreditCalculator
from src.uncertainty import UncertaintyConfig

DEMO_CSV = Path(__file__).resolve().parent.parent / "demo" / "sample_data.csv"


@pytest.fixture()
def sample_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "tree_id": ["T001", "T002", "T003"],
            "species": ["Acacia mangium", "Tectona grandis", "Swietenia macrophylla"],
            "dbh_cm": [20.0, 25.0, 18.0],
            "height_m": [14.0, 16.5, 12.5],
        }
    )


@pytest.fixture()
def calc() -> CarbonCreditCalculator:
    return CarbonCreditCalculator()


class TestValidation:
    def test_valid_df_passes(self, calc, sample_df):
        assert calc.validate(sample_df) is True

    def test_empty_df_raises(self, calc):
        with pytest.raises(ValueError, match="empty"):
            calc.validate(pd.DataFrame())

    def test_missing_required_column_raises(self, calc):
        with pytest.raises(ValueError, match="Missing required columns"):
            calc.validate(pd.DataFrame({"species": ["x"]}))

    def test_negative_dbh_raises(self, calc):
        df = pd.DataFrame({"species": ["x"], "dbh_cm": [-5.0]})
        with pytest.raises(ValueError, match="dbh_cm"):
            calc.validate(df)

    def test_nan_dbh_raises(self, calc):
        df = pd.DataFrame({"species": ["x"], "dbh_cm": [float("nan")]})
        with pytest.raises(ValueError):
            calc.validate(df)

    def test_missing_species_raises(self, calc):
        df = pd.DataFrame({"species": [None], "dbh_cm": [10.0]})
        with pytest.raises(ValueError):
            calc.validate(df)

    def test_negative_height_raises(self, calc):
        df = pd.DataFrame(
            {"species": ["Acacia mangium"], "dbh_cm": [10.0], "height_m": [-1.0]}
        )
        with pytest.raises(ValueError, match="height_m"):
            calc.validate(df)


class TestComputeTreeCo2e:
    def test_columns_added(self, calc, sample_df):
        out = calc.compute_tree_co2e(sample_df)
        for col in ("agb_kg", "carbon_kg", "co2e_tonnes"):
            assert col in out.columns

    def test_carbon_equals_biomass_times_047(self, calc, sample_df):
        out = calc.compute_tree_co2e(sample_df)
        # total biomass (AGB + roots) times 0.47 should equal carbon_kg.
        from src.allometrics import get_species_params

        for _, row in out.iterrows():
            rs = get_species_params(row["species"])["rs"]
            expected = row["agb_kg"] * (1 + rs) * 0.47
            assert row["carbon_kg"] == pytest.approx(expected, rel=1e-3)

    def test_co2e_equals_carbon_times_44_12_divided_1000(self, calc, sample_df):
        out = calc.compute_tree_co2e(sample_df)
        for _, row in out.iterrows():
            expected = row["carbon_kg"] * (44.0 / 12.0) / 1000.0
            assert row["co2e_tonnes"] == pytest.approx(expected, rel=1e-3)

    def test_chave_method_requires_height(self, sample_df):
        calc = CarbonCreditCalculator({"allometry_method": "chave2014"})
        df = sample_df.drop(columns=["height_m"])
        with pytest.raises(ValueError, match="height_m"):
            calc.compute_tree_co2e(df)

    def test_chave_method_runs_with_height(self, sample_df):
        calc = CarbonCreditCalculator({"allometry_method": "chave2014"})
        out = calc.compute_tree_co2e(sample_df)
        assert (out["co2e_tonnes"] > 0).all()


class TestComputeCredits:
    def test_totals_positive(self, calc, sample_df):
        result = calc.compute_credits(sample_df)
        assert result["total_trees"] == 3
        assert result["total_tco2e"] > 0
        assert result["net_credits_tco2e"] == result["total_tco2e"]  # no buffer

    def test_buffer_pool_reduces_net_credits(self, sample_df):
        calc = CarbonCreditCalculator({"buffer_pool_pct": 0.20})
        result = calc.compute_credits(sample_df)
        # Rounding to 6 decimals is applied on output; use abs tolerance.
        assert result["buffer_tco2e"] == pytest.approx(
            result["total_tco2e"] * 0.20, abs=1e-5
        )
        assert result["net_credits_tco2e"] == pytest.approx(
            result["total_tco2e"] - result["buffer_tco2e"], abs=1e-5
        )

    def test_bad_buffer_pct_raises(self, sample_df):
        calc = CarbonCreditCalculator({"buffer_pool_pct": 1.5})
        with pytest.raises(ValueError):
            calc.compute_credits(sample_df)

    def test_per_species_breakdown(self, calc, sample_df):
        result = calc.compute_credits(sample_df)
        species_in_breakdown = {row["species"] for row in result["per_species"]}
        assert species_in_breakdown == set(sample_df["species"])

    def test_uncertainty_block_included_when_requested(self, calc, sample_df):
        cfg = UncertaintyConfig(n_iterations=500, random_seed=1)
        result = calc.compute_credits(sample_df, uncertainty=True, uncertainty_config=cfg)
        assert "uncertainty" in result
        u = result["uncertainty"]
        assert u["ci_lower_tco2e"] < u["mean_tco2e"] < u["ci_upper_tco2e"]

    def test_empty_dataframe_raises(self, calc):
        with pytest.raises(ValueError):
            calc.compute_credits(pd.DataFrame())

    def test_does_not_mutate_input(self, calc, sample_df):
        before = sample_df.copy()
        calc.compute_credits(sample_df)
        pd.testing.assert_frame_equal(sample_df, before)


class TestLoadData:
    def test_loads_demo_csv(self, calc):
        df = calc.load_data(str(DEMO_CSV))
        assert len(df) >= 10
        assert "species" in df.columns
        assert "dbh_cm" in df.columns

    def test_missing_file_raises(self, calc):
        with pytest.raises(FileNotFoundError):
            calc.load_data("/nonexistent/path/data.csv")

    def test_unsupported_extension_raises(self, calc, tmp_path):
        bad = tmp_path / "foo.txt"
        bad.write_text("x")
        with pytest.raises(ValueError, match="Unsupported"):
            calc.load_data(str(bad))


class TestEndToEndWithDemoData:
    def test_demo_pipeline_produces_reasonable_credits(self, calc):
        result = calc.run(str(DEMO_CSV))
        # 20 trees, mixed ages — expect between 1 and 200 tCO2e.
        assert 1.0 < result["total_tco2e"] < 200.0
        assert result["total_trees"] == 20

    def test_demo_pipeline_with_uncertainty(self, calc):
        cfg = UncertaintyConfig(n_iterations=1000, random_seed=42)
        result = calc.run(str(DEMO_CSV), uncertainty=True, uncertainty_config=cfg)
        u = result["uncertainty"]
        assert u["ci_lower_tco2e"] > 0
        assert u["ci_upper_tco2e"] > u["ci_lower_tco2e"]
