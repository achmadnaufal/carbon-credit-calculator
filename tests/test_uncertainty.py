"""Tests for Monte Carlo uncertainty propagation."""
import numpy as np
import pytest

from src.uncertainty import (
    UncertaintyConfig,
    UncertaintyResult,
    propagate_population,
    propagate_tree,
)


class TestUncertaintyConfig:
    def test_defaults_are_valid(self):
        cfg = UncertaintyConfig()
        assert cfg.n_iterations >= 100
        assert 0 < cfg.confidence_level < 1

    def test_rejects_too_few_iterations(self):
        with pytest.raises(ValueError):
            UncertaintyConfig(n_iterations=10)

    def test_rejects_bad_confidence_level(self):
        with pytest.raises(ValueError):
            UncertaintyConfig(confidence_level=1.5)

    def test_rejects_negative_cv(self):
        with pytest.raises(ValueError):
            UncertaintyConfig(dbh_cv=-0.1)


class TestPropagateTree:
    def test_result_fields_populated(self):
        result = propagate_tree("Acacia mangium", 20.0)
        assert isinstance(result, UncertaintyResult)
        assert result.mean_tco2e > 0
        assert result.ci_lower < result.mean_tco2e < result.ci_upper
        assert result.std_tco2e > 0
        assert result.n_iterations == 10_000

    def test_mean_roughly_matches_point_estimate(self):
        # Monte Carlo mean should be within ~20% of a deterministic point.
        from src.allometrics import tree_co2e

        point = tree_co2e("Tectona grandis", 25.0)
        mc = propagate_tree("Tectona grandis", 25.0)
        assert mc.mean_tco2e == pytest.approx(point, rel=0.25)

    def test_reproducible_with_seed(self):
        cfg = UncertaintyConfig(random_seed=123)
        r1 = propagate_tree("Acacia mangium", 20.0, cfg)
        r2 = propagate_tree("Acacia mangium", 20.0, cfg)
        assert r1.mean_tco2e == r2.mean_tco2e
        assert r1.ci_upper == r2.ci_upper

    def test_higher_cv_widens_confidence_interval(self):
        narrow = propagate_tree(
            "Acacia mangium", 20.0,
            UncertaintyConfig(allometric_cv=0.05, random_seed=7),
        )
        wide = propagate_tree(
            "Acacia mangium", 20.0,
            UncertaintyConfig(allometric_cv=0.40, random_seed=7),
        )
        narrow_width = narrow.ci_upper - narrow.ci_lower
        wide_width = wide.ci_upper - wide.ci_lower
        assert wide_width > narrow_width

    def test_samples_are_all_positive(self):
        result = propagate_tree("Acacia mangium", 20.0)
        assert np.all(result.samples >= 0)

    def test_as_dict_has_expected_keys(self):
        d = propagate_tree("Acacia mangium", 20.0).as_dict()
        for key in (
            "mean_tco2e",
            "median_tco2e",
            "std_tco2e",
            "ci_lower_tco2e",
            "ci_upper_tco2e",
            "confidence_level",
            "n_iterations",
        ):
            assert key in d

    def test_rejects_negative_dbh(self):
        with pytest.raises(ValueError):
            propagate_tree("Acacia mangium", -5.0)

    def test_rejects_nan_dbh(self):
        with pytest.raises(ValueError):
            propagate_tree("Acacia mangium", float("nan"))


class TestPropagatePopulation:
    def test_population_total_greater_than_single_tree(self):
        trees = [
            ("Acacia mangium", 20.0),
            ("Tectona grandis", 25.0),
            ("Swietenia macrophylla", 18.0),
        ]
        single = propagate_tree("Acacia mangium", 20.0)
        pop = propagate_population(trees)
        assert pop.mean_tco2e > single.mean_tco2e

    def test_population_approximates_sum_of_means(self):
        trees = [("Acacia mangium", 20.0)] * 5
        individual_sum = 5 * propagate_tree(
            "Acacia mangium", 20.0, UncertaintyConfig(random_seed=1)
        ).mean_tco2e
        pop = propagate_population(trees)
        # Within 10% of the deterministic-sum approximation.
        assert pop.mean_tco2e == pytest.approx(individual_sum, rel=0.10)

    def test_empty_population_raises(self):
        with pytest.raises(ValueError):
            propagate_population([])

    def test_invalid_tuple_raises(self):
        with pytest.raises(ValueError):
            propagate_population([("Acacia mangium",)])  # type: ignore[list-item]
