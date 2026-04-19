"""Basic usage: compute carbon credits from the bundled demo data."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src import CarbonCreditCalculator, UncertaintyConfig

DEMO = Path(__file__).parent.parent / "demo" / "sample_data.csv"


def main() -> None:
    calc = CarbonCreditCalculator({"buffer_pool_pct": 0.20})
    cfg = UncertaintyConfig(n_iterations=5000, random_seed=42)
    result = calc.run(str(DEMO), uncertainty=True, uncertainty_config=cfg)

    print("=== Carbon Credit Project Summary ===")
    print(f"Trees inventoried:    {result['total_trees']}")
    print(f"Total AGB (kg):       {result['total_agb_kg']:,.2f}")
    print(f"Total carbon (kg):    {result['total_carbon_kg']:,.2f}")
    print(f"Gross tCO2e:          {result['total_tco2e']:,.4f}")
    print(f"Buffer pool (20%):    {result['buffer_tco2e']:,.4f}")
    print(f"Net tradable credits: {result['net_credits_tco2e']:,.4f}")

    if "uncertainty" in result:
        u = result["uncertainty"]
        print("\n=== Monte Carlo Uncertainty (90% CI) ===")
        print(f"Mean tCO2e:   {u['mean_tco2e']:,.4f}")
        print(f"Median tCO2e: {u['median_tco2e']:,.4f}")
        print(f"Std. dev.:    {u['std_tco2e']:,.4f}")
        print(f"CI:           [{u['ci_lower_tco2e']:,.4f}, {u['ci_upper_tco2e']:,.4f}]")

    print("\n=== Per-species breakdown ===")
    for row in result["per_species"]:
        print(
            f"  {row['species']:<28} n={row['n_trees']:>3}  "
            f"tCO2e={row['tco2e']:>8.4f}  per_tree={row['tco2e_per_tree']:.4f}"
        )


if __name__ == "__main__":
    main()
