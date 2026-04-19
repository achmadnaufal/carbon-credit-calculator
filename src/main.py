"""Carbon credit calculator for tree-planting projects.

Loads tree inventory data (CSV or Excel), validates it, computes
per-tree and project-level CO2e sequestration with optional Monte
Carlo uncertainty, and reports tradable carbon credits.

Author: github.com/achmadnaufal
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

from .allometrics import (
    CARBON_FRACTION,
    CO2_TO_C,
    biomass_to_co2e,
    estimate_agb,
    estimate_agb_chave2014,
    get_species_params,
)
from .uncertainty import (
    UncertaintyConfig,
    UncertaintyResult,
    propagate_population,
)

REQUIRED_COLUMNS: tuple = ("species", "dbh_cm")
OPTIONAL_COLUMNS: tuple = (
    "tree_id",
    "height_m",
    "age_years",
    "plot_id",
    "date_measured",
)


class CarbonCreditCalculator:
    """Estimate carbon credits from a tree-inventory DataFrame.

    One tCO2e of net sequestration equals one tradable carbon credit
    under most voluntary standards (VCS/Verra, Gold Standard). Results
    are gross sequestration; buyers should deduct buffer pool
    contributions and baselines before trading.

    Attributes:
        config: Optional configuration dict. Recognised keys:

            * ``allometry_method`` (``"power"`` or ``"chave2014"``,
              default ``"power"``).
            * ``carbon_fraction`` (default 0.47).
            * ``buffer_pool_pct`` (fraction held in VCS buffer pool,
              default 0.0).

    Example:
        >>> calc = CarbonCreditCalculator()
        >>> df = pd.DataFrame({
        ...     "species": ["Acacia mangium"],
        ...     "dbh_cm": [20.0],
        ...     "height_m": [12.0],
        ... })
        >>> result = calc.compute_credits(df)
        >>> result["total_tco2e"] > 0
        True
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None) -> None:
        self.config: Dict[str, Any] = dict(config) if config else {}

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    def load_data(self, filepath: str) -> pd.DataFrame:
        """Load an inventory from CSV or Excel.

        Args:
            filepath: Path to a ``.csv``, ``.xlsx``, or ``.xls`` file.

        Returns:
            A pandas DataFrame with the raw data.

        Raises:
            FileNotFoundError: If ``filepath`` does not exist.
            ValueError: If ``filepath`` has an unsupported extension or
                the file is empty.

        Example:
            >>> calc = CarbonCreditCalculator()
            >>> df = calc.load_data("demo/sample_data.csv")  # doctest: +SKIP
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        suffix = path.suffix.lower()
        if suffix in (".xlsx", ".xls"):
            df = pd.read_excel(path)
        elif suffix == ".csv":
            df = pd.read_csv(path)
        else:
            raise ValueError(
                f"Unsupported file type {suffix!r}; use .csv, .xlsx, or .xls"
            )
        if df.empty:
            raise ValueError(f"Loaded file is empty: {filepath}")
        return df

    # ------------------------------------------------------------------
    # Validation & preprocessing
    # ------------------------------------------------------------------
    def validate(self, df: pd.DataFrame) -> bool:
        """Validate an inventory DataFrame for the credit pipeline.

        Checks schema, emptiness, numeric types, and value ranges.

        Args:
            df: Tree-inventory DataFrame.

        Returns:
            ``True`` if the DataFrame passes validation.

        Raises:
            ValueError: If the DataFrame is empty, missing required
                columns, or contains non-positive / non-finite numeric
                values in ``dbh_cm`` or ``height_m``.

        Example:
            >>> calc = CarbonCreditCalculator()
            >>> calc.validate(pd.DataFrame({"species": ["x"], "dbh_cm": [10]}))
            True
        """
        if df is None:
            raise ValueError("DataFrame must not be None")
        if df.empty:
            raise ValueError("Input DataFrame is empty")
        normalized = {c.lower().strip() for c in df.columns}
        missing = [c for c in REQUIRED_COLUMNS if c not in normalized]
        if missing:
            raise ValueError(
                f"Missing required columns: {missing}. "
                f"Required: {list(REQUIRED_COLUMNS)}"
            )
        dbh = pd.to_numeric(df["dbh_cm"], errors="coerce")
        if dbh.isna().any():
            raise ValueError("dbh_cm contains non-numeric or missing values")
        if (dbh <= 0).any():
            raise ValueError("dbh_cm must be > 0 for every row")
        if "height_m" in df.columns:
            h = pd.to_numeric(df["height_m"], errors="coerce")
            if h.isna().any():
                raise ValueError("height_m contains non-numeric or missing values")
            if (h <= 0).any():
                raise ValueError("height_m must be > 0 for every row")
        if df["species"].isna().any():
            raise ValueError("species contains missing values")
        return True

    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalize column names and drop blank rows.

        Returns a new DataFrame; the input is never mutated.

        Args:
            df: Raw inventory DataFrame.

        Returns:
            A new DataFrame with lower-case column names and no fully
            empty rows.
        """
        cleaned = df.copy()
        cleaned.dropna(how="all", inplace=True)
        cleaned.columns = [c.lower().strip().replace(" ", "_") for c in cleaned.columns]
        return cleaned

    # ------------------------------------------------------------------
    # Core carbon calculations
    # ------------------------------------------------------------------
    def compute_tree_co2e(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute per-tree AGB, carbon, and CO2e.

        Adds the columns ``agb_kg``, ``carbon_kg``, and ``co2e_tonnes``
        to a copy of ``df``.

        Args:
            df: Inventory DataFrame (already preprocessed).

        Returns:
            A new DataFrame with the three computed columns appended.

        Raises:
            ValueError: If validation fails or the configured allometry
                method requires a column that is missing.
        """
        self.validate(df)
        method = self.config.get("allometry_method", "power")
        carbon_fraction = float(self.config.get("carbon_fraction", CARBON_FRACTION))
        if not (0 < carbon_fraction <= 1):
            raise ValueError(
                f"config.carbon_fraction must be in (0, 1], got {carbon_fraction}"
            )
        if method == "chave2014" and "height_m" not in df.columns:
            raise ValueError("allometry_method='chave2014' requires a height_m column")

        out = df.copy()
        agb_values: List[float] = []
        rs_values: List[float] = []
        for _, row in out.iterrows():
            species = str(row["species"])
            dbh = float(row["dbh_cm"])
            params = get_species_params(species)
            if method == "power":
                agb = estimate_agb(species, dbh)
            elif method == "chave2014":
                agb = estimate_agb_chave2014(dbh, float(row["height_m"]), params["wd"])
            else:
                raise ValueError(
                    f"Unknown allometry_method {method!r}; "
                    f"use 'power' or 'chave2014'"
                )
            agb_values.append(agb)
            rs_values.append(params["rs"])
        agb_arr = np.asarray(agb_values)
        rs_arr = np.asarray(rs_values)
        total_biomass = agb_arr * (1.0 + rs_arr)
        carbon_kg = total_biomass * carbon_fraction
        co2e_tonnes = carbon_kg * CO2_TO_C / 1000.0
        out["agb_kg"] = np.round(agb_arr, 4)
        out["carbon_kg"] = np.round(carbon_kg, 4)
        out["co2e_tonnes"] = np.round(co2e_tonnes, 6)
        return out

    def compute_credits(
        self,
        df: pd.DataFrame,
        uncertainty: bool = False,
        uncertainty_config: Optional[UncertaintyConfig] = None,
    ) -> Dict[str, Any]:
        """Compute project-level carbon credit totals.

        Args:
            df: Tree inventory DataFrame.
            uncertainty: If ``True``, run Monte Carlo propagation to
                obtain a confidence interval for total tCO2e.
            uncertainty_config: Optional :class:`UncertaintyConfig`
                overriding Monte Carlo defaults.

        Returns:
            Dict with keys: ``total_trees``, ``total_agb_kg``,
            ``total_carbon_kg``, ``total_tco2e``, ``buffer_tco2e``,
            ``net_credits_tco2e``, ``per_species``, and optionally
            ``uncertainty`` (an :class:`UncertaintyResult` dict).

        Raises:
            ValueError: If validation fails.

        Example:
            >>> calc = CarbonCreditCalculator({"buffer_pool_pct": 0.2})
            >>> df = pd.DataFrame({"species": ["Acacia mangium"], "dbh_cm": [20.0]})
            >>> r = calc.compute_credits(df)
            >>> round(r["net_credits_tco2e"], 4)
            0.2892
        """
        clean = self.preprocess(df)
        per_tree = self.compute_tree_co2e(clean)
        total_tco2e = float(per_tree["co2e_tonnes"].sum())
        buffer_pct = float(self.config.get("buffer_pool_pct", 0.0))
        if not (0.0 <= buffer_pct <= 1.0):
            raise ValueError(
                f"config.buffer_pool_pct must be in [0, 1], got {buffer_pct}"
            )
        buffer_tco2e = total_tco2e * buffer_pct
        net_credits = total_tco2e - buffer_tco2e

        per_species = (
            per_tree.groupby("species")["co2e_tonnes"]
            .agg(["count", "sum", "mean"])
            .rename(columns={"count": "n_trees", "sum": "tco2e", "mean": "tco2e_per_tree"})
            .round(6)
            .reset_index()
            .to_dict(orient="records")
        )

        result: Dict[str, Any] = {
            "total_trees": int(len(per_tree)),
            "total_agb_kg": round(float(per_tree["agb_kg"].sum()), 2),
            "total_carbon_kg": round(float(per_tree["carbon_kg"].sum()), 2),
            "total_tco2e": round(total_tco2e, 6),
            "buffer_pool_pct": buffer_pct,
            "buffer_tco2e": round(buffer_tco2e, 6),
            "net_credits_tco2e": round(net_credits, 6),
            "per_species": per_species,
        }

        if uncertainty:
            trees: Iterable = (
                (str(row["species"]), float(row["dbh_cm"]))
                for _, row in clean.iterrows()
            )
            mc: UncertaintyResult = propagate_population(
                list(trees), config=uncertainty_config
            )
            result["uncertainty"] = mc.as_dict()
        return result

    # ------------------------------------------------------------------
    # Legacy compatibility
    # ------------------------------------------------------------------
    def analyze(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Carbon-aware summary (replaces the legacy generic analyze).

        Returns the same dict as :meth:`compute_credits` without
        uncertainty. Kept for backward compatibility with earlier
        versions of this tool.

        Args:
            df: Tree inventory DataFrame.

        Returns:
            Summary dict with per-species and project totals.
        """
        return self.compute_credits(df, uncertainty=False)

    def run(
        self,
        filepath: str,
        uncertainty: bool = False,
        uncertainty_config: Optional[UncertaintyConfig] = None,
    ) -> Dict[str, Any]:
        """End-to-end pipeline: load -> validate -> compute credits.

        Args:
            filepath: Path to an inventory CSV/Excel file.
            uncertainty: Whether to include Monte Carlo propagation.
            uncertainty_config: Optional override for MC configuration.

        Returns:
            The dict returned by :meth:`compute_credits`.
        """
        df = self.load_data(filepath)
        return self.compute_credits(
            df, uncertainty=uncertainty, uncertainty_config=uncertainty_config
        )

    def to_dataframe(self, result: Dict[str, Any]) -> pd.DataFrame:
        """Flatten a credit-result dict to a long-format DataFrame.

        Args:
            result: Dict as returned by :meth:`compute_credits`.

        Returns:
            DataFrame with columns ``metric`` and ``value``; nested
            dict values are expanded with dotted keys.
        """
        rows: List[Dict[str, Any]] = []
        for key, value in result.items():
            if isinstance(value, dict):
                for sub_key, sub_val in value.items():
                    rows.append({"metric": f"{key}.{sub_key}", "value": sub_val})
            elif isinstance(value, list):
                rows.append({"metric": key, "value": f"<{len(value)} records>"})
            else:
                rows.append({"metric": key, "value": value})
        return pd.DataFrame(rows)
