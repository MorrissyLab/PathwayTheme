"""Roll pathway-level results up into broad categories.

Given a long result table (a :class:`~pathwaytheme.contracts.DiffResult` table,
or PCA top-pathways/loadings) and a term -> category mapping, aggregate a numeric
column per category.  This is the downstream summarization step from the ASPS
workflow (e.g. mean t-statistic per GO category, significant vs not).
"""

from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd


def _as_mapping(mapping: Union[pd.Series, dict]) -> dict:
    if isinstance(mapping, pd.Series):
        return {str(k): v for k, v in mapping.items()}
    return {str(k): v for k, v in dict(mapping).items()}


def summarize_by_category(table: pd.DataFrame,
                          mapping: Union[pd.Series, dict], *,
                          key_col: str = "pathway",
                          value_col: str = "effect",
                          group_cols: Optional[list[str]] = None,
                          stat: str = "mean",
                          unmapped_label: str = "Other",
                          significance_col: Optional[str] = None,
                          alpha: float = 0.05) -> pd.DataFrame:
    """Aggregate ``value_col`` per category (optionally within ``group_cols``).

    Returns a table with one row per (group..., [significance,] category):
    ``<stat>_<value_col>``, ``n_pathways``, ``n_up``, ``n_down``.

    If ``significance_col`` is given (e.g. ``"fdr"``), each category is split
    into ``significant`` (value < ``alpha``) vs ``non_significant`` rows, adding
    a ``significance`` column — the sig-vs-nonsig-per-category summary from the
    reference ASPS GO-slim analysis.  NaN thresholds count as non-significant.
    """
    if key_col not in table.columns:
        raise KeyError(f"key_col {key_col!r} not in table columns {list(table.columns)}")
    if value_col not in table.columns:
        raise KeyError(f"value_col {value_col!r} not in table columns {list(table.columns)}")
    if stat not in ("mean", "median", "sum"):
        raise ValueError(f"stat must be mean|median|sum, got {stat!r}")

    m = _as_mapping(mapping)
    df = table.copy()
    df["category"] = df[key_col].astype(str).map(m).fillna(unmapped_label)

    keys = [c for c in (group_cols or []) if c in df.columns]
    if significance_col is not None:
        if significance_col not in df.columns:
            raise KeyError(f"significance_col {significance_col!r} not in table "
                           f"columns {list(table.columns)}")
        is_sig = df[significance_col].to_numpy(dtype=float) < alpha  # NaN -> False
        df["significance"] = np.where(is_sig, "significant", "non_significant")
        keys = keys + ["significance"]
    keys = keys + ["category"]

    df["_is_up"] = (df[value_col] > 0).astype(int)
    df["_is_down"] = (df[value_col] < 0).astype(int)

    summary = df.groupby(keys, dropna=False).agg(
        stat_value=(value_col, stat),
        n_pathways=(value_col, "size"),
        n_up=("_is_up", "sum"),
        n_down=("_is_down", "sum"),
    ).reset_index().rename(columns={"stat_value": f"{stat}_{value_col}"})
    return summary.sort_values(keys).reset_index(drop=True)
