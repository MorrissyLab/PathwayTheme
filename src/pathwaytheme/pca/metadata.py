"""Associate every principal component with every metadata variable.

A PCA of pathway space says *how much* variance each component carries and
*which pathways* define it, but not *what* it corresponds to.  This module
answers the last part by testing each component's scores against each annotated
sample attribute, so a component can be attributed to a biological source
(diagnosis, tumour versus normal) or a technical one (sequencing batch, library
preparation, input quality) rather than eyeballed from a coloured scatter plot.

One test per (component, variable) cell:

* **categorical** variable -> Kruskal-Wallis across its levels, with
  eta-squared as the effect size, ``(H - k + 1) / (n - k)``, the fraction of the
  component's rank variance attributable to the variable.
* **continuous** variable -> Spearman correlation, whose square is reported in
  the same ``effect`` column so a heatmap of the grid is on one scale.

*p*-values are Benjamini-Hochberg corrected across the whole grid, because the
grid is the family actually being scanned.

The returned table is deliberately long rather than wide: one row per tested
cell, carrying the test used, the effect size, the raw and adjusted *p*, and for
a categorical variable the level with the highest mean score, which is what makes
the result readable ("PC1 separates CellLine from the rest").
"""
from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, NaN-safe."""
    out = np.full(p.shape, np.nan, dtype=float)
    ok = ~np.isnan(p)
    if not ok.any():
        return out
    vals = p[ok]
    n = vals.size
    order = np.argsort(vals)
    ranked = vals[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n, dtype=float)
    adj[order] = np.clip(ranked, 0.0, 1.0)
    out[ok] = adj
    return out


def _is_continuous(s: pd.Series, max_levels: int) -> bool:
    if not pd.api.types.is_numeric_dtype(s):
        return False
    return s.dropna().nunique() > max_levels


def associate_metadata(scores: pd.DataFrame,
                       metadata: pd.DataFrame,
                       *,
                       columns: Optional[Iterable[str]] = None,
                       min_group_size: int = 3,
                       max_levels_categorical: int = 12,
                       alpha: float = 0.05) -> pd.DataFrame:
    """Test every component in ``scores`` against every variable in ``metadata``.

    Parameters
    ----------
    scores
        Observations x components, as produced by the PCA stage
        (``PCAResult.scores``).
    metadata
        Observations x attributes.  Only rows shared with ``scores`` are used.
    columns
        Restrict to these metadata columns; default is every column that has at
        least two usable levels (or is continuous).
    min_group_size
        A categorical level with fewer observations than this is dropped before
        testing, and a variable left with fewer than two levels is skipped.
        Prevents a singleton level from producing a significant result.
    max_levels_categorical
        A numeric column with more distinct values than this is treated as
        continuous; at or below it, as categorical.  Lets a coded variable be
        recognised as a grouping.
    alpha
        Reported in the ``significant`` column; does not affect the statistics.

    Returns
    -------
    A long table, one row per tested (component, variable), sorted by adjusted
    *p*.  Empty if nothing was testable.
    """
    from scipy import stats

    shared = scores.index.intersection(metadata.index)
    if len(shared) < min_group_size * 2:
        return pd.DataFrame()
    S = scores.loc[shared]
    M = metadata.loc[shared]

    cols = list(columns) if columns is not None else list(M.columns)
    rows: list[dict] = []

    for var in cols:
        if var not in M.columns:
            continue
        raw = M[var]
        continuous = _is_continuous(raw, max_levels_categorical)

        if continuous:
            x = pd.to_numeric(raw, errors="coerce")
            usable = x.notna()
            if int(usable.sum()) < min_group_size * 2:
                continue
            for pc in S.columns:
                y = S.loc[usable, pc]
                rho, p = stats.spearmanr(x[usable], y)
                if not np.isfinite(rho):
                    continue
                rows.append({
                    "component": pc, "variable": var, "kind": "continuous",
                    "test": "spearman", "n": int(usable.sum()),
                    "n_levels": pd.NA, "statistic": float(rho),
                    "effect": float(rho ** 2), "effect_name": "rho_squared",
                    "p_value": float(p), "top_level": pd.NA,
                })
            continue

        # categorical: drop levels too small to test, then require >= 2 levels
        lab = raw.astype("object").where(raw.notna(), other=np.nan)
        counts = lab.value_counts(dropna=True)
        keep_levels = counts[counts >= min_group_size].index
        if len(keep_levels) < 2:
            continue
        mask = lab.isin(keep_levels)
        lab = lab[mask]
        k, n = len(keep_levels), int(mask.sum())
        for pc in S.columns:
            y = S.loc[mask, pc]
            groups = [y[lab == lv].to_numpy() for lv in keep_levels]
            if any(g.size < min_group_size for g in groups):
                continue
            try:
                H, p = stats.kruskal(*groups)
            except ValueError:            # identical values in every group
                continue
            eta2 = (H - k + 1) / (n - k) if n > k else np.nan
            means = {str(lv): float(y[lab == lv].mean()) for lv in keep_levels}
            rows.append({
                "component": pc, "variable": var, "kind": "categorical",
                "test": "kruskal", "n": n, "n_levels": k,
                "statistic": float(H),
                "effect": float(np.clip(eta2, 0.0, 1.0)),
                "effect_name": "eta_squared",
                "p_value": float(p),
                "top_level": max(means, key=means.get),
            })

    if not rows:
        return pd.DataFrame()

    out = pd.DataFrame(rows)
    out["fdr"] = _bh(out["p_value"].to_numpy(dtype=float))
    out["significant"] = out["fdr"] < alpha
    return out.sort_values(["fdr", "p_value"]).reset_index(drop=True)


def associate_levels(scores: pd.DataFrame,
                     metadata: pd.DataFrame,
                     *,
                     columns: Optional[Iterable[str]] = None,
                     min_group_size: int = 3,
                     max_levels_categorical: int = 12,
                     alpha: float = 0.05) -> pd.DataFrame:
    """One-versus-rest test of every *level* of every categorical variable.

    :func:`associate_metadata` runs one omnibus Kruskal-Wallis per
    (component, variable): it answers "does this component separate the
    diagnoses at all", not "which diagnosis".  A large eta-squared can be
    driven by a single level while the others are interchangeable, and the
    ``top_level`` column only names the level with the highest mean -- it is a
    direction hint, not a test.

    This function closes that gap.  For each component and each retained level
    it runs a Mann-Whitney U test of that level against every other retained
    sample of the same variable, and reports the rank-biserial correlation
    ``2U/(n1 n2) - 1`` as a signed effect size: ``+1`` means every sample of the
    level scores above every other sample, ``-1`` the reverse, ``0`` no
    separation.  ``auc``, the same quantity rescaled to ``[0, 1]``, is the
    probability that a randomly drawn sample of the level outscores a randomly
    drawn sample of the rest.

    A binary variable yields two mirror-image rows; only one of them enters the
    Benjamini-Hochberg family (the other is flagged ``redundant`` and inherits
    its adjusted *p*), so a two-level variable is not counted twice.

    Read the effect sizes, not only the flags.  The comparison is against the
    *pooled* remainder, matching the one-versus-all design used for the
    differential contrasts, so if one level is extreme the others can also come
    out significant purely because that level is inside their "rest".  The level
    with the largest ``|effect|`` on a component is the one the component is
    actually an axis for.

    Returns
    -------
    A long table, one row per tested (component, variable, level), sorted by
    adjusted *p*.  Empty if no categorical variable was testable.
    """
    from scipy import stats

    shared = scores.index.intersection(metadata.index)
    if len(shared) < min_group_size * 2:
        return pd.DataFrame()
    S = scores.loc[shared]
    M = metadata.loc[shared]

    cols = list(columns) if columns is not None else list(M.columns)
    rows: list[dict] = []

    for var in cols:
        if var not in M.columns:
            continue
        raw = M[var]
        if _is_continuous(raw, max_levels_categorical):
            continue                      # nothing to split into levels

        lab = raw.astype("object").where(raw.notna(), other=np.nan)
        counts = lab.value_counts(dropna=True)
        keep_levels = list(counts[counts >= min_group_size].index)
        if len(keep_levels) < 2:
            continue
        mask = lab.isin(keep_levels)
        lab = lab[mask]
        binary = len(keep_levels) == 2

        for pc in S.columns:
            y = S.loc[mask, pc]
            for j, lv in enumerate(keep_levels):
                inside = lab == lv
                a = y[inside].to_numpy()
                b = y[~inside].to_numpy()
                if a.size < min_group_size or b.size < min_group_size:
                    continue
                try:
                    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                except ValueError:        # all values tied
                    continue
                auc = float(U) / (a.size * b.size)
                rows.append({
                    "component": pc, "variable": var, "level": str(lv),
                    "test": "mannwhitney", "n_level": int(a.size),
                    "n_rest": int(b.size), "statistic": float(U),
                    "effect": 2.0 * auc - 1.0, "effect_name": "rank_biserial",
                    "auc": auc,
                    "mean_level": float(a.mean()), "mean_rest": float(b.mean()),
                    "p_value": float(p),
                    # a two-level variable's second row is the first one
                    # mirrored: identical p, opposite sign.  Kept for
                    # readability, excluded from the multiple-testing family.
                    "redundant": bool(binary and j == 1),
                })

    if not rows:
        return pd.DataFrame()

    out = pd.DataFrame(rows)
    fam = out.loc[~out["redundant"]].copy()
    fam["fdr"] = _bh(fam["p_value"].to_numpy(dtype=float))
    key = ["component", "variable"]
    # one fam row per (component, variable) for a binary variable, so first()
    # is that row; for a multi-level variable no row is redundant and this
    # lookup is never consulted.
    mirror = fam.groupby(key)["fdr"].first()
    out["fdr"] = np.where(
        out["redundant"].to_numpy(),
        pd.MultiIndex.from_frame(out[key]).map(mirror).to_numpy(dtype=float),
        out.index.map(fam["fdr"]).to_numpy(dtype=float))
    out["significant"] = out["fdr"] < alpha
    out["direction"] = np.where(out["effect"] >= 0, "high", "low")
    return out.sort_values(["fdr", "p_value"]).reset_index(drop=True)


def level_matrix(levels: pd.DataFrame, *, value: str = "effect",
                 variable: Optional[str] = None,
                 signed: bool = True) -> pd.DataFrame:
    """Pivot :func:`associate_levels` output to components x levels.

    Column labels are ``variable: level`` unless ``variable`` selects a single
    variable, in which case the levels are used bare.  ``signed=False`` takes
    the magnitude, for when only the strength of separation matters.
    """
    if levels.empty:
        return pd.DataFrame()
    df = levels if variable is None else levels[levels["variable"] == variable]
    if df.empty:
        return pd.DataFrame()
    df = df.copy()
    if not signed:
        df[value] = df[value].abs()
    df["_col"] = (df["level"] if variable is not None
                  else df["variable"] + ": " + df["level"])
    wide = df.pivot_table(index="component", columns="_col", values=value,
                          aggfunc="first")
    wide.columns.name = None
    order = sorted(wide.index, key=lambda s: int(str(s).lstrip("PC") or 0))
    return wide.loc[order]


def attribution_matrix(assoc: pd.DataFrame, *, value: str = "effect",
                       signed: bool = False) -> pd.DataFrame:
    """Pivot :func:`associate_metadata` output to components x variables.

    ``signed=True`` keeps the sign of a continuous association (Spearman rho)
    and leaves categorical cells positive, which is useful when the direction
    of a correlation matters; the default returns effect sizes only, all
    non-negative, so a single colour scale covers the grid.
    """
    if assoc.empty:
        return pd.DataFrame()
    df = assoc.copy()
    if signed:
        df[value] = np.where(df["kind"] == "continuous",
                             np.sign(df["statistic"]) * df[value], df[value])
    wide = df.pivot_table(index="component", columns="variable", values=value,
                          aggfunc="first")
    # keep PC order numeric rather than lexicographic (PC2 before PC10)
    order = sorted(wide.index, key=lambda s: int(str(s).lstrip("PC") or 0))
    return wide.loc[order]
