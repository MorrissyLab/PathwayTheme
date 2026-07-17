"""Per-pathway two-group statistics used by the differential stage.

Three methods, all operating on one pathway x sample score matrix split into a
*case* and a *reference* set of columns:

* ``welch``        — Welch's unequal-variance t-test (scipy).
* ``mannwhitney``  — Mann-Whitney U rank test (scipy); matches the
  ``wilcox_pathway_stats`` helper in the reference ASPS R workflow.
* ``moderated_t``  — limma's empirical-Bayes moderated t-statistic (Smyth 2004).
  Python is missing limma, so ``fit_fdist`` / ``trigamma_inverse`` below are
  direct ports of ``limma::fitFDist`` / ``limma::trigammaInverse``; the two-group
  case reduces to a pooled linear model with variance shrinkage.

Each method returns three equal-length arrays over the pathways (rows):
``effect`` (mean_case - mean_reference), ``statistic``, ``p_value``.
"""

from __future__ import annotations

import numpy as np
from scipy import special, stats


# ── limma variance-model helpers (ports of limma R functions) ──────────────
def trigamma_inverse(x: float) -> float:
    """Solve ``trigamma(y) = x`` for ``y`` (port of ``limma::trigammaInverse``)."""
    x = float(x)
    if not np.isfinite(x):
        return np.nan
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)              # trigamma(y)
        dif = tri * (1.0 - tri / x) / special.polygamma(2, y)  # /tetragamma
        y = y + dif
        if -dif / y < 1e-8:
            break
    return y


def fit_fdist(var: np.ndarray, df1: np.ndarray) -> tuple[float, float]:
    """Estimate the scaled-F prior for a set of sample variances.

    Port of ``limma::fitFDist``.  Given per-pathway residual variances ``var``
    each on ``df1`` degrees of freedom, returns ``(d0, s0_sq)`` — the prior
    degrees of freedom and prior variance for empirical-Bayes shrinkage.
    """
    var = np.asarray(var, dtype=float)
    df1 = np.broadcast_to(np.asarray(df1, dtype=float), var.shape)

    ok = np.isfinite(var) & (var > 0) & np.isfinite(df1) & (df1 > 0)
    x = var[ok]
    g1 = df1[ok]
    if x.size == 0:
        return np.inf, np.nan

    z = np.log(x)
    e = z - special.digamma(g1 / 2.0) + np.log(g1 / 2.0)
    emean = float(np.mean(e))
    n = e.size
    if n <= 1:
        return np.inf, float(np.exp(emean))

    evar = float(np.sum((e - emean) ** 2) / (n - 1))
    evar = evar - float(np.mean(special.polygamma(1, g1 / 2.0)))  # subtract mean trigamma

    if evar > 0:
        d0 = 2.0 * trigamma_inverse(evar)
        s0_sq = float(np.exp(emean + special.digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    else:
        d0 = np.inf
        s0_sq = float(np.exp(emean))
    return d0, s0_sq


# ── the three tests ────────────────────────────────────────────────────────
def welch(case: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Welch's t-test per row.  ``case``/``ref`` are ``pathways x samples``."""
    effect = np.nanmean(case, axis=1) - np.nanmean(ref, axis=1)
    res = stats.ttest_ind(case, ref, axis=1, equal_var=False, nan_policy="omit")
    stat = np.asarray(res.statistic, dtype=float)
    pval = np.asarray(res.pvalue, dtype=float)
    return effect, stat, pval


def mannwhitney(case: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mann-Whitney U per row (two-sided)."""
    n = case.shape[0]
    effect = np.nanmean(case, axis=1) - np.nanmean(ref, axis=1)
    stat = np.full(n, np.nan)
    pval = np.full(n, np.nan)
    for i in range(n):
        x = case[i][~np.isnan(case[i])]
        y = ref[i][~np.isnan(ref[i])]
        if x.size < 1 or y.size < 1:
            continue
        try:
            u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            stat[i], pval[i] = float(u), float(p)
        except ValueError:
            # all values identical -> undefined; leave as NaN
            continue
    return effect, stat, pval


def moderated_t(case: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """limma empirical-Bayes moderated t-test for a two-group comparison.

    Fits the pooled two-group linear model per pathway, shrinks the residual
    variances toward a scaled-F prior (``fit_fdist``), and computes the
    moderated t-statistic on ``d0 + d`` degrees of freedom.
    """
    n1 = np.sum(~np.isnan(case), axis=1).astype(float)
    n2 = np.sum(~np.isnan(ref), axis=1).astype(float)
    m1 = np.nanmean(case, axis=1)
    m2 = np.nanmean(ref, axis=1)
    effect = m1 - m2

    ss1 = np.nansum((case - m1[:, None]) ** 2, axis=1)
    ss2 = np.nansum((ref - m2[:, None]) ** 2, axis=1)
    dfr = n1 + n2 - 2.0                              # residual df per pathway
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = (ss1 + ss2) / dfr                       # ordinary residual variance
        stdev_unscaled = np.sqrt(1.0 / n1 + 1.0 / n2)

    # empirical-Bayes shrinkage of the residual variances
    good = np.isfinite(s2) & (dfr > 0)
    d0, s0_sq = fit_fdist(s2[good], dfr[good])

    if np.isinf(d0):
        s2_post = np.full_like(s2, s0_sq)
        df_total = np.full_like(dfr, np.inf)
    else:
        s2_post = (d0 * s0_sq + dfr * s2) / (d0 + dfr)
        df_total = dfr + d0

    with np.errstate(divide="ignore", invalid="ignore"):
        stat = effect / (stdev_unscaled * np.sqrt(s2_post))

    pval = np.full_like(stat, np.nan)
    finite = np.isfinite(stat)
    inf_df = finite & ~np.isfinite(df_total)
    fin_df = finite & np.isfinite(df_total) & (df_total > 0)
    pval[fin_df] = 2.0 * stats.t.sf(np.abs(stat[fin_df]), df_total[fin_df])
    pval[inf_df] = 2.0 * stats.norm.sf(np.abs(stat[inf_df]))
    return effect, stat, pval


_METHODS = {
    "welch": welch,
    "mannwhitney": mannwhitney,
    "moderated_t": moderated_t,
}


def run_test(method: str, case: np.ndarray, ref: np.ndarray):
    try:
        fn = _METHODS[method]
    except KeyError:
        raise ValueError(f"Unknown diff method {method!r}; expected one of "
                         f"{sorted(_METHODS)}")
    return fn(case, ref)


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """BH-FDR adjusted p-values; NaNs pass through and are excluded from ranking."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * m / (np.arange(1, m + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out
