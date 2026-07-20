"""Tests for the differential stage, category roll-up, sample filter, QC heatmap."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import special, stats

import pathwaytheme as pt
from pathwaytheme.config import DiffConfig
from pathwaytheme.contracts import Grouping
from pathwaytheme.diff import run_diff
from pathwaytheme.diff.stats import (benjamini_hochberg, fit_fdist, moderated_t,
                                     mannwhitney, trigamma_inverse, welch)
from pathwaytheme.categories import summarize_by_category


# ── low-level stats ─────────────────────────────────────────────────────────
def test_trigamma_inverse_roundtrip():
    for x in [0.01, 0.1, 0.5, 1.0, 5.0, 50.0]:
        y = trigamma_inverse(x)
        assert np.isclose(special.polygamma(1, y), x, rtol=1e-6)


def test_benjamini_hochberg_matches_reference():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.9])
    # BH by hand: sort, p*m/rank, cumulative-min from the top
    m = len(p)
    expect = np.minimum.accumulate((np.sort(p) * m / np.arange(1, m + 1))[::-1])[::-1]
    got = benjamini_hochberg(p)
    assert np.allclose(np.sort(got), np.sort(expect))
    assert np.all(got >= p - 1e-12)  # adjusted >= raw


def test_benjamini_hochberg_handles_nan():
    p = np.array([0.01, np.nan, 0.5])
    got = benjamini_hochberg(p)
    assert np.isnan(got[1])
    assert np.isfinite(got[0]) and np.isfinite(got[2])


def test_welch_matches_scipy(rng):
    case = rng.normal(1.0, 1.0, size=(5, 8))
    ref = rng.normal(0.0, 1.5, size=(5, 10))
    effect, stat, pval = welch(case, ref)
    ref_res = stats.ttest_ind(case, ref, axis=1, equal_var=False)
    assert np.allclose(stat, ref_res.statistic)
    assert np.allclose(pval, ref_res.pvalue)
    assert np.allclose(effect, case.mean(1) - ref.mean(1))


def test_mannwhitney_matches_scipy(rng):
    case = rng.normal(2.0, 1.0, size=(3, 9))
    ref = rng.normal(0.0, 1.0, size=(3, 11))
    _, stat, pval = mannwhitney(case, ref)
    for i in range(3):
        u, p = stats.mannwhitneyu(case[i], ref[i], alternative="two-sided")
        assert np.isclose(stat[i], u)
        assert np.isclose(pval[i], p)


def test_fit_fdist_recovers_prior():
    """Draw residual variances from a known scaled-F prior; fit_fdist recovers it."""
    rng = np.random.default_rng(7)
    d0_true, s0_true, dfr = 8.0, 2.0, 6.0
    n = 4000
    # s2 ~ s0^2 * (chi2_{dfr}/dfr) marginalised over prior scale chi2_{d0}
    prior_scale = s0_true * d0_true / rng.chisquare(d0_true, size=n)
    s2 = prior_scale * rng.chisquare(dfr, size=n) / dfr
    d0_hat, s0_hat = fit_fdist(s2, np.full(n, dfr))
    assert np.isclose(d0_hat, d0_true, rtol=0.30)
    assert np.isclose(s0_hat, s0_true, rtol=0.20)


def test_moderated_t_sign_and_power():
    rng = np.random.default_rng(1)
    # 30 null pathways + 10 with a real +3 shift in the case group
    case = rng.normal(0, 1, size=(40, 8))
    ref = rng.normal(0, 1, size=(40, 8))
    case[:10] += 3.0
    effect, stat, pval = moderated_t(case, ref)
    assert np.all(np.sign(stat[:10]) == np.sign(effect[:10]))
    # the shifted pathways should be far more significant than the nulls
    assert np.median(pval[:10]) < np.median(pval[10:])
    assert (pval[:10] < 0.05).mean() >= 0.8


def test_moderated_t_shrinks_toward_pooled_t():
    """With many pathways sharing a variance scale, moderated t is finite and
    correlates strongly with the ordinary pooled t-statistic."""
    rng = np.random.default_rng(3)
    case = rng.normal(0, 1, size=(200, 6))
    ref = rng.normal(0, 1, size=(200, 6))
    case[:50] += 2.0
    _, mod_stat, _ = moderated_t(case, ref)
    # ordinary pooled Student t
    pooled = stats.ttest_ind(case, ref, axis=1, equal_var=True).statistic
    ok = np.isfinite(mod_stat) & np.isfinite(pooled)
    r = np.corrcoef(mod_stat[ok], pooled[ok])[0, 1]
    assert r > 0.95


# ── run_diff orchestration ──────────────────────────────────────────────────
@pytest.mark.parametrize("method", ["welch", "mannwhitney", "moderated_t"])
def test_run_diff_one_vs_rest(score_matrix, method):
    res = run_diff(score_matrix, DiffConfig(method=method, group_col="grp"))
    # two groups -> one-vs-rest gives 2 comparisons
    assert set(res.comparisons) == {"hi_vs_rest", "lo_vs_rest"}
    cols = {"comparison", "pathway", "effect", "p_value", "fdr", "direction"}
    assert cols <= set(res.table.columns)
    # the 10 elevated terms should surface as significant in hi_vs_rest
    hi = res.table[res.table["comparison"] == "hi_vs_rest"]
    top10 = set(hi.nsmallest(10, "fdr")["pathway"])
    elevated = {f"TERM_{i:02d}" for i in range(10)}
    assert len(top10 & elevated) >= 8


def test_run_diff_reference_and_contrasts(score_matrix):
    ref = run_diff(score_matrix, DiffConfig(method="welch", group_col="grp",
                                            reference="lo"))
    assert ref.comparisons == ["hi_vs_lo"]
    con = run_diff(score_matrix, DiffConfig(method="welch", group_col="grp",
                                            contrasts=[["hi", "lo"]]))
    assert con.comparisons == ["hi_vs_lo"]
    assert np.allclose(
        ref.table.sort_values("pathway")["statistic"].to_numpy(),
        con.table.sort_values("pathway")["statistic"].to_numpy(), equal_nan=True)


def test_run_diff_uses_grouping_when_no_group_col(score_matrix):
    labels = score_matrix.metadata.get("grp")
    g = Grouping(labels=labels)
    res = run_diff(score_matrix, DiffConfig(method="welch"), g)
    assert set(res.comparisons) == {"hi_vs_rest", "lo_vs_rest"}


# ── categories ──────────────────────────────────────────────────────────────
def test_summarize_by_category():
    table = pd.DataFrame({
        "comparison": ["c"] * 4,
        "pathway": ["p1", "p2", "p3", "p4"],
        "effect": [2.0, -1.0, 4.0, 0.5],
    })
    mapping = {"p1": "A", "p2": "A", "p3": "B", "p4": "B"}
    out = summarize_by_category(table, mapping, value_col="effect",
                                group_cols=["comparison"], stat="mean")
    a = out[out["category"] == "A"].iloc[0]
    assert np.isclose(a["mean_effect"], 0.5)
    assert a["n_pathways"] == 2
    assert a["n_up"] == 1 and a["n_down"] == 1


def test_summarize_unmapped_go_to_other():
    table = pd.DataFrame({"pathway": ["p1", "pX"], "effect": [1.0, 2.0]})
    out = summarize_by_category(table, {"p1": "A"}, value_col="effect",
                                unmapped_label="Other")
    assert "Other" in set(out["category"])


def test_summarize_significance_split():
    table = pd.DataFrame({
        "pathway": ["p1", "p2", "p3", "p4"],
        "effect": [2.0, 1.0, -1.0, 0.5],
        "fdr": [0.001, 0.20, 0.01, 0.90],   # p1,p3 sig; p2,p4 not
    })
    mapping = {"p1": "A", "p2": "A", "p3": "B", "p4": "B"}
    out = summarize_by_category(table, mapping, value_col="effect",
                                significance_col="fdr", alpha=0.05)
    assert "significance" in out.columns
    assert set(out["significance"]) == {"significant", "non_significant"}
    a_sig = out[(out["category"] == "A") & (out["significance"] == "significant")]
    a_ns = out[(out["category"] == "A") & (out["significance"] == "non_significant")]
    assert a_sig.iloc[0]["n_pathways"] == 1        # only p1
    assert a_ns.iloc[0]["n_pathways"] == 1         # only p2


def test_summarize_significance_nan_is_nonsignificant():
    table = pd.DataFrame({"pathway": ["p1"], "effect": [1.0], "fdr": [np.nan]})
    out = summarize_by_category(table, {"p1": "A"}, value_col="effect",
                                significance_col="fdr", alpha=0.05)
    assert out.iloc[0]["significance"] == "non_significant"


def test_summarize_categories_api_significance(score_matrix):
    d = pt.diff(score_matrix, group_col="grp", method="welch", reference="lo")
    mapping = {t: ("early" if i < 20 else "late")
               for i, t in enumerate(f"TERM_{i:02d}" for i in range(40))}
    summ = pt.summarize_categories(d, mapping, value_col="effect",
                                   significance_col="fdr", alpha=0.05)
    assert "significance" in summ.columns


# ── api helpers ─────────────────────────────────────────────────────────────
def test_filter_samples_matrix(feature_matrix):
    sub = pt.filter_samples(feature_matrix, "grp", ["hi"])
    assert sub.shape[1] == 6
    assert set(sub.metadata.get("grp")) == {"hi"}


def test_filter_samples_bad_value(feature_matrix):
    with pytest.raises(ValueError):
        pt.filter_samples(feature_matrix, "grp", ["nope"])


def test_sanity_heatmap_writes_pdf(score_matrix, tmp_path):
    out = tmp_path / "sanity.pdf"
    p = pt.sanity_heatmap(score_matrix, str(out), label_col="grp", max_pathways=0)
    assert p is not None and out.exists() and out.stat().st_size > 0


def test_diff_api_and_category_roundtrip(score_matrix):
    d = pt.diff(score_matrix, group_col="grp", method="welch", reference="lo")
    assert isinstance(d.table, pd.DataFrame) and not d.table.empty
    mapping = {t: ("early" if i < 20 else "late")
               for i, t in enumerate(f"TERM_{i:02d}" for i in range(40))}
    summ = pt.summarize_categories(d, mapping, value_col="effect")
    assert {"category", "mean_effect", "n_pathways"} <= set(summ.columns)
