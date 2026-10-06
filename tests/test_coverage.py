"""Gene-set coverage accounting + the missing-value policy.

These cover the two things that must not happen silently when the same gene
sets are run over a modality with sparser feature coverage than RNA -- a set
disappearing from the output, and a missing value being imputed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathwaytheme.config import EnrichmentConfig
from pathwaytheme.contracts import FeatureMatrix, SampleMetadata
from pathwaytheme.enrichment import run_enrichment
from pathwaytheme.enrichment.coverage import (
    NO_OVERLAP, SCORED, TOO_LARGE, TOO_SMALL, apply_missing_policy,
    geneset_coverage)


def test_coverage_classifies_every_fate():
    features = ["A", "B", "C", "D", "E", "F"]
    gs = {
        "FULL": ["A", "B", "C", "D", "E"],       # 5 present -> scored
        "PARTIAL": ["A", "B", "Z1", "Z2"],       # 2 present -> below min 5
        "ABSENT": ["Z3", "Z4", "Z5"],            # 0 present -> no overlap
        "HUGE": features + ["Z6"],               # 6 present -> above max 5
    }
    cov = geneset_coverage(features, gs, min_size=5, max_size=5)
    got = dict(zip(cov["pathway"], cov["status"]))
    assert got == {"FULL": SCORED, "PARTIAL": TOO_SMALL,
                   "ABSENT": NO_OVERLAP, "HUGE": TOO_LARGE}
    row = cov.set_index("pathway").loc["PARTIAL"]
    assert (row["n_genes"], row["n_present"]) == (4, 2)
    assert row["fraction_present"] == pytest.approx(0.5)


def test_coverage_matches_case_insensitively():
    cov = geneset_coverage(["abc", "def"], {"S": ["ABC", "DEF"]},
                           min_size=2, max_size=100)
    assert cov.loc[0, "n_present"] == 2


def test_dropped_sets_are_reported_not_silent(feature_matrix, synthetic_gmt,
                                              tmp_path):
    """A set whose members are absent must show up in coverage, with a warning."""
    gmt, genes = synthetic_gmt
    extra = tmp_path / "with_absent.gmt"
    text = gmt.read_text(encoding="utf-8")
    text += "\nABSENT_SET\tdesc\t" + "\t".join(f"NOSUCHGENE{i}" for i in range(9))
    extra.write_text(text, encoding="utf-8")

    cfg = EnrichmentConfig(backend="ssgsea", geneset="SYN", gmt_path=str(extra),
                           min_gene_set_size=5, threads=2)
    with pytest.warns(UserWarning, match="ABSENT_SET"):
        sm = run_enrichment(feature_matrix, cfg)

    assert sm.coverage is not None
    fate = sm.coverage.set_index("pathway")["status"]
    assert fate["ABSENT_SET"] == NO_OVERLAP
    # the set is genuinely absent from the scores, and the table says why
    assert "ABSENT_SET" not in sm.data.index
    n_scored = int((sm.coverage["status"] == SCORED).sum())
    assert sm.shape[0] == n_scored


def test_coverage_survives_the_cache(feature_matrix, synthetic_gmt, tmp_path):
    gmt, _ = synthetic_gmt
    cfg = EnrichmentConfig(backend="ssgsea", geneset="SYN", gmt_path=str(gmt),
                           min_gene_set_size=5, threads=2,
                           cache_dir=str(tmp_path / "cache"))
    first = run_enrichment(feature_matrix, cfg)
    second = run_enrichment(feature_matrix, cfg)      # served from cache
    assert second.coverage is not None
    pd.testing.assert_frame_equal(first.coverage, second.coverage)


# ── missing-value policy ──────────────────────────────────────────────────
def _mat_with_gaps():
    d = pd.DataFrame(np.arange(1.0, 25.0).reshape(6, 4),
                     index=[f"F{i}" for i in range(6)],
                     columns=list("WXYZ"))
    d.iloc[0, 0] = np.nan            # F0 observed in 3 of 4
    d.iloc[1, :] = np.nan            # F1 observed in none
    d.iloc[2, :2] = np.nan           # F2 observed in 2 of 4
    return d


def test_policy_as_is_leaves_values_but_warns():
    d = _mat_with_gaps()
    with pytest.warns(UserWarning, match="as_is"):
        out = apply_missing_policy(d, "as_is", 0.0, verbose=True)
    assert "F1" not in out.index          # all-missing row always goes
    assert out.isna().sum().sum() == 3    # the rest are untouched
    assert out.loc["F0", "X"] == d.loc["F0", "X"]


def test_policy_zero_fills_explicitly():
    out = apply_missing_policy(_mat_with_gaps(), "zero", 0.0, verbose=False)
    assert out.isna().sum().sum() == 0
    assert out.loc["F0", "W"] == 0.0


def test_policy_low_ranks_missing_below_every_measured_value():
    out = apply_missing_policy(_mat_with_gaps(), "low", 0.0, verbose=False)
    assert out.isna().sum().sum() == 0
    for col in out.columns:
        filled = out.loc[["F0", "F2"], col]
        observed = out[col].drop(index=["F0", "F2"])
        assert (filled < observed.min()).all() or filled.isna().all()
    # and the fill is per-sample, not one global constant
    assert out.loc["F2", "W"] != out.loc["F2", "X"]


def test_min_observed_fraction_drops_under_observed_features():
    out = apply_missing_policy(_mat_with_gaps(), "zero", 0.75, verbose=False)
    assert list(out.index) == ["F0", "F3", "F4", "F5"]   # F1 none, F2 only 50%


def test_unknown_policy_is_rejected():
    with pytest.raises(ValueError, match="must be one of"):
        apply_missing_policy(_mat_with_gaps(), "impute_knn", 0.0, verbose=False)


def test_policy_default_is_a_noop_on_a_complete_matrix(feature_matrix,
                                                       synthetic_gmt):
    """The new options must not perturb a matrix that has no missing values."""
    gmt, _ = synthetic_gmt
    assert feature_matrix.data.isna().sum().sum() == 0
    cfg = EnrichmentConfig(backend="ssgsea", geneset="SYN", gmt_path=str(gmt),
                           min_gene_set_size=5, threads=2)
    a = run_enrichment(feature_matrix, cfg)
    b = run_enrichment(feature_matrix, EnrichmentConfig(
        backend="ssgsea", geneset="SYN", gmt_path=str(gmt),
        min_gene_set_size=5, threads=2, missing="low",
        min_observed_fraction=0.5))
    pd.testing.assert_frame_equal(a.data, b.data)


def test_fill_actually_reaches_the_backend(feature_matrix, synthetic_gmt):
    """Regression: filling changes no shape, so a shape-based 'did anything
    change?' test discards the filled matrix and lets gseapy's own zero-fill
    win.  The policy must alter the scores."""
    from pathwaytheme.enrichment.base import _prepare_features

    gmt, _ = synthetic_gmt
    gappy = feature_matrix.data.copy()
    gappy.iloc[:60, :4] = np.nan          # holes, but no empty row or column
    fm = FeatureMatrix(gappy, feature_matrix.metadata)

    prepared = _prepare_features(
        fm, EnrichmentConfig(missing="low"), verbose=False)
    assert prepared.data.shape == fm.data.shape      # shape is unchanged ...
    assert prepared.data.isna().sum().sum() == 0     # ... but values are not
    assert prepared is not fm

    base = EnrichmentConfig(backend="ssgsea", geneset="SYN", gmt_path=str(gmt),
                            min_gene_set_size=5, threads=2)
    with pytest.warns(UserWarning):
        as_is = run_enrichment(fm, base)
    low = run_enrichment(fm, EnrichmentConfig(
        backend="ssgsea", geneset="SYN", gmt_path=str(gmt),
        min_gene_set_size=5, threads=2, missing="low"))
    common = as_is.data.index.intersection(low.data.index)
    assert len(common) > 0
    assert not np.allclose(as_is.data.loc[common].values,
                           low.data.loc[common].values)


def test_missing_policy_changes_the_cache_key(tmp_path):
    from pathwaytheme.enrichment.base import _cache_path
    base = dict(backend="ssgsea", geneset="SYN", cache_dir=str(tmp_path))
    p_default = _cache_path(EnrichmentConfig(**base))
    p_low = _cache_path(EnrichmentConfig(**base, missing="low"))
    p_obs = _cache_path(EnrichmentConfig(**base, min_observed_fraction=0.75))
    assert p_default != p_low != p_obs and p_default != p_obs
    assert p_default.name == "SYN_ssgsea_scores.tsv"   # unchanged for defaults
