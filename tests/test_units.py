"""Unit tests for contracts, config, gene sets, io, grouping, and PCA."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathwaytheme.config import PipelineConfig, GroupingConfig, PCAConfig, EnrichmentConfig
from pathwaytheme.contracts import Grouping
from pathwaytheme.io import parse_gmt, gene_universe, load_input
from pathwaytheme.config import InputConfig
from pathwaytheme.grouping import resolve_grouping
from pathwaytheme.pca.core import run_pca, run_pca_scopes


# ── contracts ─────────────────────────────────────────────────────────────
def test_metadata_align_and_missing(score_matrix):
    md = score_matrix.metadata
    sub = md.align_to(["S02", "S00"])
    assert sub.samples == ["S02", "S00"]
    with pytest.raises(KeyError):
        md.align_to(["nope"])


def test_grouping_scopes():
    idx = ["a", "b", "c"]
    g = Grouping(labels=pd.Series(["x", "x", "y"], index=idx),
                 scope=pd.Series(["s1", "s1", "s2"], index=idx))
    assert g.scopes() == {"s1": ["a", "b"], "s2": ["c"]}
    g2 = Grouping(labels=pd.Series(["x", "y"], index=["a", "b"]))
    assert g2.scopes() == {"__all__": ["a", "b"]}


def test_scorematrix_subset(score_matrix):
    sub = score_matrix.subset_samples(["S00", "S01"])
    assert sub.samples == ["S00", "S01"]
    assert sub.data.shape == (40, 2)


# ── config ────────────────────────────────────────────────────────────────
def test_config_yaml_roundtrip(tmp_path):
    cfg = PipelineConfig()
    cfg.enrichment.backend = "enrichr"
    cfg.pca.k_max = 7
    p = tmp_path / "c.yaml"
    cfg.to_yaml(p)
    back = PipelineConfig.from_yaml(p)
    assert back.enrichment.backend == "enrichr"
    assert back.pca.k_max == 7


# ── gene sets ───────────────────────────────────────────────────────────────
def test_parse_gmt(synthetic_gmt):
    path, _ = synthetic_gmt
    gs = parse_gmt(path)
    assert len(gs) == 15
    assert all(len(v) >= 8 for v in gs.values())
    assert len(gene_universe(gs)) <= 200


# ── io matrix adapter ───────────────────────────────────────────────────────
def test_matrix_adapter_roundtrip(feature_matrix, tmp_path):
    mp, mdp = tmp_path / "m.tsv", tmp_path / "md.tsv"
    feature_matrix.data.to_csv(mp, sep="\t")
    feature_matrix.metadata.table.to_csv(mdp, sep="\t")
    fm = load_input(InputConfig(kind="matrix", matrix_path=str(mp), metadata_path=str(mdp)))
    assert fm.shape == feature_matrix.shape
    assert np.allclose(fm.data.values, feature_matrix.data.values)


def test_matrix_adapter_missing_metadata_raises(feature_matrix, tmp_path):
    mp, mdp = tmp_path / "m.tsv", tmp_path / "md.tsv"
    feature_matrix.data.to_csv(mp, sep="\t")
    feature_matrix.metadata.table.iloc[:3].to_csv(mdp, sep="\t")  # incomplete
    with pytest.raises(KeyError):
        load_input(InputConfig(kind="matrix", matrix_path=str(mp), metadata_path=str(mdp)))


# ── grouping ────────────────────────────────────────────────────────────────
def test_grouping_target(score_matrix):
    g = resolve_grouping(score_matrix, GroupingConfig(mode="target", target_col="grp",
                                                      scope_col="sample_id"))
    assert set(g.labels.unique()) == {"hi", "lo"}
    assert set(g.scopes()) == {"A", "B"}


def test_grouping_target_missing_col(score_matrix):
    with pytest.raises(KeyError):
        resolve_grouping(score_matrix, GroupingConfig(mode="target", target_col="ghost"))


def test_grouping_auto(score_matrix):
    g = resolve_grouping(score_matrix, GroupingConfig(mode="auto", k_min=2, k_max=4))
    assert g.labels.nunique() >= 2


# ── pca ─────────────────────────────────────────────────────────────────────
def test_run_pca_shapes_and_variance(score_matrix):
    labels = score_matrix.metadata.get("grp")
    res = run_pca(score_matrix.data, labels, PCAConfig())
    assert res is not None
    K = res.K
    assert res.scores.shape == (12, K)
    assert res.loadings.shape == (K, 40)
    assert 0 < res.variance_explained.sum() <= 1.0 + 1e-9
    # signatures: one block per observation, top-25 (or fewer) each
    assert res.signatures["cluster"].nunique() == 12


def test_run_pca_too_few_obs_returns_none(score_matrix):
    sub = score_matrix.subset_samples(["S00", "S01"])
    res = run_pca(sub.data, sub.metadata.get("grp"), PCAConfig(min_observations=3))
    assert res is None


def test_run_pca_scopes(score_matrix):
    g = resolve_grouping(score_matrix, GroupingConfig(mode="target", target_col="grp",
                                                      scope_col="sample_id"))
    results = run_pca_scopes(score_matrix, g, PCAConfig(), size_col="n_cells")
    assert {r.scope for r in results} == {"A", "B"}
    for r in results:
        assert r.sizes is not None and len(r.sizes) == 6


def test_run_pca_deterministic(score_matrix):
    labels = score_matrix.metadata.get("grp")
    r1 = run_pca(score_matrix.data, labels, PCAConfig())
    r2 = run_pca(score_matrix.data, labels, PCAConfig())
    assert np.allclose(r1.variance_explained, r2.variance_explained)
    assert np.allclose(r1.scores.values, r2.scores.values)
