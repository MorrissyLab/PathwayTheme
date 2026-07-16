"""Enrichment-backend tests + a synthetic end-to-end pipeline run.

All offline / deterministic: ssGSEA and the over-representation backends run
against the synthetic .gmt; no network access is required.
"""

from __future__ import annotations

import numpy as np
import pytest

from pathwaytheme.config import (PipelineConfig, InputConfig, EnrichmentConfig,
                                 GroupingConfig, PCAConfig, VizConfig)
from pathwaytheme.enrichment import run_enrichment, get_backend
from pathwaytheme.pipeline import run


def test_backend_registry():
    for name in ("ssgsea", "enrichr", "goslim"):
        assert get_backend(name).name == name
    with pytest.raises(ValueError):
        get_backend("does_not_exist")


def test_ssgsea_backend(feature_matrix, synthetic_gmt):
    gmt, _ = synthetic_gmt
    sm = run_enrichment(feature_matrix, EnrichmentConfig(
        backend="ssgsea", geneset="SYN", gmt_path=str(gmt),
        min_gene_set_size=5, threads=2))
    assert sm.shape[1] == feature_matrix.shape[1]
    assert sm.shape[0] >= 1
    assert list(sm.samples) == list(feature_matrix.samples)


def test_enrichr_offline_backend(feature_matrix, synthetic_gmt):
    gmt, _ = synthetic_gmt
    sm = run_enrichment(feature_matrix, EnrichmentConfig(
        backend="enrichr", geneset="SYN_ORA", gmt_path=str(gmt),
        gene_selection="top_n", top_n=60))
    assert sm.shape[1] >= 1
    assert (sm.data.values >= 0).all()  # -log10(padj) is non-negative


def test_goslim_gmt_fraction(feature_matrix, synthetic_gmt):
    gmt, _ = synthetic_gmt
    sm = run_enrichment(feature_matrix, EnrichmentConfig(
        backend="goslim", geneset="SLIM", gmt_path=str(gmt),
        gene_selection="top_n", top_n=60, goslim_score="fraction"))
    assert sm.shape[1] == feature_matrix.shape[1]
    assert (sm.data.values >= 0).all() and (sm.data.values <= 1).all()


def test_enrichment_cache(feature_matrix, synthetic_gmt, tmp_path):
    gmt, _ = synthetic_gmt
    cfg = EnrichmentConfig(backend="goslim", geneset="SLIM", gmt_path=str(gmt),
                           top_n=60, cache_dir=str(tmp_path / "cache"))
    sm1 = run_enrichment(feature_matrix, cfg)
    cache_file = tmp_path / "cache" / "SLIM_goslim_scores.tsv"
    assert cache_file.exists()
    sm2 = run_enrichment(feature_matrix, cfg)  # served from cache
    assert np.allclose(sm1.data.values, sm2.data.loc[sm1.terms, sm1.samples].values)


def test_pipeline_end_to_end(feature_matrix, synthetic_gmt, tmp_path, monkeypatch):
    """Matrix input -> goslim -> target grouping -> per-scope PCA -> figs+tables."""
    gmt, _ = synthetic_gmt
    # write the feature matrix out so the matrix adapter loads it
    mp, mdp = tmp_path / "m.tsv", tmp_path / "md.tsv"
    feature_matrix.data.to_csv(mp, sep="\t")
    feature_matrix.metadata.table.to_csv(mdp, sep="\t")

    out = tmp_path / "out"
    cfg = PipelineConfig(
        input=InputConfig(kind="matrix", matrix_path=str(mp), metadata_path=str(mdp)),
        enrichment=EnrichmentConfig(backend="goslim", geneset="SLIM",
                                    gmt_path=str(gmt), top_n=80),
        grouping=GroupingConfig(mode="target", target_col="cluster_id",
                                scope_col="sample_id", size_col="n_cells"),
        pca=PCAConfig(min_features=5),
        viz=VizConfig(make_figures=True, make_tables=True),
        output_dir=str(out),
    )
    result = run(cfg, verbose=False)
    assert {r.scope for r in result.pca_results} == {"A", "B"}
    root = out / "SLIM" / "sample_pca"
    for scope in ("A", "B"):
        d = root / scope
        assert sum(p.suffix == ".tsv" for p in d.iterdir()) == 3
        assert sum(p.suffix == ".pdf" for p in d.iterdir()) == 11


def test_pipeline_tables_only(feature_matrix, synthetic_gmt, tmp_path):
    gmt, _ = synthetic_gmt
    mp, mdp = tmp_path / "m.tsv", tmp_path / "md.tsv"
    feature_matrix.data.to_csv(mp, sep="\t")
    feature_matrix.metadata.table.to_csv(mdp, sep="\t")
    cfg = PipelineConfig(
        input=InputConfig(kind="matrix", matrix_path=str(mp), metadata_path=str(mdp)),
        enrichment=EnrichmentConfig(backend="goslim", geneset="SLIM",
                                    gmt_path=str(gmt), top_n=80),
        grouping=GroupingConfig(mode="target", target_col="cluster_id",
                                scope_col="sample_id"),
        pca=PCAConfig(min_features=5),
        viz=VizConfig(make_figures=False, make_tables=True),
        output_dir=str(tmp_path / "out2"),
    )
    result = run(cfg, verbose=False)
    for scope in ("A", "B"):
        d = tmp_path / "out2" / "SLIM" / "sample_pca" / scope
        assert sum(p.suffix == ".pdf" for p in d.iterdir()) == 0
        assert sum(p.suffix == ".tsv" for p in d.iterdir()) == 3
