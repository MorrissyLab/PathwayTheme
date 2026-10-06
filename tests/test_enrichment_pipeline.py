"""Enrichment-backend tests + a synthetic end-to-end pipeline run.

All offline / deterministic: ssGSEA and the over-representation backends run
against the synthetic .gmt; no network access is required.
"""

from __future__ import annotations

import numpy as np
import pytest

from pathwaytheme.config import (PipelineConfig, InputConfig, EnrichmentConfig,
                                 GroupingConfig, PCAConfig, VizConfig)
from pathwaytheme.contracts import FeatureMatrix
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


def test_cache_partial_hit_scores_missing_samples(feature_matrix, synthetic_gmt,
                                                  tmp_path):
    """A sample absent from the cache must be scored, never silently dropped.

    The cache key is backend+geneset, so two configs selecting different
    samples share one file.  Priming it with a subset and then asking for
    everything used to return only the subset.
    """
    gmt, _ = synthetic_gmt
    cfg = EnrichmentConfig(backend="goslim", geneset="SLIM", gmt_path=str(gmt),
                           top_n=60, cache_dir=str(tmp_path / "cache"))
    all_samples = feature_matrix.samples
    head = all_samples[:2]
    subset = FeatureMatrix(feature_matrix.data.loc[:, head],
                           feature_matrix.metadata.align_to(head))

    primed = run_enrichment(subset, cfg)          # cache holds 2 samples
    assert primed.shape[1] == 2

    full = run_enrichment(feature_matrix, cfg)    # ask for all of them
    assert full.samples == all_samples
    assert full.shape[1] == len(all_samples)
    # the samples already cached keep exactly their cached values
    assert np.allclose(primed.data.values,
                       full.data.loc[primed.terms, head].values)
    # and a fresh run with no cache agrees with the filled-in one
    clean = run_enrichment(feature_matrix, EnrichmentConfig(
        backend="goslim", geneset="SLIM", gmt_path=str(gmt), top_n=60))
    assert np.allclose(clean.data.loc[full.terms, full.samples].values,
                       full.data.values)


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
        # 3 PCA summary tables + the per-sample component scores
        assert sum(p.suffix == ".tsv" for p in d.iterdir()) == 4
        assert (d / f"SLIM_{scope}_pca_scores.tsv").exists()
        # 11 exploratory panels + the scree; no metadata scatter, because
        # grouping.target_col is set but a scope here holds one sample
        assert sum(p.suffix == ".pdf" for p in d.iterdir()) >= 12
        assert (d / f"SLIM_{scope}_pca_scree.pdf").exists()


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
        assert sum(p.suffix == ".tsv" for p in d.iterdir()) == 4
        assert (d / f"SLIM_{scope}_pca_scores.tsv").exists()


def test_exploratory_panels_can_be_switched_off(feature_matrix, synthetic_gmt,
                                                tmp_path):
    """The 11 per-observation panels are skippable; the diagnostics are not.

    Their canvas grows with the number of observations, so on a cohort of
    hundreds of samples they are unreadable as well as slow.  Turning them off
    must still leave the scree and the tables.
    """
    gmt, _ = synthetic_gmt
    mp, mdp = tmp_path / "m.tsv", tmp_path / "md.tsv"
    feature_matrix.data.to_csv(mp, sep="\t")
    feature_matrix.metadata.table.to_csv(mdp, sep="\t")
    out = tmp_path / "no_panels"
    cfg = PipelineConfig(
        input=InputConfig(kind="matrix", matrix_path=str(mp), metadata_path=str(mdp)),
        enrichment=EnrichmentConfig(backend="goslim", geneset="SLIM",
                                    gmt_path=str(gmt), top_n=80),
        grouping=GroupingConfig(mode="target", target_col="cluster_id",
                                scope_col="sample_id"),
        pca=PCAConfig(min_features=5),
        viz=VizConfig(make_figures=True, make_tables=True,
                      exploratory_panels=False),
        output_dir=str(out),
    )
    run(cfg, verbose=False)
    d = out / "SLIM" / "sample_pca" / "A"
    pdfs = {p.name for p in d.iterdir() if p.suffix == ".pdf"}
    # the scree stays; none of the 11 per-observation panels is written
    assert "SLIM_A_pca_scree.pdf" in pdfs
    assert not any("cluster_scores_heatmap" in n or "pairwise" in n
                   or "signatures" in n or "loading" in n for n in pdfs)
    assert (d / "SLIM_A_pca_scores.tsv").exists()


def test_observation_cap_skips_the_panels(feature_matrix, synthetic_gmt, tmp_path):
    """The same skip happens automatically above the observation cap."""
    gmt, _ = synthetic_gmt
    mp, mdp = tmp_path / "m.tsv", tmp_path / "md.tsv"
    feature_matrix.data.to_csv(mp, sep="\t")
    feature_matrix.metadata.table.to_csv(mdp, sep="\t")
    out = tmp_path / "capped"
    cfg = PipelineConfig(
        input=InputConfig(kind="matrix", matrix_path=str(mp), metadata_path=str(mdp)),
        enrichment=EnrichmentConfig(backend="goslim", geneset="SLIM",
                                    gmt_path=str(gmt), top_n=80),
        grouping=GroupingConfig(mode="target", target_col="cluster_id",
                                scope_col="sample_id"),
        pca=PCAConfig(min_features=5),
        viz=VizConfig(make_figures=True, make_tables=False,
                      exploratory_max_observations=1),
        output_dir=str(out),
    )
    run(cfg, verbose=False)
    d = out / "SLIM" / "sample_pca" / "A"
    pdfs = {p.name for p in d.iterdir() if p.suffix == ".pdf"}
    assert "SLIM_A_pca_scree.pdf" in pdfs
    assert not any("cluster_scores_heatmap" in n for n in pdfs)


def test_figures_survive_missing_group_labels(feature_matrix, synthetic_gmt,
                                              tmp_path):
    """A metadata column with gaps must not break the exploratory panels.

    NaN is truthy and unorderable against str, so an unlabelled sample used to
    crash the palette lookup rather than being drawn grey.
    """
    from pathwaytheme.viz.pca_figures import _as_labels
    assert _as_labels([np.nan, "A", None, 1]) == ["", "A", "", "1"]

    gmt, _ = synthetic_gmt
    md = feature_matrix.metadata.table.copy()
    md.loc[md.index[0], "cluster_id"] = np.nan   # one unlabelled observation
    mp, mdp = tmp_path / "m.tsv", tmp_path / "md.tsv"
    feature_matrix.data.to_csv(mp, sep="\t")
    md.to_csv(mdp, sep="\t")
    out = tmp_path / "gappy"
    cfg = PipelineConfig(
        input=InputConfig(kind="matrix", matrix_path=str(mp), metadata_path=str(mdp)),
        enrichment=EnrichmentConfig(backend="goslim", geneset="SLIM",
                                    gmt_path=str(gmt), top_n=80),
        grouping=GroupingConfig(mode="target", target_col="cluster_id"),
        pca=PCAConfig(min_features=5),
        viz=VizConfig(make_figures=True, make_tables=False),
        output_dir=str(out),
    )
    result = run(cfg, verbose=False)
    assert result.pca_results
    written = [p for paths in result.written.values() for p in paths]
    assert any(p.suffix == ".pdf" for p in written)
