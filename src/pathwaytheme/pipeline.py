"""End-to-end orchestrator: io -> enrichment -> grouping -> pca -> viz.

``run(config)`` executes all stages and writes, per PCA scope, the 11 figures
and 3 tables into ``<output_dir>/<geneset>/sample_pca/<scope>/``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from .config import PipelineConfig
from .contracts import FeatureMatrix, ScoreMatrix, Grouping, PCAResult
from .io import load_input
from .enrichment import run_enrichment
from .grouping import resolve_grouping
from .pca import run_pca_scopes
from .viz import render_pca_figures, write_pca_tables


@dataclass
class PipelineResult:
    feature_matrix: FeatureMatrix
    score_matrix: ScoreMatrix
    grouping: Grouping
    pca_results: list[PCAResult]
    output_dir: Path
    written: dict[str, list[Path]] = field(default_factory=dict)


def _build_display_labels(sm: ScoreMatrix, grouping_cfg) -> Optional[pd.Series]:
    """Observation display labels: explicit column, else `cl{cluster_id}`, else raw id."""
    meta = sm.metadata.table
    if grouping_cfg.display_col and grouping_cfg.display_col in meta.columns:
        return meta[grouping_cfg.display_col].astype(str)
    if "cluster_id" in meta.columns:
        return ("cl" + meta["cluster_id"].astype(str))
    return None


def run(config: PipelineConfig, *, verbose: bool = True) -> PipelineResult:
    def log(msg: str):
        if verbose:
            print(msg, flush=True)

    # 1 — input
    log(f"[1/5] loading input (kind={config.input.kind}) ...")
    fm = load_input(config.input)
    log(f"      FeatureMatrix: {fm.shape[0]} features x {fm.shape[1]} samples")

    # 2 — enrichment
    log(f"[2/5] enrichment (backend={config.enrichment.backend}, "
        f"geneset={config.enrichment.geneset}) ...")
    sm = run_enrichment(fm, config.enrichment)
    log(f"      ScoreMatrix: {sm.shape[0]} terms x {sm.shape[1]} samples")

    # 3 — grouping
    log(f"[3/5] grouping (mode={config.grouping.mode}) ...")
    grouping = resolve_grouping(sm, config.grouping)
    scopes = grouping.scopes()
    log(f"      {grouping.labels.nunique()} labels, {len(scopes)} scope(s)")

    # 4 — PCA per scope
    log("[4/5] PCA per scope ...")
    display = _build_display_labels(sm, config.grouping)
    # display labels are keyed by sample id; run_pca_scopes expects a metadata
    # column name, so stash the computed series on the metadata table.
    disp_col = None
    if display is not None:
        sm.metadata.table = sm.metadata.table.assign(__display__=display.reindex(sm.samples))
        disp_col = "__display__"
    results = run_pca_scopes(
        sm, grouping, config.pca,
        display_col=disp_col,
        sublabel_col=config.grouping.secondary_col,
        size_col=config.grouping.size_col,
    )
    log(f"      {len(results)} PCA result(s) (scopes with >= "
        f"{config.pca.min_observations} observations)")

    # 5 — figures + tables
    out_root = Path(config.output_dir) / config.enrichment.geneset / "sample_pca"
    written: dict[str, list[Path]] = {}
    if config.viz.make_figures or config.viz.make_tables:
        log("[5/5] writing figures + tables ...")
        for res in results:
            out_dir = out_root / res.scope
            prefix = f"{config.enrichment.geneset}_{res.scope}"
            title = f"{config.enrichment.geneset} / {res.scope}"
            paths: list[Path] = []
            if config.viz.make_tables:
                paths += write_pca_tables(res, out_dir, prefix, config.pca)
            if config.viz.make_figures:
                paths += render_pca_figures(res, out_dir, prefix, title, config.viz)
            written[res.scope] = paths
            log(f"      {res.scope}: {len(paths)} file(s)")
    else:
        log("[5/5] viz disabled")

    return PipelineResult(fm, sm, grouping, results, Path(config.output_dir), written)
