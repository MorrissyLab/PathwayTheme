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
from .contracts import FeatureMatrix, ScoreMatrix, Grouping, PCAResult, DiffResult
from .io import load_input
from .enrichment import run_enrichment
from .grouping import resolve_grouping
from .pca import run_pca_scopes
from .diff import run_diff
from .categories import load_category_map, summarize_by_category
from .viz import (render_pca_figures, write_pca_tables, render_sanity_heatmap,
                  render_diff_figures, write_diff_table, write_category_table,
                  render_scree, render_pca_scatter,
                  render_metadata_attribution, render_level_attribution,
                  render_theme_heatmap)
from .pca.metadata import (associate_metadata, associate_levels,
                           attribution_matrix)


@dataclass
class PipelineResult:
    feature_matrix: FeatureMatrix
    score_matrix: ScoreMatrix
    grouping: Grouping
    pca_results: list[PCAResult]
    output_dir: Path
    written: dict[str, list[Path]] = field(default_factory=dict)
    diff_result: Optional[DiffResult] = None
    category_summary: Optional[pd.DataFrame] = None


def _apply_sample_filter(fm: FeatureMatrix, cfg) -> FeatureMatrix:
    """Keep only samples whose metadata[filter_col] is in filter_values."""
    from .api import filter_samples
    if not cfg.filter_col or not cfg.filter_values:
        return fm
    return filter_samples(fm, cfg.filter_col, cfg.filter_values)


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
    if config.input.filter_col and config.input.filter_values:
        fm = _apply_sample_filter(fm, config.input)
        log(f"      filtered to {config.input.filter_col} in "
            f"{config.input.filter_values}: {fm.shape[1]} samples")

    # 2 — enrichment
    log(f"[2/5] enrichment (backend={config.enrichment.backend}, "
        f"geneset={config.enrichment.geneset}) ...")
    sm = run_enrichment(fm, config.enrichment)
    log(f"      ScoreMatrix: {sm.shape[0]} terms x {sm.shape[1]} samples")
    if sm.coverage is not None:
        n_req = len(sm.coverage)
        n_drop = int((sm.coverage["status"] != "scored").sum())
        log(f"      gene sets: {n_req - n_drop} of {n_req} scored"
            + (f", {n_drop} dropped for insufficient coverage" if n_drop else ""))

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
            # created here, not by whichever writer happens to run first: with
            # tables off and the exploratory panels skipped, the scree would be
            # the first write into a directory nobody had made
            out_dir.mkdir(parents=True, exist_ok=True)
            prefix = f"{config.enrichment.geneset}_{res.scope}"
            title = f"{config.enrichment.geneset} / {res.scope}"
            paths: list[Path] = []
            if config.viz.make_tables:
                paths += write_pca_tables(res, out_dir, prefix, config.pca)
                # per-sample component scores: needed to relate a component to
                # anything outside the PCA, and not recoverable from the others
                sp = out_dir / f"{prefix}_pca_scores.tsv"
                res.scores.to_csv(sp, sep="\t")
                paths.append(sp)
            if config.viz.make_figures:
                cap = config.viz.exploratory_max_observations
                n_obs = res.scores.shape[0]
                if not config.viz.exploratory_panels:
                    log(f"      {res.scope}: exploratory panels off")
                elif cap and n_obs > cap:
                    # skipping loudly: these panels label every observation, so
                    # at this size they would be unreadable rather than slow
                    log(f"      {res.scope}: {n_obs} observations exceeds "
                        f"viz.exploratory_max_observations={cap} — the 11 "
                        "per-observation panels are skipped, the diagnostics "
                        "below are not")
                else:
                    paths += render_pca_figures(res, out_dir, prefix, title,
                                                config.viz)
            # the scree is a diagnostic, not one of the 11 exploratory panels,
            # so it is emitted whenever figures are on -- it is cheap and it is
            # what decides how many components anything downstream should use
            if config.viz.make_figures:
                paths.append(render_scree(
                    res.variance_explained,
                    out_dir / f"{prefix}_pca_scree.pdf",
                    elbow=config.pca.elbow_components,
                    title=f"{title} — variance per component",
                    dpi=config.viz.dpi))
                if config.grouping.target_col:
                    # reindex, not set_axis: a scope holds a subset of samples,
                    # and a display_col makes scores.index differ from the
                    # sample ids.  Reindexing yields NaN on a mismatch, which the
                    # plot skips, rather than silently pairing the wrong rows.
                    md = sm.metadata.table.reindex(res.scores.index)
                    paths.append(render_pca_scatter(
                        res.scores, res.variance_explained,
                        out_dir / f"{prefix}_pca_scatter_metadata.pdf",
                        colour_by=md.get(config.grouping.target_col),
                        marker_by=(md.get(config.grouping.secondary_col)
                                   if config.grouping.secondary_col else None),
                        group_means=True, dpi=config.viz.dpi,
                        title=f"{title} — colour = {config.grouping.target_col}"))

            # 4b — associate every component with every metadata variable
            if config.metadata.enabled:
                assoc = associate_metadata(
                    res.scores,
                    sm.metadata.table.reindex(res.scores.index),
                    columns=config.metadata.columns or None,
                    min_group_size=config.metadata.min_group_size,
                    alpha=config.metadata.alpha)
                if not assoc.empty:
                    n_sig = int(assoc["significant"].sum())
                    log(f"      {res.scope}: metadata association — "
                        f"{len(assoc)} test(s), {n_sig} at FDR<{config.metadata.alpha}")
                    if config.viz.make_tables:
                        ap = out_dir / f"{prefix}_metadata_pc_association.tsv"
                        assoc.to_csv(ap, sep="\t", index=False)
                        paths.append(ap)
                        mp = out_dir / f"{prefix}_metadata_pc_matrix.tsv"
                        attribution_matrix(assoc).to_csv(mp, sep="\t")
                        paths.append(mp)
                    if config.viz.make_figures:
                        paths.append(render_metadata_attribution(
                            assoc,
                            out_dir / f"{prefix}_metadata_pc_attribution.pdf",
                            variance_explained=res.variance_explained,
                            technical=config.metadata.technical,
                            n_components=config.metadata.plot_components,
                            title=f"{title} — variance attribution",
                            dpi=config.viz.dpi))

            # 4c — which component is the axis of one particular level?  The
            # omnibus test above cannot say; this splits every categorical
            # variable into level-versus-rest.
            if config.metadata.enabled and config.metadata.per_level:
                cols = (config.metadata.per_level_columns
                        or config.metadata.columns or None)
                lv = associate_levels(
                    res.scores,
                    sm.metadata.table.reindex(res.scores.index),
                    columns=cols,
                    min_group_size=config.metadata.min_group_size,
                    alpha=config.metadata.alpha)
                if not lv.empty:
                    n_sig = int(lv["significant"].sum())
                    log(f"      {res.scope}: level association — "
                        f"{len(lv)} test(s), {n_sig} at FDR<{config.metadata.alpha}")
                    if config.viz.make_tables:
                        lp = out_dir / f"{prefix}_metadata_pc_levels.tsv"
                        lv.to_csv(lp, sep="\t", index=False)
                        paths.append(lp)
                    if config.viz.make_figures:
                        paths.append(render_level_attribution(
                            lv,
                            out_dir / f"{prefix}_metadata_pc_levels.pdf",
                            variance_explained=res.variance_explained,
                            n_components=config.metadata.plot_components,
                            title=f"{title} — level vs rest",
                            dpi=config.viz.dpi))
            written[res.scope] = paths
            log(f"      {res.scope}: {len(paths)} file(s)")
    else:
        log("[5/5] viz disabled")

    analysis_root = Path(config.output_dir) / config.enrichment.geneset
    base_prefix = config.enrichment.geneset

    # gene-set coverage: which requested sets were scored, and why not
    if sm.coverage is not None and config.viz.make_tables:
        analysis_root.mkdir(parents=True, exist_ok=True)
        cov_path = analysis_root / f"{base_prefix}_geneset_coverage.tsv"
        sm.coverage.to_csv(cov_path, sep="\t", index=False)
        written.setdefault("_coverage", []).append(cov_path)

    # optional — full-matrix sanity heatmap (QC)
    if config.viz.sanity_heatmap and config.viz.make_figures:
        log("[qc] writing sanity heatmap ...")
        p = render_sanity_heatmap(
            sm, analysis_root / f"{base_prefix}_sanity_heatmap.pdf",
            label_col=config.grouping.target_col,
            max_pathways=config.viz.sanity_max_pathways, dpi=config.viz.dpi)
        if p is not None:
            written.setdefault("_qc", []).append(p)

    # optional — differential pathway analysis (parallel to PCA)
    diff_result: Optional[DiffResult] = None
    category_summary: Optional[pd.DataFrame] = None
    if config.diff.enabled:
        log(f"[diff] differential analysis (method={config.diff.method}) ...")
        diff_result = run_diff(sm, config.diff, grouping)
        n_sig = int((diff_result.table["fdr"] < 0.05).sum()) if not diff_result.table.empty else 0
        log(f"       {len(diff_result.comparisons)} comparison(s), {n_sig} pathway(s) FDR<0.05")
        dpaths: list[Path] = []
        if config.viz.make_tables:
            dpaths += write_diff_table(diff_result, analysis_root, base_prefix)
        if config.viz.make_figures:
            dpaths += render_diff_figures(diff_result, analysis_root, base_prefix,
                                          top_n=config.diff.top_n, dpi=config.viz.dpi)

        # optional — roll differential results up into broad categories
        if config.categories.enabled and config.categories.map_path:
            log("[categories] summarizing by category ...")
            mapping = load_category_map(
                config.categories.map_path, key_col=config.categories.key_col,
                category_col=config.categories.category_col)
            sig_col = (config.categories.significance_col
                       if config.categories.split_significance else None)
            category_summary = summarize_by_category(
                diff_result.table, mapping,
                key_col="pathway", value_col="effect",
                group_cols=["comparison"], stat=config.categories.stat,
                unmapped_label=config.categories.unmapped_label,
                significance_col=sig_col, alpha=config.categories.alpha)
            if config.viz.make_tables:
                dpaths += write_category_table(category_summary, analysis_root, base_prefix)
            if config.viz.make_figures:
                dpaths.append(render_theme_heatmap(
                    category_summary,
                    analysis_root / f"{base_prefix}_theme_heatmap.pdf",
                    top_n=config.categories.plot_top_themes,
                    title=f"{base_prefix} — GO-slim themes",
                    dpi=config.viz.dpi))
        written.setdefault("_diff", []).extend(dpaths)

    return PipelineResult(fm, sm, grouping, results, Path(config.output_dir),
                          written, diff_result=diff_result,
                          category_summary=category_summary)
