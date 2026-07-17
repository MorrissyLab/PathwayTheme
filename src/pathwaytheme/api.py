"""Simple high-level API — one call per stage, or one call for everything.

    import pathwaytheme as pt

    # ---- full pipeline ----
    pt.run_pipeline("examples/moh_sm.yaml")

    # ---- or stage by stage ----
    fm      = pt.load_matrix("expr.tsv", metadata="meta.tsv")
    scores  = pt.enrich(fm, backend="ssgsea", gmt="go_bp.gmt", geneset="GO_BP")
    groups  = pt.group(scores, mode="target", target_col="grp", scope_col="sample_id")
    results = pt.pca(scores, groups, size_col="n_cells")
    pt.figures(results, "out", geneset="GO_BP")

Every wrapper builds the underlying config for you; drop down to
``pathwaytheme.config`` + the stage modules only when you need full control.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd

from .config import (PipelineConfig, InputConfig, EnrichmentConfig,
                     GroupingConfig, PCAConfig, DiffConfig, VizConfig)
from .contracts import (FeatureMatrix, ScoreMatrix, Grouping, PCAResult,
                        DiffResult, SampleMetadata)
from .io import load_input
from .enrichment import run_enrichment
from .grouping import resolve_grouping
from .pca import run_pca_scopes
from .diff import run_diff
from .categories import load_category_map, summarize_by_category
from .viz import (render_pca_figures, write_pca_tables, render_sanity_heatmap,
                  render_diff_figures, write_diff_table, write_category_table)


# ── stage 1: input ────────────────────────────────────────────────────────
def load_matrix(matrix_path: str, metadata: Optional[str] = None,
                sep: str = "\t") -> FeatureMatrix:
    """Load a pre-quantified features x samples table (+ optional metadata)."""
    return load_input(InputConfig(kind="matrix", matrix_path=matrix_path,
                                  metadata_path=metadata, sep=sep))


def load_h5ad(directory: Optional[str] = None, files: Optional[list[str]] = None,
              cluster_col: str = "seurat_clusters", sample_col: str = "Sample_name",
              aggregate: str = "mean") -> FeatureMatrix:
    """Aggregate single-cell .h5ad files into per-cluster pseudobulk."""
    return load_input(InputConfig(kind="h5ad", h5ad_dir=directory,
                                  h5ad_paths=files or [], cluster_col=cluster_col,
                                  sample_col=sample_col, aggregate=aggregate))


# ── optional: select samples by a metadata value ──────────────────────────
def filter_samples(obj: Union[FeatureMatrix, ScoreMatrix], column: str,
                   keep) -> Union[FeatureMatrix, ScoreMatrix]:
    """Keep only samples whose ``metadata[column]`` is in ``keep``.

    Works on a FeatureMatrix (before enrichment) or a ScoreMatrix (after), e.g.
    ``fm = pt.filter_samples(fm, "OncoTree_Code", ["ASPS"])`` to run the rest of
    the pipeline on one tumour type.
    """
    keep = {str(v) for v in ([keep] if isinstance(keep, str) else keep)}
    meta = obj.metadata.table
    if column not in meta.columns:
        raise KeyError(f"metadata column {column!r} not found; available: "
                       f"{list(meta.columns)}")
    cols = [s for s in obj.data.columns if str(meta.loc[s, column]) in keep]
    if not cols:
        raise ValueError(f"no samples have {column} in {sorted(keep)}")
    sub = obj.data.loc[:, cols]
    new_meta = SampleMetadata(meta.loc[cols])
    if isinstance(obj, FeatureMatrix):
        return FeatureMatrix(sub, new_meta)
    return ScoreMatrix(sub, new_meta, backend=obj.backend, geneset=obj.geneset)


# ── stage 2: enrichment ────────────────────────────────────────────────────
def enrich(fm: FeatureMatrix, backend: str = "ssgsea", *, geneset: str = "geneset",
           gmt: Optional[str] = None, cache_dir: Optional[str] = None,
           **kwargs) -> ScoreMatrix:
    """Run an enrichment backend -> ScoreMatrix (terms x samples).

    ``backend`` is "ssgsea" | "enrichr" | "goslim".  Extra keyword args map to
    :class:`~pathwaytheme.config.EnrichmentConfig` fields (e.g. ``top_n=200``,
    ``threads=8``, ``enrichr_library=...``, ``goslim_score="count"``).
    """
    cfg = EnrichmentConfig(backend=backend, geneset=geneset, gmt_path=gmt,
                           cache_dir=cache_dir, **kwargs)
    return run_enrichment(fm, cfg)


# ── stage 3: grouping ──────────────────────────────────────────────────────
def group(sm: ScoreMatrix, mode: str = "target", *, target_col: Optional[str] = None,
          scope_col: Optional[str] = None, **kwargs) -> Grouping:
    """Resolve a comparison Grouping.

    ``mode`` is "target" | "existing" | "auto".  Extra kwargs map to
    :class:`~pathwaytheme.config.GroupingConfig` (e.g. ``n_clusters=4``,
    ``method="hierarchical"``).
    """
    cfg = GroupingConfig(mode=mode, target_col=target_col, scope_col=scope_col, **kwargs)
    return resolve_grouping(sm, cfg)


# ── stage 4: PCA ───────────────────────────────────────────────────────────
def pca(sm: ScoreMatrix, grouping: Optional[Grouping] = None, *,
        display_col: Optional[str] = None, secondary_col: Optional[str] = None,
        size_col: Optional[str] = None, **kwargs) -> list[PCAResult]:
    """Run one PCA per scope.  If ``grouping`` is omitted, all samples form a
    single unlabelled scope.  Extra kwargs map to
    :class:`~pathwaytheme.config.PCAConfig` (e.g. ``k_max=8``)."""
    if grouping is None:
        grouping = Grouping(labels=pd.Series("all", index=sm.samples, name="group"))
    return run_pca_scopes(sm, grouping, PCAConfig(**kwargs),
                          display_col=display_col, sublabel_col=secondary_col,
                          size_col=size_col)


# ── optional stage: differential pathway analysis ─────────────────────────
def diff(sm: ScoreMatrix, grouping: Optional[Grouping] = None, *,
         method: str = "welch", group_col: Optional[str] = None,
         reference: Optional[str] = None, contrasts=None,
         **kwargs) -> DiffResult:
    """Compare pathway scores between groups -> DiffResult (per-pathway p/FDR).

    ``method`` is "welch" | "mannwhitney" | "moderated_t".  Groups come from
    ``group_col`` (a metadata column) or from ``grouping``.  With no
    ``reference``/``contrasts`` this runs one-vs-rest for every group; pass
    ``reference="X"`` for every-group-vs-X or ``contrasts=[["A","B"], ...]``.
    """
    cfg = DiffConfig(enabled=True, method=method, group_col=group_col,
                     reference=reference, contrasts=contrasts, **kwargs)
    return run_diff(sm, cfg, grouping)


def summarize_categories(table, mapping, *, key_col: str = "pathway",
                         value_col: str = "effect",
                         group_cols: Optional[list] = None,
                         stat: str = "mean"):
    """Roll a result table up into broad categories.

    ``table`` is a DiffResult, its ``.table``, or any long DataFrame with a
    pathway column.  ``mapping`` is a term->category dict/Series or a path to a
    TSV (loaded via :func:`load_category_map`).
    """
    if isinstance(table, DiffResult):
        table = table.table
    if isinstance(mapping, (str, Path)):
        mapping = load_category_map(mapping, key_col=key_col)
    default_groups = ["comparison"] if "comparison" in table.columns else None
    return summarize_by_category(
        table, mapping, key_col=key_col, value_col=value_col,
        group_cols=group_cols if group_cols is not None else default_groups,
        stat=stat)


def sanity_heatmap(sm: ScoreMatrix, out_path: str, *, label_col: Optional[str] = None,
                   max_pathways: int = 200, dpi: int = 120):
    """Write a z-scored pathway x sample QC heatmap of the whole score matrix."""
    return render_sanity_heatmap(sm, out_path, label_col=label_col,
                                 max_pathways=max_pathways, dpi=dpi)


# ── stage 5: figures + tables ──────────────────────────────────────────────
def figures(results: Union[PCAResult, list[PCAResult]], out_dir: str, *,
            geneset: str = "analysis", tables: bool = True, figures: bool = True,
            dpi: int = 120) -> dict[str, list[Path]]:
    """Write the 11 PDFs + 3 TSVs per PCA result into
    ``<out_dir>/<geneset>/sample_pca/<scope>/``."""
    if isinstance(results, PCAResult):
        results = [results]
    viz = VizConfig(make_figures=figures, make_tables=tables, dpi=dpi)
    pcfg = PCAConfig()
    root = Path(out_dir) / geneset / "sample_pca"
    written: dict[str, list[Path]] = {}
    for res in results:
        d = root / res.scope
        prefix = f"{geneset}_{res.scope}"
        title = f"{geneset} / {res.scope}"
        paths: list[Path] = []
        if tables:
            paths += write_pca_tables(res, d, prefix, pcfg)
        if figures:
            paths += render_pca_figures(res, d, prefix, title, viz)
        written[res.scope] = paths
    return written


# ── all-in-one ─────────────────────────────────────────────────────────────
def run_pipeline(config: Union[str, Path, PipelineConfig, None] = None, **overrides):
    """Run the full pipeline.

    ``config`` may be a YAML path, a :class:`PipelineConfig`, or omitted (build
    entirely from ``overrides``).  ``overrides`` are dotted keys, e.g.
    ``run_pipeline(**{"input.kind": "matrix", "input.matrix_path": "x.tsv"})``.
    """
    if isinstance(config, (str, Path)):
        cfg = PipelineConfig.from_yaml(config)
    elif isinstance(config, PipelineConfig):
        cfg = config
    else:
        cfg = PipelineConfig()
    for dotted, value in overrides.items():
        section, _, field = dotted.partition(".")
        target = getattr(cfg, section) if field else cfg
        setattr(target, field or section, value)
    from .pipeline import run
    return run(cfg)
