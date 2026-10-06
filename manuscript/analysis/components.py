"""Figs. S1-S3 and Table S15: what the principal components are made of.

    figure_s1   the strongest signed loadings of PC1-PC4 in both collections
    figure_s2   the GO:BP loadings as a clustered pathway x component heat map,
                drawn by the package's own loadings panel
    figure_s3   PC1-PC10 of the GO:BP decomposition as a pairs plot coloured by
                sample type -- the exhaustive visual scan that the attribution
                grid of Fig. 1d replaces
    table_s15   the numbers behind Fig. S1

Table S7 takes the pipeline's top loadings by magnitude, which on a component
whose largest loadings are all one sign shows only that sign; Table S15 takes
the top ten at each end separately, which is what Fig. S1 draws.

The full loadings matrix is not among the pipeline's written outputs (only the
top pathways per component are), so the PCA is refitted here from the cached
ssGSEA score matrix with the configuration the run used.  The refit is checked
against the scores the pipeline wrote before anything is drawn, and a refit
that does not reproduce them raises rather than drawing a second, slightly
different decomposition.
"""
from __future__ import annotations

import shutil
from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from pathwaytheme.config import PipelineConfig
from pathwaytheme.pca import run_pca
from pathwaytheme.viz import draw_pc_loadings, draw_pc_pairs

from .paths import CONFIGS, PROC, RESULTS, SCORES, SFIG, STAB, pca_dir

TOP_N = 10
# (run directory, gene-set collection, column label in the figure)
RUNS = [("hallmark_sampletype", "HALLMARK", "50 hallmark processes"),
        ("gobp_sampletype", "GOBP", "GO:BP")]
COMPONENTS = ["PC1", "PC2", "PC3", "PC4"]
TYPE_LABEL = {"CELL_LINE": "Cell line", "TUMOR": "Tumour", "NORMAL": "Normal"}
TYPE_COLOUR = {"Cell line": "#d73027", "Tumour": "#4575b4", "Normal": "#1a9850"}

# where the package's loadings panels are drawn before Fig. S2 is copied out
PANELS = RESULTS / "pca_panels"
RAW_FIG = SFIG / "raw"


def _save(fig, stem: str) -> None:
    SFIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(SFIG / f"{stem}.{ext}", bbox_inches="tight",
                    dpi=200 if ext == "png" else None)


def refit_pca(run: str, geneset: str):
    """Refit a run's PCA from its cached scores; return (result, variance explained).

    Raises if the refit does not reproduce the component scores the run wrote.
    """
    d = pca_dir(run, geneset)
    written = pd.read_csv(d / f"{geneset}___all___pca_scores.tsv",
                          sep="\t", index_col=0)
    sm = pd.read_csv(SCORES / f"{geneset}_ssgsea_scores.tsv", sep="\t", index_col=0)
    # the same samples in the same order the run used, which the written
    # scores record exactly
    sm = sm.loc[:, list(written.index)]
    labels = pd.Series("__all__", index=sm.columns)
    config = PipelineConfig.from_yaml(CONFIGS / f"{run}.yaml").pca
    res = run_pca(sm, labels, config)
    if res is None:
        raise RuntimeError(f"PCA returned nothing for {run}")
    got = res.scores.loc[written.index, list(written.columns)]
    if not np.allclose(got.to_numpy(), written.to_numpy(), atol=1e-6):
        raise RuntimeError(f"refitted PCA does not reproduce {run}'s written scores")
    ve = pd.read_csv(d / f"{geneset}___all___pca_variance_explained.tsv", sep="\t")
    return res, ve["variance_explained"].to_numpy()


def refit_all() -> dict:
    """Run directory -> (PCA result, variance explained), for both collections."""
    return {run: refit_pca(run, geneset) for run, geneset, _ in RUNS}


def _label(term: str) -> str:
    head, _, rest = term.partition("_")
    body = rest if head.isupper() and len(head) > 1 else term
    return body.replace("_", " ").lower()


def table_s15(fits: dict | None = None) -> pd.DataFrame:
    """Table S15: the top ten positive and top ten negative loadings of PC1-PC4,
    both collections -- the numbers Fig. S1 draws, one row per bar."""
    fits = fits or refit_all()
    rows = []
    for run, geneset, _ in RUNS:
        res, ve = fits[run]
        for i, pc in enumerate(COMPONENTS):
            lv = res.loadings.loc[pc].astype(float)
            for axis, sel in (("positive",
                               lv[lv > 0].sort_values(ascending=False).head(TOP_N)),
                              ("negative",
                               lv[lv < 0].sort_values().head(TOP_N))):
                for rank, (path, v) in enumerate(sel.items(), 1):
                    rows.append({"run": run, "geneset": geneset,
                                 "component": pc,
                                 "variance_explained": round(float(ve[i]), 5),
                                 "axis": axis, "rank": rank, "pathway": path,
                                 "label": _label(path),
                                 "loading": round(float(v), 5)})
    out = pd.DataFrame(rows)
    STAB.mkdir(parents=True, exist_ok=True)
    out.to_csv(STAB / "tableS15_pc_loadings_signed.tsv", sep="\t", index=False)
    return out


def figure_s1(fits: dict | None = None) -> plt.Figure:
    """Fig. S1: the strongest signed pathway loadings of PC1-PC4, one column per
    collection (hallmark, GO:BP), one row per component."""
    fits = fits or refit_all()
    fig = plt.figure(figsize=(11.6, 15.2))
    gs = GridSpec(len(COMPONENTS), len(RUNS), figure=fig,
                  hspace=0.20, wspace=0.78,
                  left=0.175, right=0.985, top=0.945, bottom=0.052)
    for j, (run, geneset, column_title) in enumerate(RUNS):
        res, ve = fits[run]
        # a column header, not a title on the first panel: the first panel
        # already carries its component and its variance
        ax0 = fig.add_subplot(gs[0, j])
        pos = ax0.get_position()
        fig.text(pos.x0 + pos.width / 2, 0.968, column_title, fontsize=10,
                 fontweight="bold", ha="center", va="bottom")
        for i, pc in enumerate(COMPONENTS):
            ax = ax0 if i == 0 else fig.add_subplot(gs[i, j])
            draw_pc_loadings(ax, res.loadings.loc[pc], component=pc,
                             top_n=TOP_N, variance_explained=float(ve[i]),
                             wrap=34 if geneset == "GOBP" else 26,
                             label_fontsize=5.6)
    fig.text(0.5, 0.008, "red = the positive end of the component,   "
                         "blue = the negative end", fontsize=8, color="0.35",
             ha="center")
    _save(fig, "figS1_pc_loadings")
    return fig


@contextmanager
def _also_png():
    """Make the package's loadings panel write a PNG beside every PDF it saves."""
    # pathwaytheme exposes no public function for the clustered loadings panel
    # and no PNG option, so its private writer is called and its save patched
    from pathwaytheme.viz import pca_figures
    original = pca_figures._safe_savefig

    def both(fig, path, **kw):
        ok = original(fig, path, **kw)
        if ok and Path(path).suffix == ".pdf":
            original(fig, Path(path).with_suffix(".png"), **kw)
        return ok

    pca_figures._safe_savefig = both
    try:
        yield pca_figures._draw_pc_loadings
    finally:
        pca_figures._safe_savefig = original


def draw_loadings_panels(fits: dict | None = None) -> dict[str, Path]:
    """Draw the package's three loadings panels for each collection into
    ``results/pca_panels/<run>/``: per-component bars, the pathway x component
    heat map in loading order, and the same matrix with the pathways clustered.

    Returns run -> output directory.
    """
    fits = fits or refit_all()
    out = {}
    for run, geneset, column_title in RUNS:
        res, _ = fits[run]
        d = PANELS / run
        d.mkdir(parents=True, exist_ok=True)
        with _also_png() as draw:
            draw(res, column_title, d, geneset, dpi=200)
        if not (d / f"{geneset}_pca_pathway_pc_clustermap.pdf").exists():
            raise RuntimeError(f"the package wrote no clustermap for {run}")
        out[run] = d
    return out


def figure_s2(fits: dict | None = None) -> Path:
    """Fig. S2: the GO:BP loadings as a clustered pathway x component heat map.

    The package's own panel, copied rather than redrawn so it is the panel the
    pipeline would emit.  It reduces the GO:BP collection to the 80 terms with
    the largest loading anywhere and lets the clustering group them.  The
    hallmark collection has 50 pathways in total, all already in Table S7, so
    its copy is not a numbered figure (see :func:`raw_loadings_panels`).

    The panel is saved and closed by the package, so the PNG path is returned
    instead of a Figure.
    """
    panels = draw_loadings_panels(fits)
    src = panels["gobp_sampletype"] / "GOBP_pca_pathway_pc_clustermap.pdf"
    dest = SFIG / "figS2_pc_loadings_clustermap.pdf"
    dest.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png"):
        shutil.copyfile(src.with_suffix(suffix), dest.with_suffix(suffix))
    raw_loadings_panels()
    return dest.with_suffix(".png")


def raw_loadings_panels() -> Path:
    """Not a paper item: copy every loadings panel the package drew, both
    collections, into ``submission/figures/raw/pc_loadings/``.

    Figs. S1 and S2 are composed from these; they are part of the submission
    bundle.  Run after :func:`draw_loadings_panels`.
    """
    keep = RAW_FIG / "pc_loadings"
    keep.mkdir(parents=True, exist_ok=True)
    for run, geneset, _ in RUNS:
        for f in sorted((PANELS / run).glob(f"{geneset}_pca_*")):
            shutil.copyfile(f, keep / f.name)
    return keep


def figure_s3(fits: dict | None = None) -> plt.Figure:
    """Fig. S3: PC1-PC10 of the GO:BP decomposition as a pairs plot coloured by
    sample type."""
    fits = fits or {"gobp_sampletype": refit_pca("gobp_sampletype", "GOBP")}
    res, ve = fits["gobp_sampletype"]
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)
    colour = meta["sample_type"].map(TYPE_LABEL).reindex(res.scores.index)
    colour.name = "sample type"
    fig = plt.figure(figsize=(9.0, 9.0))
    draw_pc_pairs(fig, res.scores, ve, n_components=10, colour_by=colour,
                  colours=TYPE_COLOUR, point_size=1.5,
                  legend_loc=(0.80, 0.82))
    fig.subplots_adjust(left=0.065, right=0.99, top=0.99, bottom=0.06)
    _save(fig, "figS3_pc_pairs")
    return fig


def pc_pairs_by_diagnosis(fits: dict | None = None,
                          min_n: int = 5, size_floor: int = 10) -> list[Path]:
    """Not a paper item: Fig. S3's pairs plot once per diagnosis, that diagnosis
    against every other sample, into ``submission/figures/raw/pc_pairs_by_diagnosis/``.

    A diagnosis is one small level against a large rest, so each gets its own
    copy rather than many colours on one grid.  The one-vs-rest question they
    ask is answered numerically in Tables S8, S9 and S12.  Diagnoses below
    ``min_n`` are skipped (no density is drawn for fewer than five samples);
    those below ``size_floor``, the contrast size floor, are flagged in the
    title.

    The rest is drawn first and the diagnosis on top of it: ``draw_pc_pairs``
    takes its group order from the order the labels first appear, so the rows
    are sorted accordingly.
    """
    fits = fits or {"gobp_sampletype": refit_pca("gobp_sampletype", "GOBP")}
    out_dir = RAW_FIG / "pc_pairs_by_diagnosis"
    out_dir.mkdir(parents=True, exist_ok=True)
    res, ve = fits["gobp_sampletype"]
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)
    dx = meta["diagnosis"].reindex(res.scores.index)
    written = []
    for name, n in dx.value_counts().items():
        if n < min_n:
            continue
        is_dx = (dx == name)
        order = list(dx.index[~is_dx]) + list(dx.index[is_dx])
        colour = pd.Series(np.where(is_dx[order], name, "rest"),
                           index=order, name="diagnosis")
        floor = "" if n >= size_floor else f"   (below the n = {size_floor} floor)"
        fig = plt.figure(figsize=(9.0, 9.0))
        draw_pc_pairs(fig, res.scores.loc[order], ve, n_components=10,
                      colour_by=colour,
                      colours={name: "#c0392b", "rest": "#c8ccd2"},
                      point_size=1.5, legend_loc=(0.80, 0.86),
                      title=f"{name} (n = {n}) against the rest "
                            f"(n = {len(dx) - n}), ten leading GO:BP "
                            f"components{floor}")
        fig.subplots_adjust(left=0.065, right=0.99, top=0.965, bottom=0.06)
        for ext in ("pdf", "png"):
            fig.savefig(out_dir / f"pc_pairs_{name}_vs_rest.{ext}",
                        bbox_inches="tight", dpi=200 if ext == "png" else None)
        plt.close(fig)
        written.append(out_dir / f"pc_pairs_{name}_vs_rest.png")
    return written


def run_all() -> None:
    """Refit both decompositions once, then write Table S15, Figs. S1-S3 and the
    raw loadings and per-diagnosis panels of the submission bundle."""
    fits = refit_all()
    table_s15(fits)
    plt.close(figure_s1(fits))
    figure_s2(fits)
    plt.close(figure_s3(fits))
    pc_pairs_by_diagnosis(fits)
