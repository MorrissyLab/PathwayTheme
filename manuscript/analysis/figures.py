"""Figure 1, Figs. S5-S7 and Tables S1-S10 and S20.

Every panel is drawn by ``pathwaytheme.viz``, the same functions ``pathwaytheme
run`` uses to write one file per plot; what lives here is only the arrangement:
which panels go in Fig. 1, how they are laid out and lettered, and which of the
pipeline's outputs become supplementary tables.  Nothing is recomputed: the
inputs are the pipeline's outputs under ``results/pipeline``, the processed
cohort tables under ``data/processed``, and the summary tables of
:mod:`analysis.summaries`.

:func:`analysis.summaries.run_all` must run before the table functions here:
Tables S5, S6 and S7 are copies of the files it writes into ``results/tables``.

Writes into ``submission/``:

    figures/figure1, figS5_themes_pairwise, figS6_matched_vs_unmatched,
    figS7_themes_by_diagnosis       .pdf and .png
    figures/raw/pipeline_<run>/     every plot each pipeline run drew
    tables/tableS1..S10, tableS20   .tsv
    configs/*.yaml                  the run configurations
"""
from __future__ import annotations

import glob
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from pathwaytheme.config import PipelineConfig
from pathwaytheme.viz import (draw_workflow, draw_scree, draw_pca_scatter,
                              draw_metadata_attribution, draw_theme_heatmap,
                              draw_theme_reduction, draw_effect_concordance)

from .paths import CONFIGS, PIPE, PROC, SFIG, STAB, SUB, TABLES

# one directory per experiment behind a figure: the plots a reader would need
# to trace a panel back to what produced it
RAW_FIG = SFIG / "raw"

ALPHA = 0.05
# the smallest group a contrast is read at; the differential stage applies the
# same floor
SIZE_FLOOR = 10

TN = "TUMOR_vs_NORMAL"
CN = "CELL_LINE_vs_NORMAL"
CT = "CELL_LINE_vs_TUMOR"
MATCHED = ["ES_Tumor_vs_ES_CellLine", "NBL_Tumor_vs_NBL_CellLine",
           "RMS_Tumor_vs_RMS_CellLine"]
MATCHED_SHORT = ["Ewing\ntumour vs line", "neuroblastoma\ntumour vs line",
                 "rhabdomyosarcoma\ntumour vs line"]
MATCHED_NAME = ["Ewing sarcoma", "neuroblastoma", "rhabdomyosarcoma"]

# The stored values are the delivered vocabulary; a figure carries the words
# the paper uses.
TYPE_LABEL = {"CELL_LINE": "Cell line", "TUMOR": "Tumour", "NORMAL": "Normal"}
VAR_LABEL = {"sample_type": "sample type", "diagnosis": "diagnosis",
             "library_type": "library chemistry", "rin": "RIN",
             "seq_kit": "sequencing kit", "zero_fraction": "zero fraction"}
TYPE_COLOUR = {"Cell line": "#d73027", "Tumour": "#4575b4",
               "Normal": "#1a9850"}
TECHNICAL = ("library_type", "rin", "seq_kit", "zero_fraction")
VAR_ORDER = ["sample_type", "diagnosis", "library_type", "rin", "seq_kit",
             "zero_fraction"]


# ── the pipeline's outputs, located by name ────────────────────────────────
def _one(run: str, pattern: str) -> Path:
    hits = glob.glob(str(PIPE / run / "**" / pattern), recursive=True)
    if not hits:
        raise FileNotFoundError(f"{run}: no {pattern} — run the pipeline first")
    return Path(hits[0])


def scores(run: str) -> pd.DataFrame:
    return pd.read_csv(_one(run, "*_pca_scores.tsv"), sep="\t", index_col=0)


def variance(run: str) -> np.ndarray:
    v = pd.read_csv(_one(run, "*_pca_variance_explained.tsv"), sep="\t")
    return v["variance_explained"].to_numpy(dtype=float)


def assoc(run: str) -> pd.DataFrame:
    return pd.read_csv(_one(run, "*_metadata_pc_association.tsv"), sep="\t")


def levels(run: str) -> pd.DataFrame:
    return pd.read_csv(_one(run, "*_metadata_pc_levels.tsv"), sep="\t")


def categories(run: str) -> pd.DataFrame:
    return pd.read_csv(_one(run, "*_category_summary.tsv"), sep="\t")


def differential(run: str, contrast: str | None = None) -> pd.DataFrame:
    d = pd.read_csv(_one(run, "*_differential_*.tsv"), sep="\t")
    return d if contrast is None else d[d.comparison == contrast]


def meta() -> pd.DataFrame:
    return pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)


def elbow(run: str) -> int:
    """The number of components the run's configuration marks as the elbow."""
    return PipelineConfig.from_yaml(CONFIGS / f"{run}.yaml").pca.elbow_components


def group_sizes(run: str) -> dict[str, tuple[int, int]]:
    """Contrast -> (n case, n reference), as the differential stage recorded them."""
    g = differential(run).groupby("comparison")[["n_case", "n_reference"]].first()
    return {c: (int(r.n_case), int(r.n_reference)) for c, r in g.iterrows()}


def diagnosis_contrasts(floor: int = SIZE_FLOOR) -> list[tuple[str, int]]:
    """The one-vs-rest diagnosis contrasts that clear the size floor, largest first.

    Returns (contrast name, n case) pairs.
    """
    d = differential("gobp_diagnosis")
    g = (d.groupby("comparison")["n_case"].first()
         .sort_values(ascending=False))
    g = g[g >= floor]
    return [(c, int(n)) for c, n in g.items()]


def _label(ax, letter: str, dx: float = -0.16, dy: float = 1.06) -> None:
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left")


def _save(fig, path: Path, dpi: int = 300) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=dpi, bbox_inches="tight")


# ── components and diagnoses ───────────────────────────────────────────────
def diagnosis_components(run: str, floor: int = SIZE_FLOOR) -> pd.DataFrame:
    """For each diagnosis, the component that separates it best.

    Read from the level-versus-rest grid the run itself wrote: every
    (component, diagnosis) pair carries a one-vs-rest rank-biserial r and its
    AUC, and the row kept per diagnosis is the one with the largest |r|.
    Diagnoses below the ``floor`` group size are dropped, the same floor the
    differential stage applies.

    This is what Fig. 1b marks and what Table S20 lists, so neither can drift
    from Tables S8 and S9.
    """
    lv = levels(run)
    d = lv[(lv.variable == "diagnosis") & (lv.n_level >= floor)].copy()
    d["mag"] = d["effect"].abs()
    best = d.sort_values("mag", ascending=False).groupby("level").head(1)
    best = best.assign(pc=best["component"].str.replace("PC", "").astype(int))
    return (best.rename(columns={"level": "diagnosis", "n_level": "n"})
            [["diagnosis", "n", "component", "pc", "effect", "auc", "fdr",
              "direction", "mag"]]
            .sort_values(["pc", "mag"], ascending=[True, False])
            .reset_index(drop=True))


def named_components(run: str, floor: int = SIZE_FLOOR) -> dict[int, str]:
    """Component number -> the diagnoses that separate best on it.

    A scree plot read alone says "keep the first few".  The diagnoses peak on
    components well beyond the elbow, so the scree panel names them rather than
    leaving the reader to infer from a bar height that nothing is out there.
    Where more than one diagnosis peaks on the same component all of them are
    named, strongest first.
    """
    d = diagnosis_components(run, floor)
    return {int(pc): ", ".join(b["diagnosis"])
            for pc, b in d.groupby("pc", sort=True)}


# ── figure 1 ──────────────────────────────────────────────────────────────
def figure1() -> plt.Figure:
    """Figure 1: the pipeline (a), scree (b), PC1-PC2 scatter (c), attribution grid (d)
    and tumour-vs-normal themes (e).

    Panels b-d are the hallmark decomposition of the 901-sample cohort; e is
    the GO:BP tumour-vs-normal contrast at term and theme grain.
    """
    md = meta()
    hs, hv = scores("hallmark_sampletype"), variance("hallmark_sampletype")
    m = md.reindex(hs.index)
    # authored at the width of a journal page.  Three gridspecs rather than one
    # because the three bands want different margins: the schematic is full
    # bleed, the scree and the scatter are a pair, and the attribution grid
    # sits in a band of its own.
    fig = plt.figure(figsize=(7.4, 8.1))

    # (a) the pipeline as a band across the top: the two branches sit side by
    # side under the score matrix they share
    band = fig.add_gridspec(1, 1, left=0.010, right=0.995,
                            top=0.928, bottom=0.760)
    ax = fig.add_subplot(band[0])
    draw_workflow(ax)
    _label(ax, "a", dx=0.002, dy=1.10)

    # the gap under this row is wide because the scree's diagnosis names are
    # set vertically beneath its axis: they live in the gap, not in the panel
    pair = fig.add_gridspec(1, 2, left=0.115, right=0.975,
                            top=0.690, bottom=0.470, wspace=0.60)

    # (b) how much variance, and how far out the diagnoses reach
    ax = fig.add_subplot(pair[0])
    k = elbow("hallmark_sampletype")
    marked = named_components("hallmark_sampletype")
    n_diag = sum(len(v.split(", ")) for v in marked.values())
    draw_scree(ax, hv, elbow=k, mark=marked,
               mark_label="separates a diagnosis best",
               title=f"variance concentrates early — PC1–PC{k} carry "
                     f"{hv[:k].sum() * 100:.1f}%,\n"
                     f"but the {n_diag} diagnoses peak on {len(marked)} "
                     f"different components")
    _label(ax, "b", dx=-0.20)

    # (c) the biological axis
    ax = fig.add_subplot(pair[1])
    draw_pca_scatter(ax, hs, hv, colour_by=m["sample_type"].map(TYPE_LABEL),
                     colours=TYPE_COLOUR, group_means=True,
                     legend_loc="lower right",
                     title="PC1: cell line < tumour < normal")
    ax.set_xlabel(f"{ax.get_xlabel()}  — malignancy axis", fontsize=8)
    _label(ax, "c", dx=-0.20)

    # (d) what every component is: the test that assigns each one to a
    # biological or technical variable
    grid = fig.add_gridspec(1, 2, left=0.115, right=0.975,
                            top=0.285, bottom=0.110, wspace=0.55,
                            width_ratios=[1.10, 1.00])
    ax = fig.add_subplot(grid[0])
    draw_metadata_attribution(
        ax, assoc("hallmark_sampletype"), variance_explained=hv,
        technical=TECHNICAL, n_components=6, column_order=VAR_ORDER,
        title="what each component corresponds to\n"
              "$\\eta^2$ / $\\rho^2$, * FDR<0.05; red = technical")
    ax.set_xticklabels([VAR_LABEL[v] for v in VAR_ORDER])
    _label(ax, "d", dx=-0.235, dy=1.12)

    # (e) the supervised branch that panel (a) promises.  Fig. S5 draws the
    # themes of all three pairwise contrasts; what is shown here is the
    # reduction itself, both grains of one contrast on one axis.
    ax = fig.add_subplot(grid[1])
    d_tn = differential("gobp_tumor_normal", TN)
    cat_tn = categories("gobp_tumor_normal")
    sig = cat_tn[(cat_tn.comparison == TN)
                 & (cat_tn.significance == "significant")
                 & (cat_tn.category != "unmapped")].set_index("category")
    draw_theme_reduction(
        ax, d_tn.loc[d_tn.fdr < ALPHA, "effect"], sig["mean_effect"],
        theme_sizes=sig["n_pathways"], annotate=2, label_chars=13,
        title="one contrast at both grains\n"
              "tumour vs normal; marker area = terms in the theme")
    _label(ax, "e", dx=-0.30, dy=1.12)

    fig.suptitle(f"PathwayTheme on {len(hs)} paediatric solid-tumour, normal and "
                 "cell-line transcriptomes:\nbiological and technical axes "
                 "separated and named, contrasts reduced to themes",
                 fontsize=10, y=0.988)
    _save(fig, SFIG / "figure1")
    return fig


# ── supplementary figures ──────────────────────────────────────────────────
def _themes_figure(run: str, stem: str, title: str,
                   comparisons: list[str], labels: list[str],
                   row_order: str, col_width: float = 1.65) -> plt.Figure:
    """A run's contrasts as GO-slim themes, every theme that has one.

    Every theme is drawn, not a top few: the summarisation exists to turn
    thousands of terms into a number a reader can hold, and showing a subset of
    the themes would put that out of reach again.  ``comparisons`` fixes which
    contrasts appear and in what order, and ``labels`` sets them in the paper's
    words rather than the delivered vocabulary.
    """
    c = categories(run)
    sig = c[(c.significance == "significant") & (c.category != "unmapped")]
    sig = sig[sig.comparison.isin(comparisons)]
    c = c[c.comparison.isin(comparisons)]
    n_col = len(comparisons)
    n_row = sig.category.nunique() or 1
    fig, ax = plt.subplots(figsize=(2.9 + col_width * n_col,
                                    1.5 + 0.135 * n_row))
    draw_theme_heatmap(ax, c, top_n=None, comparisons=comparisons,
                       row_order=row_order, fontsize=5.6, label_chars=60,
                       title=f"{title}\nall {n_row} GO-slim themes with a "
                             "significant term in any contrast")
    ax.set_xticklabels(labels, fontsize=6.4)
    ax.tick_params(axis="x", length=0)
    fig.tight_layout()
    _save(fig, SFIG / stem, dpi=200)
    return fig


def figure_s5() -> plt.Figure:
    """Fig. S5: GO-slim themes of the three pairwise sample-type contrasts.

    Ordered along the malignancy axis (tumour vs normal, cell line vs normal,
    cell line vs tumour) rather than alphabetically, with each contrast's group
    sizes under its name.
    """
    n = group_sizes("gobp_tumor_normal")
    names = {TN: "tumour\nvs normal", CN: "cell line\nvs normal",
             CT: "cell line\nvs tumour"}
    order = [TN, CN, CT]
    return _themes_figure(
        "gobp_tumor_normal", "figS5_themes_pairwise",
        "Fig. S5  Tumour, normal and cell line, pairwise",
        comparisons=order,
        labels=[f"{names[c]}\n{n[c][0]} vs {n[c][1]}" for c in order],
        row_order="signed")


def figure_s6() -> plt.Figure:
    """Fig. S6: the three matched tumour-vs-cell-line contrasts (a) and each against
    the pooled, diagnosis-blind contrast (b-d), at theme grain.
    """
    def wide(c):
        c = c[(c.significance == "significant") & (c.category != "unmapped")]
        return c.pivot_table(index="category", columns="comparison",
                             values="mean_effect")
    cat_m = categories("gobp_matched")
    u, m = wide(categories("gobp_tumor_normal")), wide(cat_m)
    n_heat = int(m.reindex(columns=MATCHED).dropna(how="all").shape[0])
    d_m = differential("gobp_matched")
    n_terms = int(d_m.pathway.nunique())
    n_themes = int(cat_m.loc[cat_m.category != "unmapped", "category"].nunique())

    # every theme is drawn, so the heat map is a tall column and (b-d) stack
    # beside it
    fig = plt.figure(figsize=(11.4, 10.2))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[1.02, 1.0],
                  wspace=0.42, hspace=0.60,
                  left=0.115, right=0.975, top=0.885, bottom=0.055)

    ax = fig.add_subplot(gs[:, 0])
    draw_theme_heatmap(ax, cat_m, top_n=None,
                       comparisons=MATCHED, colorbar_horizontal=True,
                       row_order="signed", fontsize=5.6, label_chars=60,
                       title="the three matched contrasts as GO-slim themes\n"
                             f"all {n_heat}; colour = mean effect, text = terms")
    ax.set_xticklabels(MATCHED_SHORT, fontsize=6.0)
    _label(ax, "a", dx=-0.52, dy=1.02)

    # the group sizes go on the axes: "matched" and "pooled" name a design, and
    # a reader cannot tell which samples either one used from the words alone
    sizes = {**group_sizes("gobp_tumor_normal"), **group_sizes("gobp_matched")}
    n_line, n_tum = sizes[CT]

    for i, (mc, name) in enumerate(zip(MATCHED, MATCHED_NAME)):
        ax = fig.add_subplot(gs[i, 1])
        ax.set_anchor("N")
        if mc not in m.columns or CT not in u.columns:
            ax.set_axis_off()
            continue
        # sign: the pooled contrast is cell line vs tumour, the matched ones are
        # tumour vs cell line, so one side is negated to put them on one axis
        n_case, n_ref = sizes[mc]
        draw_effect_concordance(
            ax, -u[CT], m[mc], title=name, unit="themes",
            poles=("higher in the cell lines", "higher in the tumours"),
            colour_by_agreement=True, annotate=4,
            annotate_disagreements=True, label_chars=22,
            xlabel=f"pooled: every diagnosis together\n"
                   f"all {n_tum} tumours vs all {n_line} cell lines",
            ylabel=f"matched: {name} only\n"
                   f"{n_case} tumours vs {n_ref} cell lines")
        _label(ax, "bcd"[i], dx=-0.26, dy=1.17)

    fig.suptitle("Fig. S6  Three contrasts that share no samples agree with "
                 "each other (a) and with the contrast that ignores the "
                 "diagnosis (b–d)\n"
                 f"a theme is a GO-slim category; a contrast's {n_terms:,} GO:BP "
                 f"terms fall into {n_themes} of them, and its effect is the mean "
                 "over the terms significant in that contrast",
                 fontsize=10, y=0.992,
                 linespacing=1.5)
    _save(fig, SFIG / "figS6_matched_vs_unmatched", dpi=200)
    return fig


def figure_s7() -> plt.Figure:
    """Fig. S7: GO-slim themes of each diagnosis against the other tumours.

    Held to the diagnoses that clear the contrast size floor, largest first.
    """
    dx = diagnosis_contrasts()
    return _themes_figure(
        "gobp_diagnosis", "figS7_themes_by_diagnosis",
        "Fig. S7  One-vs-rest per diagnosis, tumours only",
        comparisons=[c for c, _ in dx], row_order="magnitude",
        labels=[f"{c.replace('_vs_rest', '')}\nvs rest\n{n}" for c, n in dx],
        col_width=0.95)


# ── supplementary tables ──────────────────────────────────────────────────
def _copy_table(src: Path, name: str) -> pd.DataFrame:
    """Copy ``src`` byte for byte into ``submission/tables/name``."""
    STAB.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, STAB / name)
    return pd.read_csv(STAB / name, sep="\t")


def table_s1() -> pd.DataFrame:
    """Table S1: the analysed cohort by diagnosis and sample type."""
    return _copy_table(PROC / "cohort_summary.tsv", "tableS1_cohort.tsv")


def table_s2() -> pd.DataFrame:
    """Table S2: provenance of every annotation column: source, coverage and completions."""
    return _copy_table(PROC / "provenance.tsv", "tableS2_provenance.tsv")


def table_s3() -> pd.DataFrame:
    """Table S3: metadata-by-component association tests, hallmark collection."""
    return _copy_table(_one("hallmark_sampletype", "*_metadata_pc_association.tsv"),
                       "tableS3_metadata_pc_HALLMARK.tsv")


def table_s4() -> pd.DataFrame:
    """Table S4: metadata-by-component association tests, GO:BP collection."""
    return _copy_table(_one("gobp_sampletype", "*_metadata_pc_association.tsv"),
                       "tableS4_metadata_pc_GOBP.tsv")


def table_s5() -> pd.DataFrame:
    """Table S5: terms tested, significant, up and down per contrast."""
    return _copy_table(TABLES / "contrast_counts.tsv", "tableS5_contrast_counts.tsv")


def table_s6() -> pd.DataFrame:
    """Table S6: every GO-slim theme per contrast, significant and not."""
    return _copy_table(TABLES / "theme_summary_full.tsv", "tableS6_theme_summary.tsv")


def table_s7() -> pd.DataFrame:
    """Table S7: the strongest positive and negative pathways per component."""
    return _copy_table(TABLES / "pc_loadings_named.tsv", "tableS7_pc_loadings.tsv")


def table_s8() -> pd.DataFrame:
    """Table S8: each level of each variable against the rest per component, hallmark."""
    return _copy_table(_one("hallmark_sampletype", "*_metadata_pc_levels.tsv"),
                       "tableS8_pc_levels_HALLMARK.tsv")


def table_s9() -> pd.DataFrame:
    """Table S9: each level of each variable against the rest per component, GO:BP."""
    return _copy_table(_one("gobp_sampletype", "*_metadata_pc_levels.tsv"),
                       "tableS9_pc_levels_GOBP.tsv")


def table_s10() -> pd.DataFrame:
    """Table S10: per-sample quality control."""
    return _copy_table(PROC / "qc_samples.tsv", "tableS10_sample_qc.tsv")


def table_s20() -> pd.DataFrame:
    """Table S20: which component separates which diagnosis, in both collections.

    The lookup behind the marks on Fig. 1b, with the GO:BP answer beside the
    hallmark one so a reader can see which assignments replicate.
    """
    stats = ("component", "pc", "effect", "auc", "fdr", "direction")
    h = diagnosis_components("hallmark_sampletype").drop(columns="mag")
    g = diagnosis_components("gobp_sampletype").drop(columns=["mag", "n"])
    t = (h.rename(columns={c: f"hallmark_{c}" for c in stats})
         .merge(g.rename(columns={c: f"gobp_{c}" for c in stats}),
                on="diagnosis", how="left"))
    t["same_component"] = t["hallmark_component"] == t["gobp_component"]
    t = t.sort_values(["hallmark_pc", "hallmark_effect"],
                      key=lambda s: -s.abs() if s.name == "hallmark_effect"
                      else s)
    t = t.drop(columns=["hallmark_pc", "gobp_pc"])
    STAB.mkdir(parents=True, exist_ok=True)
    t.to_csv(STAB / "tableS20_pc_diagnosis_map.tsv", sep="\t", index=False)
    return t


# ── the rest of the submission bundle ──────────────────────────────────────
def copy_configs() -> list[Path]:
    """Copy every run configuration into ``submission/configs/``."""
    cfg = SUB / "configs"
    cfg.mkdir(parents=True, exist_ok=True)
    out = []
    for c in sorted(CONFIGS.glob("*.yaml")):
        shutil.copy2(c, cfg / c.name)
        out.append(cfg / c.name)
    return out


def _mirror(dest: Path, src: list[str]) -> None:
    """Copy ``src`` into ``dest``, removing only what was copied last time.

    ``figures/raw/`` holds material written by other steps too, so a blanket
    wipe is never safe; each mirror keeps a manifest of its own files.
    """
    dest.mkdir(parents=True, exist_ok=True)
    manifest = dest / "_mirrored.txt"
    if manifest.exists():
        for name in manifest.read_text(encoding="utf-8").split():
            (dest / name).unlink(missing_ok=True)
    for f in src:
        shutil.copy2(f, dest / Path(f).name)
    manifest.write_text("\n".join(sorted(Path(f).name for f in src)) + "\n",
                        encoding="utf-8")


def mirror_pipeline_figures() -> list[Path]:
    """Mirror every plot each pipeline run drew into ``figures/raw/pipeline_<run>/``.

    These let a reader check that arranging a panel into Fig. 1 did not change
    it, which only works if they are the plots the current run emitted.  Runs
    that draw nothing (the GO:BP runs) get no directory.
    """
    written = []
    for run in sorted(p.name for p in PIPE.iterdir() if p.is_dir()):
        src = sorted(glob.glob(str(PIPE / run / "**" / "*.pdf"), recursive=True))
        if not src:
            continue
        _mirror(RAW_FIG / f"pipeline_{run}", src)
        written.append(RAW_FIG / f"pipeline_{run}")
    return written


def run_all() -> None:
    """Draw Fig. 1 and Figs. S5-S7, mirror the pipeline plots, write Tables S1-S10
    and S20 and copy the configurations.  Run :func:`analysis.summaries.run_all`
    first."""
    for draw in (figure1, figure_s5, figure_s6, figure_s7):
        plt.close(draw())
    mirror_pipeline_figures()
    for table in (table_s1, table_s2, table_s3, table_s4, table_s5, table_s6,
                  table_s7, table_s8, table_s9, table_s10, table_s20):
        table()
    copy_configs()
