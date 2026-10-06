"""Fig. S8 and Tables S17-S19 -- MYCN-amplified against non-amplified neuroblastoma.

Brohl et al. report, over the neuroblastoma tumours of this same cohort, that
tumour-infiltrating-lymphocyte enrichment and immune-checkpoint expression are
higher in MYCN-not-amplified (MYCN.NA) than in MYCN-amplified (MYCN.A)
tumours.  Their TIL enrichment is ssGSEA, which is what PathwayTheme computes,
so the two are directly comparable rather than merely consistent.

**The grouping is theirs, not ours.**  MYCN status is already in the delivered
annotation: Brohl Table S1A carries it in `Diagnosis Abbreviation`, joined into
``clinical_annotation.tsv`` as `brohl_diagnosis_code`, and over the 206
neuroblastoma tumours it reads 47 `NB.MYCN.A`, 147 `NB.MYCN.NA` and 12
`NB.Unknown` -- their Table 1, sample for sample.  Nothing is inferred, and the
twelve they could not type are dropped rather than guessed at.

:func:`write_mycn_metadata` writes the label and the metadata files that five
shipped configurations read; the package runs themselves are made from those
configurations, outside this module:

    gobp_nbl_mycn             194 tumours, GO:BP, GO-slim themes on
    hallmark_nbl_mycn         the same 194 over the 50 hallmark processes
    gobp_nbl_mycn_celllines   the 39 neuroblastoma cell lines, same contrast
    gobp_nbl_mycn_4way        all 233 typed samples as a 2 x 2 (MYCN state x
    hallmark_nbl_mycn_4way    tumour / cell line), four contrasts in one run

The cell-line run is a specificity control rather than a replication: a cell
line has no immune infiltrate, so an immune theme that is infiltrate has to
weaken there, while a MYCN programme intrinsic to the tumour cell need not.
The 2 x 2 puts the tumour contrast back among the other three it contains --
amplification within each sample type, and culture within each MYCN state --
from one enrichment pass, one standardisation and one test.

The figure functions read the runs' outputs under ``nbl_mycn/results`` and
``nbl_mycn/four_way/results``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from pathwaytheme.viz import draw_theme_heatmap, draw_top_bars

from .paths import MYCN, PROC, SFIG, STAB

FOUR_WAY_DIR = MYCN / "four_way"

# Brohl et al. Table 1, over the 206 neuroblastoma tumours: the join is checked
# against these counts before the label is used.
BROHL_TUMOURS = {"MYCN.A": 47, "MYCN.NA": 147, "Unknown": 12}

# (metadata file stem, the sample type it keeps)
SUBSETS = [("tumours", "TUMOR"), ("celllines", "CELL_LINE")]

# The four comparisons of the 2 x 2, in the order Fig. S8c shows them.
FOUR_WAY = ["amplified vs not, TUMOURS", "amplified vs not, CELL LINES",
            "cell line vs tumour, AMPLIFIED",
            "cell line vs tumour, NOT AMPLIFIED"]
FOUR_WAY_CONTRAST = {
    "MYCN.A_Tumor_vs_MYCN.NA_Tumor": FOUR_WAY[0],
    "MYCN.A_CellLine_vs_MYCN.NA_CellLine": FOUR_WAY[1],
    "MYCN.A_CellLine_vs_MYCN.A_Tumor": FOUR_WAY[2],
    "MYCN.NA_CellLine_vs_MYCN.NA_Tumor": FOUR_WAY[3],
}
FOUR_WAY_GROUPS = ["MYCN.A_Tumor", "MYCN.NA_Tumor",
                   "MYCN.A_CellLine", "MYCN.NA_CellLine"]
IMMUNE_THEMES = ["immune system process", "inflammatory response",
                 "defense response to other organism"]

# how many of the strongest effects each way the two bar panels show
THEMES_PER_SIDE, HALLMARK_PER_SIDE = 10, 7
_WORDS = {7: "seven", 10: "ten"}


# ── the label and the metadata the five runs read ─────────────────────────
def mycn_labels() -> pd.DataFrame:
    """Brohl's MYCN call per neuroblastoma sample, as it arrives through the join."""
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)
    ann = pd.read_csv(PROC / "clinical_annotation.tsv", sep="\t", index_col=0)
    nbl = meta.index[meta["diagnosis"] == "NBL"]
    code = ann.loc[nbl, "brohl_diagnosis_code"]
    if not code.notna().all():
        raise RuntimeError(f"{int(code.isna().sum())} neuroblastoma samples "
                           f"carry no brohl_diagnosis_code")
    out = pd.DataFrame({
        "sample_type": meta.loc[nbl, "sample_type"],
        "mycn_status": code.str.replace("NB.", "", regex=False),
    })
    got = (out[out["sample_type"] == "TUMOR"]["mycn_status"]
           .value_counts().to_dict())
    if got != BROHL_TUMOURS:
        raise RuntimeError(f"the joined MYCN labels are {got}, Brohl Table 1 "
                           f"reports {BROHL_TUMOURS}")
    return out


def write_mycn_metadata() -> pd.DataFrame:
    """The MYCN label and the metadata files the five MYCN configurations read.

    Writes ``nbl_mycn/mycn_status.tsv``, ``nbl_mycn/metadata_tumours.tsv``,
    ``nbl_mycn/metadata_celllines.tsv`` and
    ``nbl_mycn/four_way/metadata_4way.tsv``, and returns the label.

    The matrix adapter needs metadata for every column of the matrix, so each
    file keeps all 901 rows and marks anything outside its run's subset
    `not_applicable` in the grouping column -- the same mechanism
    gobp_diagnosis.yaml uses to run on tumours only.  The cohort's own
    metadata.tsv is left untouched.
    """
    call = mycn_labels()
    MYCN.mkdir(parents=True, exist_ok=True)
    call.to_csv(MYCN / "mycn_status.tsv", sep="\t")

    for stem, subset in SUBSETS:
        meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)
        keep = call.index[(call["sample_type"] == subset)
                          & (call["mycn_status"] != "Unknown")]
        meta["mycn_status"] = pd.Series("not_applicable", index=meta.index)
        meta.loc[keep, "mycn_status"] = call.loc[keep, "mycn_status"]
        meta.to_csv(MYCN / f"metadata_{stem}.tsv", sep="\t")

    # the 2 x 2: one column holding the four groups
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)
    typed = call[call["mycn_status"] != "Unknown"]
    role = typed["sample_type"].map({"TUMOR": "Tumor", "CELL_LINE": "CellLine"})
    group = typed["mycn_status"] + "_" + role
    if sorted(group.unique()) != sorted(FOUR_WAY_GROUPS):
        raise RuntimeError(f"expected the four groups {FOUR_WAY_GROUPS}, "
                           f"built {sorted(group.unique())}")
    meta["mycn_group"] = group.reindex(meta.index).fillna("not_applicable")
    FOUR_WAY_DIR.mkdir(parents=True, exist_ok=True)
    meta.to_csv(FOUR_WAY_DIR / "metadata_4way.tsv", sep="\t")
    return call


# ── reading the runs ──────────────────────────────────────────────────────
def _read(path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; run the MYCN "
                                f"configurations first")
    return pd.read_csv(path, sep="\t")


def themes(run: str) -> pd.DataFrame:
    """The GO-slim theme summary one of the two-group MYCN runs wrote."""
    return _read(MYCN / "results" / run / "GOBP" / "GOBP_category_summary.tsv")


def hallmark() -> pd.DataFrame:
    """The hallmark moderated-t table of the tumour contrast."""
    return _read(MYCN / "results" / "hallmark_nbl_mycn" / "HALLMARK"
                 / "HALLMARK_differential_moderated_t.tsv")


def four_way_themes() -> pd.DataFrame:
    """The 2 x 2 run's theme summary, relabelled for Fig. S8c.

    The four contrast names are mapped to the readable labels the figure
    prints; an unmapped name would mean the design changed underneath the
    panel, so it is an error rather than a passthrough.
    """
    d = _read(FOUR_WAY_DIR / "results" / "gobp_nbl_mycn_4way" / "GOBP"
              / "GOBP_category_summary.tsv")
    unknown = set(d["comparison"]) - set(FOUR_WAY_CONTRAST)
    if unknown:
        raise RuntimeError(f"the 2 x 2 run holds contrasts Fig. S8c does not "
                           f"name: {sorted(unknown)}")
    d["comparison"] = d["comparison"].map(FOUR_WAY_CONTRAST)
    return d


def theme_table() -> pd.DataFrame:
    """Tumour and cell-line effects of every theme with a significant term."""
    def sig(run, name):
        c = themes(run)
        c = c[(c["significance"] == "significant") & (c["category"] != "unmapped")]
        return c.set_index("category")["mean_effect"].rename(name)

    t = sig("gobp_nbl_mycn", "tumours")
    c = sig("gobp_nbl_mycn_celllines", "cell_lines")
    out = pd.concat([t, c], axis=1)
    out["rank_tumours"] = out["tumours"].rank()
    out["rank_cell_lines"] = out["cell_lines"].rank()
    out["immune"] = out.index.isin(IMMUNE_THEMES)
    return out.sort_values("tumours")


def _terms_per_theme() -> dict:
    return themes("gobp_nbl_mycn").set_index("category")["n_pathways"].to_dict()


# ── the paper items ───────────────────────────────────────────────────────
def figure_s8() -> plt.Figure:
    """Fig. S8: MYCN-amplified against non-amplified neuroblastoma.

    (a) GO-slim themes and (b) the 50 hallmark processes for the tumour
    contrast -- hallmark is an independent collection with its own control:
    the grouping is a MYCN call, so the MYC target sets have to come out up.
    (c) puts that contrast back among the other three the 2 x 2 contains, so
    whether a theme moves with amplification, with culture, or with both is a
    column comparison rather than an assertion.  All three panels are the
    package's own primitives; the theme summary is renamed into the two
    columns the bar panel reads rather than redrawn.
    """
    th, hm, four = theme_table(), hallmark(), four_way_themes()
    counts = _terms_per_theme()
    bars = (th.dropna(subset=["tumours"]).reset_index()
            .rename(columns={"category": "pathway", "tumours": "effect"}))

    n_a, n_na = int(hm["n_case"].iloc[0]), int(hm["n_reference"].iloc[0])
    n_hallmark = int(hm["pathway"].nunique())
    n_themes = int(themes("gobp_nbl_mycn").query("category != 'unmapped'")
                   ["category"].nunique())
    n_four = int(four[(four["significance"] == "significant")
                      & (four["category"] != "unmapped")]["category"].nunique())
    w_th, w_hm = _WORDS[THEMES_PER_SIDE], _WORDS[HALLMARK_PER_SIDE]

    # the two bar panels are stacked on the left so the heat map gets a column
    # of its own and can show every theme that has a significant term
    fig = plt.figure(figsize=(13.6, 10.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.06, 1.0],
                          height_ratios=[1.08, 1.0],
                          wspace=0.46, hspace=0.40,
                          left=0.135, right=0.975, top=0.905, bottom=0.055)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]),
            fig.add_subplot(gs[:, 1])]

    ax = axes[0]
    draw_top_bars(ax, bars, per_side=THEMES_PER_SIDE, label_chars=34,
                  fontsize=8.0)
    ax.set_xlabel("mean effect  (positive = higher in the amplified tumours)",
                  fontsize=8)
    # the count sits outside the end of its own bar, so the axis has to make
    # room on both sides or the longest negative bars write over the theme names
    span = max(abs(v) for v in ax.get_xlim())
    ax.set_xlim(-span * 1.22, span * 1.22)
    pad = span * 0.022
    effect = bars.set_index("pathway")["effect"]
    for i, lbl in enumerate(ax.get_yticklabels()):
        name = lbl.get_text().replace("\n", " ")
        n = counts.get(name)
        if n is None:
            continue
        v = float(effect.get(name, 0.0))
        ax.text(v + (pad if v >= 0 else -pad), i, str(int(n)),
                va="center", ha="left" if v >= 0 else "right", fontsize=5.8,
                color="0.35")
    # every label says which sample type this is: the same contrast exists
    # among the cell lines and the two are not interchangeable
    ax.set_title(f"a  GO-slim themes, GO:BP.  {n_a} amplified against {n_na}\n"
                 f"non-amplified TUMOURS.  {w_th.capitalize()} strongest each "
                 f"way of {n_themes};\nthe number is the terms in the theme",
                 fontsize=8.5, loc="left")

    ax = axes[1]
    draw_top_bars(ax, hm, per_side=HALLMARK_PER_SIDE, short_labels=True,
                  label_chars=30, fontsize=8.0)
    ax.set_xlabel("effect  (positive = higher in the amplified tumours)",
                  fontsize=8)
    ax.set_title(f"b  The {n_hallmark} hallmark processes, the same "
                 f"{n_a + n_na} TUMOURS.\n"
                 "MYC targets are the control: the grouping is a MYCN call",
                 fontsize=8.5, loc="left")

    ax = axes[2]
    draw_theme_heatmap(ax, four, top_n=None, comparisons=FOUR_WAY,
                       colorbar_horizontal=True, row_order="signed",
                       fontsize=5.6, label_chars=60,
                       title="c  the same summarisation over all four "
                             "comparisons the 2 × 2 contains.\n"
                             f"All {n_four} themes; "
                             "column 1 is panel a, the number is the terms")
    # the full labels are wider than a column, so they are numbered and cut to
    # the two facts that separate them: what varies, and within what
    ax.set_xticklabels(["1  amplified\nvs not\nTUMOURS",
                        "2  amplified\nvs not\nCELL LINES",
                        "3  cell line\nvs tumour\nAMPLIFIED",
                        "4  cell line\nvs tumour\nNOT AMPL."], fontsize=5.8)

    fig.suptitle("Fig. S8  MYCN in neuroblastoma.  (a, b) amplified against "
                 f"non-amplified TUMOURS — {n_a} against {n_na}, cell lines "
                 f"excluded — the {w_th} and\n{w_hm} strongest effects each "
                 "way: immune programmes down, proliferative programmes up.  "
                 "(c) the same themes across the whole 2 × 2.",
                 fontsize=9.5, y=0.985)
    SFIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(SFIG / f"figS8_nbl_mycn.{ext}", bbox_inches="tight",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)
    return fig


def table_s17() -> pd.DataFrame:
    """Table S17: GO-slim themes of the MYCN contrast, tumours beside cell lines."""
    out = theme_table()
    out["n_terms_tumours"] = out.index.map(_terms_per_theme())
    out.index.name = "theme"
    STAB.mkdir(parents=True, exist_ok=True)
    out.to_csv(STAB / "tableS17_nbl_mycn_themes.tsv", sep="\t")
    return out


def table_s18() -> pd.DataFrame:
    """Table S18: the 50 hallmark processes, MYCN.A against MYCN.NA tumours."""
    hm = hallmark()
    STAB.mkdir(parents=True, exist_ok=True)
    hm.to_csv(STAB / "tableS18_nbl_mycn_hallmark.tsv", sep="\t", index=False)
    return hm


def table_s19() -> pd.DataFrame:
    """Table S19: theme effects across the four comparisons of the MYCN 2 x 2."""
    four = four_way_themes()
    wide = (four[(four["significance"] == "significant")
                 & (four["category"] != "unmapped")]
            .pivot_table(index="category", columns="comparison",
                         values="mean_effect").reindex(columns=FOUR_WAY))
    wide.index.name = "theme"
    STAB.mkdir(parents=True, exist_ok=True)
    wide.to_csv(STAB / "tableS19_nbl_mycn_four_way.tsv", sep="\t")
    return wide


def run_all() -> None:
    """The MYCN metadata, then Fig. S8 and Tables S17-S19.

    The figure and tables read the outputs of the five MYCN runs, which are
    made from the shipped configurations after :func:`write_mycn_metadata` and
    before the rest; this function assumes they exist.
    """
    write_mycn_metadata()
    figure_s8()
    table_s17()
    table_s18()
    table_s19()
