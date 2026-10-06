"""The cohort annotation as an object of analysis: Figs. S4 and S9-S11, Tables S11-S14.

Four questions about the annotation itself, answered from the delivered run
sheet, the external clinical annotation and the pipeline's own PCA output:

    figure_s9   / table_s11    every annotated field, one strip per field over
                               the 901 samples (Supplementary Methods 8.1)
    figure_s4   / table_s12    per diagnosis, is there one GO:BP component that
                               separates it from the rest?  (Methods 6)
    figure_s10  / table_s13a   every annotated variable against every other,
                               with the support behind each cell (Methods 8.2)
    figure_s11  / table_s13b   where every biological level large enough to
                               contrast sits in every technical variable,
                               against the subset its own field is recorded in
                               (Methods 8.3)
    table_s14                  sample type on PC1/PC4 within each chemistry
                               stratum (Methods 8.4)

Nothing is recomputed that the pipeline already computed: the one-vs-rest
AUCs of Fig. S4 are read from the level-versus-rest grids of Tables S8-S9.

Demographic fields (sex, race, ethnicity, age / disease at diagnosis), the
normals' tissue of origin and the alignment metrics come from
``data/processed/clinical_annotation.tsv``, written by the build step.
"""
from __future__ import annotations

import textwrap
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from scipy import stats

from pathwaytheme.diff.stats import benjamini_hochberg

from .paths import PROC, SFIG, STAB, SUB, pca_dir
from .style import GRID, INK, INK_2, MUTED, styled

CLINICAL = PROC / "clinical_annotation.tsv"

# Fig. S4 reads the whole-cohort level-vs-rest grids (Tables S8-S9), not the
# tumours-only gobp_diagnosis run, whose PCA is a different fit and numbers its
# components differently.
LEVELS_GOBP = STAB / "tableS9_pc_levels_GOBP.tsv"
LEVELS_HALLMARK = STAB / "tableS8_pc_levels_HALLMARK.tsv"
ALL_PCA = pca_dir("gobp_sampletype", "GOBP")
HALL_PCA = pca_dir("hallmark_sampletype", "HALLMARK")

# An AUC this far from 0.5 is what we are willing to call "separated by a single
# component".  Stated on the figure; nothing downstream depends on it.
AUC_STRONG = 0.90
AUC_PARTIAL = 0.75

# Smallest diagnosis summarised in Fig. S4: the contrast size floor.
MIN_CONTRAST_N = 10

# A biological level counts as confounded with a technical variable when this
# much of it sits in one stratum of that variable -- i.e. you cannot contrast
# the level against the rest without also contrasting the chemistry.
CONFOUND_SHARE = 0.90

# Fewest pairwise-complete observations we are willing to report an association
# from.  Several fields are recorded on a small, non-random slice of the cohort
# (tissue of origin on the 117 annotated normals, disease at diagnosis on 85
# osteosarcomas), so some pairs have almost no overlap.  An effect size computed
# on a handful of samples is not a weak finding, it is not a finding, and
# printing one next to an effect computed on 900 samples invites them to be
# read as comparable.
MIN_PAIR_N = 30

# Smallest biological level worth asking the confounding question of.  Below
# this the "share in the largest stratum" statistic is dominated by how few
# samples there are: six samples land on one flowcell often enough by chance
# that the observation would mean nothing.
MIN_LEVEL_N = 8

# A standardised mean difference this large is the numeric counterpart of the
# CONFOUND_SHARE rule: at |d| >= 0.8 the level and the rest of its own annotated
# subset barely overlap on that technical variable.
BIG_D = 0.8

# Categorical technical variables with more strata than this cannot be drawn as
# a composition -- 92 flowcells would be 92 unreadable segments -- so they are
# summarised by how concentrated the level is instead.
MAX_STRATA_TO_STACK = 8

# ---------------------------------------------------------------- palettes
# sample_type gets the manuscript's validated categorical palette.
TYPE_COLOURS = {"TUMOR": "#2a78d6", "CELL_LINE": "#eb6834",
                "NORMAL": "#1baf7a"}

# Diagnosis needs 16 slots, and tissue of origin 20, which is past the point
# where hue alone can carry identity.  Samples are therefore ORDERED by
# diagnosis so every diagnosis is a contiguous block that is direct-labelled;
# colour is the secondary channel and the block label is the primary one.  Hues
# are Okabe-Ito (CVD-safe) first, so the largest levels -- the ones a reader
# actually tracks -- are the ones with the best separation.
QUAL = ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#d55e00", "#56b4e9",
        "#8b5a00", "#332288", "#44aa99", "#999933", "#882255", "#6699cc",
        "#aa4499", "#117733", "#ddcc77", "#771122", "#4a3aa7", "#7a4b2a",
        "#2a78d6", "#eb6834", "#1baf7a", "#b35806"]
MISSING = "#e9e8e4"

HEAT_LO, HEAT_HI = "#f7f7f5", "#1c4f8f"

# ---------------------------------------------------------------- fields
LABELS = {
    "sample_type": "Sample type",
    "diagnosis": "Diagnosis (OncoTree)", "sex": "Sex", "race": "Race",
    "ethnicity": "Ethnicity", "disease_at_diagnosis": "Disease at diagnosis",
    "age_at_diagnosis_years": "Age at diagnosis (y)",
    "vital_status": "Vital status",
    "library_type": "Library chemistry",
    "seq_kit": "Sequencing kit", "flowcell": "Flowcell", "rin": "RIN",
    "zero_fraction": "Zero fraction",
    "tissue_of_origin": "Tissue of origin (normals)",
    "qc_mapping_rate": "Mapping rate", "qc_unique_rate": "Unique rate",
    "qc_duplication_rate": "Duplication rate", "qc_exonic_rate": "Exonic rate",
    "qc_intronic_rate": "Intronic rate", "qc_intragenic_rate": "Intragenic rate",
    "qc_intergenic_rate": "Intergenic rate",
    "qc_fragment_length_mean": "Fragment length", "qc_read_length": "Read length",
    "qc_library_size": "Est. library size", "qc_mapped": "Mapped reads",
    "qc_mapped_unique": "Mapped unique",
    "qc_base_mismatch_rate": "Base mismatch rate",
    "qc_mean_per_base_cov": "Mean per-base cov.", "qc_rrna": "rRNA fraction",
}
NUMERIC = {"rin", "zero_fraction", "age_at_diagnosis_years",
           "qc_mapping_rate", "qc_unique_rate", "qc_duplication_rate",
           "qc_exonic_rate", "qc_intronic_rate", "qc_intragenic_rate",
           "qc_intergenic_rate", "qc_fragment_length_mean", "qc_read_length",
           "qc_library_size", "qc_mapped", "qc_mapped_unique",
           "qc_base_mismatch_rate", "qc_mean_per_base_cov", "qc_rrna"}

# Where each field comes from, and -- for the fields that are still not
# complete after every source has been used -- what the hole is.  This is
# printed beside Fig. S9, because "12% annotated" and "12% annotated, because
# the field only exists for normal tissue" are different statements and only
# the second one is useful.  Every incomplete field in this cohort is
# incomplete structurally: no field is missing at random.
SOURCE = {
    "sample_type": "run sheet; blank = normal",
    "diagnosis": "run sheet; N/A for normals",
    "sex": "Y-linked expression",
    "race": "GDC TARGET: NBL + OS only",
    "ethnicity": "GDC TARGET: NBL + OS only",
    "age_at_diagnosis_years": "GDC TARGET: NBL + OS only",
    "vital_status": "GDC TARGET: NBL + OS only",
    "disease_at_diagnosis": "Brohl S1D: osteosarcoma only",
    "tissue_of_origin": "Brohl S1A: normals only",
    "library_type": "run sheet + Brohl S1A",
    "seq_kit": "run sheet + flowcell",
    "flowcell": "run sheet + sample id",
    "rin": "run sheet + Brohl S1B",
    "zero_fraction": "derived from the matrix",
}
QC_SOURCE = "Brohl S1B"

# One recorded level is longer than any legend column can hold.  Truncating it
# to "Native Hawaiian or Other ..." loses the half that identifies it, so it is
# abbreviated here instead; the table keeps the value as recorded.
LEVEL_ALIAS = {
    "Native Hawaiian or Other Pacific Islander":
        "Native Hawaiian / Pacific Isl.",
}
LEVEL_MAX_CHARS = 30

# Recorded, and deliberately kept out of the association matrix of Fig. S10 --
# with the reason, because "we dropped a column" is not a method.  Each is a
# re-spelling of a column already in the matrix, so including it would report
# the same signal twice and put a guaranteed 1.00 in the grid.
REDUNDANT = {
    "batch": "a delivered column whose values are literally "
             "'<Sample Type>@<Sequencing Centre>'",
    "sequencing_centre": "constant -- every analysed sample is from the public "
                         "RNAseq Landscape centre",
}

# The processed annotation shipped with the submission, so that submission/ is
# a standalone bundle and cannot go stale against data/processed.
SHIPPED_DATA = ("metadata.tsv", "clinical_annotation.tsv", "provenance.tsv",
                "qc_samples.tsv", "cohort_summary.tsv",
                "supplement_link_report.tsv")


# ---------------------------------------------------------------- helpers
def _seq_ramp(v: np.ndarray) -> np.ndarray:
    """Light -> dark single-hue ramp; NaN renders as the missing grey."""
    lo, hi = np.array(to_rgb(HEAT_LO)), np.array(to_rgb(HEAT_HI))
    out = np.tile(np.array(to_rgb(MISSING)), (len(v), 1))
    ok = np.isfinite(v)
    if ok.sum() > 1:
        x = v[ok]
        t = (x - np.nanmin(x)) / max(np.nanmax(x) - np.nanmin(x), 1e-12)
        out[ok] = lo + t[:, None] * (hi - lo)
    return out


def _legend_name(level: str) -> str:
    name = LEVEL_ALIAS.get(level, level)
    return name if len(name) <= LEVEL_MAX_CHARS else name[:LEVEL_MAX_CHARS - 1] + "…"


def _cat_colours(levels: list[str], palette: list[str] | dict) -> dict:
    if isinstance(palette, dict):
        return palette
    return {lv: palette[i % len(palette)] for i, lv in enumerate(levels)}


def _save(fig: plt.Figure, stem: str) -> None:
    SFIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(SFIG / f"{stem}.pdf")
    fig.savefig(SFIG / f"{stem}.png")


def _write(df: pd.DataFrame, name: str, index: bool = True) -> None:
    STAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(STAB / name, sep="\t", index=index)


# ---------------------------------------------------------------- statistics
#
# Three measures, one scale.  Every pair of annotated variables falls into
# exactly one of three cases, and each case gets the measure appropriate to it.
# All three are bounded on [0, 1], 0 = independent, 1 = one variable
# determines the other, so they can share one colour bar; which one produced a
# given cell is recorded per row in Table S13a.
#
#   categorical x categorical   bias-corrected Cramer's V
#   categorical x continuous    bias-corrected eta  (= sqrt of adjusted eta^2)
#   continuous  x continuous    |Spearman rho|
#
# The bias corrections are the point.  Both Cramer's V and eta^2 grow with the
# number of levels at fixed n, so a 20-level field observed on 117 samples
# (tissue of origin) will out-score a 3-level field observed on 901 unless the
# correction is applied to both.
#
# Each measure returns (effect, p, n, reason).  `reason` is "ok" or names why
# nothing is reportable, because "the cell is blank" and "the cell is blank
# because two variables never co-occur" are different facts about the cohort and
# the figure should not render them identically:
#
#   n_below_floor          fewer than MIN_PAIR_N samples carry both variables
#   one_variable_constant  one has a single level wherever the other is recorded
#                          -- e.g. every sample with a recorded race is polyA, so
#                          race and chemistry cannot be told apart at all
#   levels_saturate_n      the bias correction has no degrees of freedom left,
#                          the categorical variable having about as many levels
#                          as there are observations.  A raw Cramer's V would
#                          return ~1 here and mean nothing by it

def _cramers_v(a: pd.Series, b: pd.Series) -> tuple[float, float, int, str]:
    """Bias-corrected Cramer's V between two categorical variables."""
    ok = a.notna() & b.notna()
    n_ok = int(ok.sum())
    if n_ok < MIN_PAIR_N:
        return np.nan, np.nan, n_ok, "n_below_floor"
    tab = pd.crosstab(a[ok], b[ok])
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return np.nan, np.nan, n_ok, "one_variable_constant"
    chi2, p, _, _ = stats.chi2_contingency(tab)
    n = int(tab.to_numpy().sum())
    phi2 = chi2 / n
    r, k = tab.shape
    phi2c = max(0.0, phi2 - (k - 1) * (r - 1) / (n - 1))
    rc = r - (r - 1) ** 2 / (n - 1)
    kc = k - (k - 1) ** 2 / (n - 1)
    denom = min(kc - 1, rc - 1)
    if denom <= 0:
        return np.nan, np.nan, n, "levels_saturate_n"
    return float(np.sqrt(phi2c / denom)), float(p), n, "ok"


def _eta(cat: pd.Series, num: pd.Series) -> tuple[float, float, int, str]:
    """Bias-corrected eta between a categorical and a continuous variable.

    The adjusted eta^2 subtracts the variance a grouping of this many levels
    would explain by chance,

        eta^2_adj = (SS_between - (k - 1) * MS_within) / SS_total,

    which is the omega^2 numerator; it is clipped at 0 and square-rooted so the
    result sits on the same 0-1 scale as Cramer's V.  Without the correction a
    20-level grouping of 117 samples reports a large eta whatever the data.
    """
    ok = cat.notna() & num.notna()
    n_ok = int(ok.sum())
    if n_ok < MIN_PAIR_N:
        return np.nan, np.nan, n_ok, "n_below_floor"
    groups = [g.to_numpy() for _, g in num[ok].groupby(cat[ok]) if len(g) >= 2]
    if len(groups) < 2:
        return np.nan, np.nan, n_ok, "one_variable_constant"
    y = np.concatenate(groups)
    n, k = len(y), len(groups)
    if n - k < 1:
        return np.nan, np.nan, n, "levels_saturate_n"
    grand = y.mean()
    ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ss_t = ((y - grand) ** 2).sum()
    if ss_t <= 0:
        return np.nan, np.nan, n, "one_variable_constant"
    ms_w = (ss_t - ss_b) / (n - k)
    eta2 = max(0.0, (ss_b - (k - 1) * ms_w) / ss_t)
    try:
        _, p = stats.kruskal(*groups)
    except ValueError:
        p = np.nan
    return float(np.sqrt(eta2)), float(p), n, "ok"


def _abs_spearman(a: pd.Series, b: pd.Series) -> tuple[float, float, int, str]:
    """|Spearman rho| between two continuous variables.

    A continuous variable must not be fed to a crosstab or a group-wise eta:
    that treats every distinct value as its own level, and grouping a few
    hundred samples by their nearly-unique ages manufactures a large effect
    out of nothing at all.
    """
    ok = a.notna() & b.notna()
    n_ok = int(ok.sum())
    if n_ok < MIN_PAIR_N:
        return np.nan, np.nan, n_ok, "n_below_floor"
    x, y = a[ok].to_numpy(float), b[ok].to_numpy(float)
    if np.nanstd(x) <= 0 or np.nanstd(y) <= 0:
        return np.nan, np.nan, n_ok, "one_variable_constant"
    rho, p = stats.spearmanr(x, y)
    if not np.isfinite(rho):
        return np.nan, np.nan, n_ok, "one_variable_constant"
    return float(abs(rho)), float(p), n_ok, "ok"


def _associate(a: pd.Series, b: pd.Series, a_num: bool,
               b_num: bool) -> tuple[float, float, int, str, str]:
    """Dispatch a variable pair to the measure its two types call for."""
    if a_num and b_num:
        v, p, n, why = _abs_spearman(a, b)
        return v, p, n, "abs_spearman_rho", why
    if a_num and not b_num:
        v, p, n, why = _eta(b, a)
        return v, p, n, "eta_adjusted", why
    if b_num and not a_num:
        v, p, n, why = _eta(a, b)
        return v, p, n, "eta_adjusted", why
    v, p, n, why = _cramers_v(a, b)
    return v, p, n, "cramers_v_corrected", why


def _eta_squared_raw(cat: pd.Series, num: pd.Series) -> tuple[float, float, int]:
    """Raw eta^2 for a categorical against a numeric, with Kruskal's p.

    Deliberately the *uncorrected* eta^2, unlike `_eta`: this quantity is
    compared against the pipeline's own attribution grid, which reports raw
    eta^2 (Methods 6), and the conditional test is a comparison of the same
    statistic between strata rather than a cross-variable ranking.
    """
    ok = cat.notna() & num.notna()
    if ok.sum() < 10:
        return np.nan, np.nan, int(ok.sum())
    groups = [g.to_numpy() for _, g in num[ok].groupby(cat[ok]) if len(g) >= 3]
    if len(groups) < 2:
        return np.nan, np.nan, int(ok.sum())
    y = np.concatenate(groups)
    grand = y.mean()
    ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ss_t = ((y - grand) ** 2).sum()
    try:
        _, p = stats.kruskal(*groups)
    except ValueError:
        p = np.nan
    return (float(ss_b / ss_t) if ss_t > 0 else np.nan, float(p),
            int(sum(len(g) for g in groups)))


# ---------------------------------------------------------------- load
def load_annotation(verbose: bool = False
                    ) -> tuple[pd.DataFrame, list[str], list[str]]:
    """The analysed annotation: the run-sheet metadata joined to the clinical one.

    Returns the per-sample table and the names of its biological and technical
    fields, each restricted to fields with more than one observed level.
    """
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", dtype={"sample_id": str})
    meta = meta.set_index("sample_id")

    for c in ("rin", "zero_fraction"):
        meta[c] = pd.to_numeric(meta[c], errors="coerce")

    # `flowcell` is built with the matrix, not read from the run sheet: the run
    # sheet's own column holds SRA run accessions for the 85 TARGET
    # osteosarcomas, one per sample, and holds nothing at all for the normals.
    biological = ["sample_type", "diagnosis"]
    technical = ["library_type", "seq_kit", "sequencing_centre", "flowcell",
                 "rin", "zero_fraction"]

    if CLINICAL.exists():
        cl = pd.read_csv(CLINICAL, sep="\t").set_index("sample_id")
        for c in cl.columns:
            if c not in meta.columns:
                meta[c] = cl[c].reindex(meta.index)
        # tissue of origin is biological, and it is the field that resolves the
        # normals: Brohl S1A records "Normal-cerebellum", "Normal-prostate", ...
        n_before = len(biological)
        for c in ("sex", "tissue_of_origin", "race", "ethnicity",
                  "disease_at_diagnosis", "age_at_diagnosis_years",
                  "vital_status"):
            if c in meta.columns and meta[c].notna().any():
                biological.append(c)
        # Brohl S1B's alignment/library QC metrics are extra *technical*
        # variables, which is what makes "is this component chemistry?" a
        # testable question over more than library type alone
        technical += [c for c in meta.columns
                      if c.startswith("qc_") and meta[c].notna().any()]
        if verbose:
            print(f"   clinical annotation merged: +{len(biological) - n_before} "
                  f"biological, +{sum(c.startswith('qc_') for c in technical)} "
                  f"QC technical")
    else:
        warnings.warn(f"{CLINICAL} absent -- sex / race / ethnicity / disease "
                      "at diagnosis / tissue of origin / alignment metrics are "
                      "left out of the annotation analyses", stacklevel=2)

    # A variable with one observed level cannot be associated with anything,
    # and `sequencing_centre` is constant by construction, so it drops out here.
    def informative(c: str) -> bool:
        return meta[c].notna().any() and (c in NUMERIC or meta[c].nunique() > 1)

    dropped = [c for c in technical + biological if not informative(c)]
    if dropped and verbose:
        print(f"   single-level / empty, not testable: {dropped}")
    technical = [c for c in technical if informative(c)]
    biological = [c for c in biological if informative(c)]
    return meta, biological, technical


# ---------------------------------------------------------------- Fig. S9 / Table S11
def _annotation_matrix(meta: pd.DataFrame, biological: list[str],
                       technical: list[str]) -> pd.DataFrame:
    """Every analysed field, samples ordered by sample type, diagnosis, chemistry."""
    order = meta.sort_values(
        ["sample_type", "diagnosis", "library_type"],
        na_position="last", kind="mergesort").index
    return meta.loc[order, biological + technical]


def table_s11() -> pd.DataFrame:
    """Table S11: the annotation matrix behind Fig. S9 -- every analysed field
    for every sample, as text, in the figure's sample order."""
    meta, biological, technical = load_annotation()
    out = _annotation_matrix(meta, biological, technical).copy()
    out.index.name = "sample_id"
    _write(out, "tableS11_annotation_matrix.tsv")
    return out


@styled
def figure_s9() -> plt.Figure:
    """Fig. S9: the cohort annotation -- one strip per analysed field over the
    901 samples, with each field's source and every categorical level listed.

    Only recorded and joined fields are drawn.  A principal component is an
    output of the analysis, not an annotation of the cohort, so it is not shown
    beside the recorded fields; and `batch` is `<sample type>@<sequencing
    centre>` with one centre, which is `sample_type` drawn a second time.
    """
    meta, biological, technical = load_annotation()
    m = _annotation_matrix(meta, biological, technical)

    blocks = [("Biological", biological), ("Technical", technical)]
    rows = [c for _, cols in blocks for c in cols]

    # per-row RGB, and the legend entries that go with it
    img, legends = [], []
    for c in rows:
        if c in NUMERIC:
            img.append(_seq_ramp(m[c].to_numpy(float)))
            continue
        s = m[c].astype("object")
        counts = s.value_counts()
        levels = list(counts.index)
        pal = TYPE_COLOURS if c == "sample_type" else QUAL
        cmap = _cat_colours(levels, pal)
        img.append(np.array([to_rgb(cmap.get(v, MISSING)) if isinstance(v, str)
                             else to_rgb(MISSING) for v in s]))
        # Every level is listed.  A truncated legend ("+7 more") hides exactly
        # the small levels a reader is checking the figure for, and the space it
        # saves is space this figure has.
        if c != "flowcell":
            legends.append((c, [(_legend_name(lv), cmap.get(lv, MISSING),
                                 int(counts[lv])) for lv in levels]))
    img = np.stack(img)

    n = len(m)
    # Lay the legend out before sizing the figure.  Every level of every field
    # is listed, which is 60-odd entries, so they go in columns -- and a column
    # break must fall *between* fields, never inside one, or a field's levels end
    # up in a different column from its name.  That constraint means the columns
    # cannot be balanced by simply dividing the entry count: the split point is
    # chosen to minimise the tallest column, and the figure is then made tall
    # enough for whichever of the strip stack and that column is taller.
    HEADER, LINE_IN, ROW_IN = 1.35, 0.108, 0.163
    heights = [len(items) + HEADER for _, items in legends]
    cuts = range(1, len(legends))
    cut = min(cuts, key=lambda c: max(sum(heights[:c]), sum(heights[c:])))
    col_of = [0 if i < cut else 1 for i in range(len(legends))]
    tallest = max(sum(heights[:cut]), sum(heights[cut:]))
    fig_h = max(ROW_IN * len(rows), LINE_IN * tallest) + 1.5

    fig = plt.figure(figsize=(8.6, fig_h))
    gs = GridSpec(1, 3, width_ratios=[1, 0.28, 0.62], wspace=0.03, figure=fig,
                  left=0.250, right=0.995, top=0.905, bottom=0.055)
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(img, aspect="auto", interpolation="nearest",
              extent=(0, n, len(rows), 0))
    ax.set_yticks(np.arange(len(rows)) + 0.5)
    # Coverage goes on the row label: without it a 12%-complete row is visually
    # indistinguishable from a broken one, since the grey dominates either way.
    ax.set_yticklabels(
        [f"{LABELS.get(c, c)}  ({int(m[c].notna().sum())})"
         if m[c].notna().sum() < n else LABELS.get(c, c) for c in rows],
        fontsize=6)
    ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    for i in range(1, len(rows)):
        ax.axhline(i, color="white", lw=0.9)
    y = 0
    for name, cols in blocks:
        if not cols:
            continue
        ax.axhline(y, color=INK_2, lw=0.8)
        ax.text(-0.405, y + len(cols) / 2, name,
                transform=ax.get_yaxis_transform(), rotation=90, va="center",
                ha="center", fontsize=6.5, color=INK_2, fontweight="bold")
        y += len(cols)
    ax.axhline(y, color=INK_2, lw=0.8)

    # direct labels for the sample-type blocks, so identity is never colour-alone
    st = m["sample_type"].to_numpy()
    start, k = 0, 0
    for i in range(1, n + 1):
        if i == n or st[i] != st[start]:
            if i - start >= 25:
                ax.text((start + i) / 2, -0.32 - 1.05 * (k % 2),
                        f"{st[start]}  n={i - start}",
                        ha="center", va="bottom", fontsize=5.5, color=INK)
                k += 1
            if i < n:
                ax.axvline(i, color=INK, lw=0.7)
            start = i
    ax.set_title(f"a  Cohort annotation, all {len(rows)} analysed fields "
                 f"({n} samples, ordered by sample type then diagnosis).\n"
                 "Grey = not annotated; (n) on a label = samples annotated",
                 fontsize=8, pad=36)

    # ---- where each field comes from, and why it stops where it stops
    #
    # The run sheet has a row for the 783 tumours and cell lines and no row at
    # all for the 118 normals, and the external clinical tables exist only for
    # the TARGET projects.  Every hole in this figure is one of those two facts
    # -- none of it is missing at random -- and a reader cannot tell that from
    # the grey.  So the origin of each field is printed against it.
    sax = fig.add_subplot(gs[0, 1])
    sax.axis("off")
    sax.set_xlim(0, 1)
    sax.set_ylim(len(rows), 0)
    for i, c in enumerate(rows):
        src = SOURCE.get(c, QC_SOURCE if c.startswith("qc_") else "")
        full = int(m[c].notna().sum()) == n
        sax.text(0.04, i + 0.5, src, va="center", ha="left", fontsize=4.9,
                 color=INK_2 if full else INK)
    sax.set_title("source", fontsize=6, pad=36, color=INK_2, loc="left")

    # ---- legend: every level of every categorical field, in two columns
    lax = fig.add_subplot(gs[0, 2])
    lax.axis("off")
    # one legend line is LINE_IN inches, expressed in this axes' own fraction,
    # so both columns use the same spacing whatever the figure height came out at
    axes_h_in = fig_h * (0.905 - 0.055)
    step = LINE_IN / axes_h_in
    ypos = [1.0, 1.0]
    for (c, items), col in zip(legends, col_of):
        x = col * 0.52
        lax.text(x, ypos[col], LABELS.get(c, c), fontsize=6, fontweight="bold",
                 color=INK, transform=lax.transAxes, va="top")
        ypos[col] -= step * HEADER
        for lv, colr, cnt in items:
            lax.add_patch(Rectangle((x + 0.005, ypos[col] - step * 0.74),
                                    0.038, step * 0.66,
                                    transform=lax.transAxes, facecolor=colr,
                                    edgecolor="white", lw=0.4, clip_on=False))
            lax.text(x + 0.052, ypos[col] - step * 0.41, f"{lv}  ({cnt})",
                     fontsize=5.2, color=INK_2, transform=lax.transAxes,
                     va="center")
            ypos[col] -= step
    lax.set_title("b  Every level of every categorical field, in full",
                  fontsize=8, pad=36)

    _save(fig, "figS9_cohort_metadata")
    incomplete = [c for c in rows if int(m[c].notna().sum()) < n]
    print(f"   figS9: {len(rows)} fields x {n} samples, "
          f"{sum(len(i) for _, i in legends)} legend entries, all shown; "
          f"{len(rows) - len(incomplete)} fields complete, "
          f"{len(incomplete)} not: {incomplete}")
    return fig


# ---------------------------------------------------------------- Fig. S4 / Table S12
def _diagnosis_separation() -> dict:
    """Diagnosis x GO:BP component one-vs-rest AUCs and the per-diagnosis summary.

    Separation is symmetric in sign -- a component that puts a diagnosis at the
    bottom separates it exactly as well as one that puts it at the top -- so the
    summary quantity is max(AUC, 1 - AUC), which runs 0.5 (no separation) to
    1.0 (complete).
    """
    lv = pd.read_csv(LEVELS_GOBP, sep="\t")
    lv = lv[lv["variable"] == "diagnosis"].copy()
    hm = pd.read_csv(LEVELS_HALLMARK, sep="\t")
    hm = hm[hm["variable"] == "diagnosis"].copy()

    pcs = sorted(lv["component"].unique(), key=lambda x: int(x[2:]))
    auc = lv.pivot(index="level", columns="component", values="auc")[pcs]
    fdr = lv.pivot(index="level", columns="component", values="fdr")[pcs]
    nlev = lv.groupby("level")["n_level"].max()

    dev = (auc - 0.5).abs()
    sep = dev + 0.5                      # max(AUC, 1 - AUC), the plotted matrix
    best_pc = dev.idxmax(axis=1)
    best_auc = auc.to_numpy()[np.arange(len(auc)), [pcs.index(p) for p in best_pc]]

    hm_pcs = sorted(hm["component"].unique(), key=lambda x: int(x[2:]))
    hm_auc = hm.pivot(index="level", columns="component", values="auc").reindex(
        index=auc.index, columns=hm_pcs)
    hm_sep = (hm_auc - 0.5).abs().max(axis=1) + 0.5

    summary = pd.DataFrame({
        "n": nlev.reindex(auc.index), "best_pathway_pc": best_pc,
        "best_pathway_auc": best_auc,
        "separation": np.maximum(best_auc, 1 - best_auc),
        "best_pathway_fdr": [fdr.loc[d, best_pc[d]] for d in auc.index],
        "hallmark_separation": hm_sep,
    }).sort_values("separation", ascending=False)
    return {"pcs": pcs, "auc": auc, "sep": sep, "summary": summary}


def table_s12() -> pd.DataFrame:
    """Table S12: diagnosis x GO:BP component one-vs-rest AUC, with the best
    component per diagnosis and the hallmark separation beside it as
    replication."""
    d = _diagnosis_separation()
    summary = d["summary"]
    tbl = d["auc"].copy()
    tbl.columns = [f"gobp_{c}_auc" for c in tbl.columns]
    tbl = tbl.join(summary)
    tbl.index.name = "diagnosis"
    tbl = tbl.sort_values("separation", ascending=False)
    _write(tbl, "tableS12_diagnosis_pc_auc.tsv")
    ok = summary[summary["n"] >= MIN_CONTRAST_N]
    nsep = int((ok["separation"] >= AUC_STRONG).sum())
    part = int(((ok["separation"] >= AUC_PARTIAL)
                & (ok["separation"] < AUC_STRONG)).sum())
    print(f"   tableS12: of {len(ok)} diagnoses at n >= {MIN_CONTRAST_N}, {nsep} "
          f"separated at >= {AUC_STRONG} and {part} partially "
          f"({AUC_PARTIAL}-{AUC_STRONG}) by a single component")
    return tbl


@styled
def figure_s4() -> plt.Figure:
    """Fig. S4: one-vs-rest separation of every diagnosis by every GO:BP
    component (a), and the separating component drawn for each diagnosis that
    clears the threshold (b).

    max(AUC, 1 - AUC) is a one-ended scale, so it gets a one-hue sequential
    ramp: a diverging map on raw AUC would spend two hues on a sign that is
    arbitrary in a PCA, and make a perfectly separated diagnosis look like two
    opposite findings depending on which way the component happened to point.
    The direction is still reported, as an arrow on the cells that clear the
    threshold.
    """
    meta, _, _ = load_annotation()
    d = _diagnosis_separation()
    pcs, auc, sep, summary = d["pcs"], d["auc"], d["sep"], d["summary"]
    sc = pd.read_csv(ALL_PCA / "GOBP___all___pca_scores.tsv", sep="\t", index_col=0)
    ve = pd.read_csv(ALL_PCA / "GOBP___all___pca_variance_explained.tsv", sep="\t")
    ve = ve.set_index("PC")["variance_explained"]

    fig = plt.figure(figsize=(7.2, 4.9))
    gs = GridSpec(2, 1, height_ratios=[1, 0.72], hspace=0.34, figure=fig,
                  left=0.155, right=0.90, top=0.90, bottom=0.09)

    # (a) diagnosis x component separation -- one hue, 0.5 -> 1.0
    ax = fig.add_subplot(gs[0])
    o = summary.index
    S = sep.loc[o]
    im = ax.imshow(S.to_numpy(), aspect="auto", cmap="Blues", vmin=0.5, vmax=1.0,
                   interpolation="nearest")
    ax.set_xticks(np.arange(len(pcs)))
    ax.set_xticklabels([f"{p}\n{ve[p] * 100:.1f}%" for p in pcs], fontsize=5.5)
    ax.set_yticks(np.arange(len(o)))
    ax.set_yticklabels([f"{d}  (n={int(summary.loc[d, 'n'])})"
                        + ("  †" if summary.loc[d, "n"] < MIN_CONTRAST_N else "")
                        for d in o], fontsize=6)
    for i, dx in enumerate(o):
        for j, p in enumerate(pcs):
            s_ij, a_ij = S.loc[dx, p], auc.loc[dx, p]
            if not np.isfinite(s_ij):
                continue
            if s_ij >= AUC_STRONG:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor=INK, lw=0.9))
                arrow = "▲" if a_ij > 0.5 else "▼"
                ax.text(j, i, f"{s_ij:.2f}{arrow}", ha="center", va="center",
                        fontsize=4.4, color="white" if s_ij > 0.80 else INK)
    ax.set_title("a  One-vs-rest separation of each diagnosis by each GO:BP "
                 f"component, max(AUC, 1−AUC); all {len(sc)} samples.\n"
                 f"Boxed ≥ {AUC_STRONG:.2f}; ▲ the diagnosis scores high on that "
                 f"component, ▼ low.  † n < {MIN_CONTRAST_N}, below the contrast "
                 "size floor",
                 fontsize=8, pad=6)
    cb = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.012)
    cb.set_label("separation, max(AUC, 1−AUC)", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    cb.outline.set_visible(False)

    # (b) the separating component, drawn, for every diagnosis that clears it
    ax2 = fig.add_subplot(gs[1])
    top = [dx for dx in o if summary.loc[dx, "separation"] >= AUC_STRONG
           and summary.loc[dx, "n"] >= MIN_CONTRAST_N]
    dg = meta["diagnosis"].reindex(sc.index)
    for i, dx in enumerate(top):
        pc = summary.loc[dx, "best_pathway_pc"]
        v = sc[pc]
        z = (v - v.mean()) / v.std()
        zi, zo = z[dg == dx], z[dg != dx]
        rng = np.random.default_rng(i)
        ax2.scatter(zo, np.full(len(zo), i) + rng.uniform(-0.13, 0.13, len(zo)),
                    s=1.6, color=GRID, lw=0, zorder=1)
        ax2.scatter(zi, np.full(len(zi), i) + rng.uniform(-0.13, 0.13, len(zi)),
                    s=3.4, color=TYPE_COLOURS["TUMOR"], lw=0.25,
                    edgecolor="white", zorder=2)
        side = "high" if summary.loc[dx, "best_pathway_auc"] > 0.5 else "low"
        ax2.text(1.005, i, f"{dx} · {pc} · {summary.loc[dx, 'separation']:.2f} "
                           f"({side})",
                 transform=ax2.get_yaxis_transform(), fontsize=5.4, color=INK,
                 va="center", ha="left")
    ax2.set_yticks([])
    ax2.set_ylim(len(top) - 0.5, -0.5)
    ax2.set_xlim(-5.5, 5.5)
    ax2.set_xlabel("component score (z)", fontsize=6)
    ax2.set_title("b  The separating component drawn, for every diagnosis "
                  f"clearing {AUC_STRONG:.2f} at n ≥ {MIN_CONTRAST_N} "
                  "(coloured = that diagnosis, grey = every other sample)",
                  fontsize=8, pad=4)
    ax2.spines["left"].set_visible(False)

    _save(fig, "figS4_diagnosis_pc")
    print(f"   figS4: {len(top)} diagnoses drawn in panel b")
    return fig


# ---------------------------------------------------------------- Fig. S10 / Table S13a
def _association(meta: pd.DataFrame, biological: list[str],
                 technical: list[str]) -> tuple[pd.DataFrame, list[str], dict]:
    """Every pair of analysed variables: measure, effect, p, support, BH-FDR.

    The grid is the full symmetric one -- biological x biological, biological x
    technical and technical x technical -- because all three blocks matter for
    reading it: whether the biological variables are separable from each other,
    and whether the twenty technical variables are twenty pieces of evidence or
    one (fifteen alignment metrics that agree with each other to |rho| > 0.9
    are not fifteen independent chances for a component to be technical).
    """
    variables = [v for v in biological + technical if v not in REDUNDANT]
    kind = {v: ("biological" if v in biological else "technical")
            for v in variables}
    is_num = {v: v in NUMERIC for v in variables}

    recs = []
    for i, a in enumerate(variables):
        for b in variables[i + 1:]:
            v, p, n, measure, why = _associate(meta[a], meta[b],
                                               is_num[a], is_num[b])
            recs.append({"variable_a": a, "variable_b": b,
                         "block": "x".join(sorted({kind[a], kind[b]})),
                         "measure": measure, "effect": v, "p_value": p,
                         "n_pairwise": n, "status": why,
                         "reportable": bool(np.isfinite(v))})
    assoc = pd.DataFrame(recs)
    assoc["fdr"] = benjamini_hochberg(assoc["p_value"].to_numpy())
    return assoc, variables, kind


def table_s13a() -> pd.DataFrame:
    """Table S13a: every pair of the annotated variables -- block, measure,
    effect, samples carrying both, raw and adjusted p, and why no effect is
    reportable where none is -- plus the recorded columns left out, with the
    reason."""
    meta, biological, technical = load_annotation()
    assoc, variables, _ = _association(meta, biological, technical)
    excl = pd.DataFrame([{"variable_a": v, "variable_b": "(excluded)",
                          "block": "excluded", "measure": "not tested",
                          "effect": np.nan, "p_value": np.nan,
                          "n_pairwise": int(meta[v].notna().sum()),
                          "reportable": False, "fdr": np.nan,
                          "status": "not_tested", "note": why}
                         for v, why in REDUNDANT.items() if v in meta.columns])
    assoc["note"] = ""
    out = pd.concat([assoc, excl], ignore_index=True)
    _write(out, "tableS13a_annotation_association.tsv", index=False)

    rep = assoc[assoc["reportable"]]
    sig = int((rep["fdr"] < 0.05).sum())
    print(f"   tableS13a: {len(variables)} variables, {len(assoc)} pairs, "
          f"{len(rep)} reportable at n >= {MIN_PAIR_N}, {sig} at FDR < 0.05")
    for blk, g in rep.groupby("block"):
        print(f"           {blk}: {len(g)} pairs, median effect "
              f"{g['effect'].median():.3f}")
    return out


@styled
def figure_s10() -> plt.Figure:
    """Fig. S10: association between every pair of the annotated variables in
    one symmetric matrix -- the effect above the diagonal, the samples it rests
    on below it."""
    meta, biological, technical = load_annotation()
    assoc, variables, kind = _association(meta, biological, technical)

    # marks for the cells that carry no effect, one per reason
    BLANK = {"n_below_floor": ("·", MUTED),
             "one_variable_constant": ("=", MUTED),
             "levels_saturate_n": ("~", MUTED)}

    k = len(variables)
    idx = {v: i for i, v in enumerate(variables)}
    M = np.full((k, k), np.nan)
    Q = np.full((k, k), np.nan)
    N = np.full((k, k), np.nan)
    W = np.empty((k, k), dtype=object)
    for r in assoc.itertuples():
        i, j = idx[r.variable_a], idx[r.variable_b]
        M[i, j] = M[j, i] = r.effect
        Q[i, j] = Q[j, i] = r.fdr
        N[i, j] = N[j, i] = r.n_pairwise
        W[i, j] = W[j, i] = r.status

    n_bio = sum(1 for v in variables if kind[v] == "biological")

    # One matrix, not two.  The grid is symmetric, so drawing it whole spends
    # half the canvas twice over; instead the effect goes above the diagonal and
    # the number of samples the effect rests on goes below it, which also puts
    # every cell next to its own support rather than on a separate panel the
    # reader has to hold in their head.
    LEFT, RIGHT, TOP, BOTTOM = 0.215, 0.855, 0.930, 0.175
    panel_in = 7.4 * (RIGHT - LEFT)
    fig_h = panel_in / (TOP - BOTTOM)
    fig = plt.figure(figsize=(7.4, fig_h))
    ax = fig.add_axes((LEFT, BOTTOM, RIGHT - LEFT, TOP - BOTTOM))

    upper = np.where(np.triu(np.ones((k, k)), 1) > 0, M, np.nan)
    lower = np.where(np.tril(np.ones((k, k)), -1) > 0, N, np.nan)
    im_e = ax.imshow(upper, aspect="auto", cmap="Blues", vmin=0.0, vmax=1.0,
                     interpolation="nearest")
    im_n = ax.imshow(lower, aspect="auto", cmap="Greys", vmin=0.0,
                     vmax=float(len(meta)), interpolation="nearest")

    ax.set_xticks(np.arange(k))
    ax.set_xticklabels([LABELS.get(v, v) for v in variables], rotation=90,
                       fontsize=5.4)
    ax.set_yticks(np.arange(k))
    ax.set_yticklabels([LABELS.get(v, v) for v in variables], fontsize=5.4)
    for i, (xl, yl) in enumerate(zip(ax.get_xticklabels(), ax.get_yticklabels())):
        colour = (TYPE_COLOURS["CELL_LINE"] if kind[variables[i]] == "technical"
                  else INK)
        xl.set_color(colour)
        yl.set_color(colour)

    for i in range(k):
        for j in range(k):
            if i == j:
                ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1,
                                       facecolor=GRID, edgecolor="white",
                                       lw=0.4, zorder=3))
                continue
            if j > i:                                   # effect, above
                v = M[i, j]
                if not np.isfinite(v):
                    mark, colour = BLANK.get(W[i, j], ("", MUTED))
                    if mark:
                        ax.text(j, i, mark, ha="center", va="center",
                                fontsize=6, color=colour)
                    continue
                star = "*" if np.isfinite(Q[i, j]) and Q[i, j] < 0.05 else ""
                ax.text(j, i, f"{v:.2f}{star}".lstrip("0"), ha="center",
                        va="center", fontsize=3.85,
                        color="white" if v > 0.55 else INK)
            else:                                       # support, below
                v = N[i, j]
                if not np.isfinite(v):
                    continue
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        fontsize=3.85, color="white" if v > 520 else INK_2)

    # block rules: biological | technical, on both axes
    ax.axhline(n_bio - 0.5, color=INK, lw=1.0, zorder=4)
    ax.axvline(n_bio - 0.5, color=INK, lw=1.0, zorder=4)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # Which half of the matrix is which block, said with a bracket on the axis
    # rather than with a sentence about the colour of the tick labels.
    ORANGE = TYPE_COLOURS["CELL_LINE"]
    for span, name, colour in (((0, n_bio - 1), "biological", INK),
                               ((n_bio, k - 1), "technical", ORANGE)):
        mid = (span[0] + span[1]) / 2
        ax.text(-0.235, mid, name, transform=ax.get_yaxis_transform(),
                rotation=90, ha="center", va="center", fontsize=7,
                fontweight="bold", color=colour)
        ax.text(mid, -0.235, name, transform=ax.get_xaxis_transform(),
                ha="center", va="center", fontsize=7, fontweight="bold",
                color=colour)

    cax_e = fig.add_axes((RIGHT + 0.020, TOP - 0.235, 0.013, 0.235))
    cb = fig.colorbar(im_e, cax=cax_e)
    cb.set_label("association\n(0 = independent, 1 = determines)", fontsize=5.6)
    cb.ax.tick_params(labelsize=5)
    cb.outline.set_visible(False)
    cax_n = fig.add_axes((RIGHT + 0.020, TOP - 0.520, 0.013, 0.235))
    cb2 = fig.colorbar(im_n, cax=cax_n)
    cb2.set_label(f"samples carrying both (of {len(meta)})", fontsize=5.6)
    cb2.ax.tick_params(labelsize=5)
    cb2.outline.set_visible(False)

    # a key, not a paragraph: every mark on the matrix, one entry each.  The
    # step follows the wrapped height, so a multi-line entry never sits on the
    # next one.
    key = [("*", "FDR < 0.05"),
           ("·", f"fewer than {MIN_PAIR_N} samples carry both"),
           ("=", "one is constant wherever the other is recorded"),
           ("~", "as many levels as samples carrying it")]
    y = TOP - 0.575
    for mark, what in key:
        lines = textwrap.wrap(what, 26)
        fig.text(RIGHT + 0.020, y, mark, fontsize=7, color=INK, ha="left",
                 va="top")
        fig.text(RIGHT + 0.044, y, "\n".join(lines), fontsize=5.4,
                 color=INK_2, ha="left", va="top", linespacing=1.35)
        y -= 0.014 * len(lines) + 0.014

    fig.suptitle(
        f"Fig. S10  Association between every pair of the {k} annotated "
        f"variables\n"
        f"{k * (k - 1) // 2} pairs in one symmetric matrix: above the diagonal "
        "the association, below it the samples the effect rests on",
        fontsize=8.6, x=LEFT - 0.02, ha="left", y=0.995, va="top",
        linespacing=1.6)

    _save(fig, "figS10_annotation_association")
    return fig


# ---------------------------------------------------------------- Fig. S11 / Table S13b
#
# The level grain of confounding.  Fig. S10 says two variables move together;
# it does not say which contrast that forbids, and the variable-level answer can
# be driven entirely by levels nobody would contrast.  Two choices shape it.
#
# *The baseline is the subset in which the biological variable is recorded,
# not the whole cohort.*  Race, ethnicity, vital status and age exist only for
# the TARGET neuroblastoma and osteosarcoma samples, and those are all polyA.
# Against a whole-cohort baseline every level of every one of those fields
# would read "100% polyA" -- a confound between *having a race recorded* and
# chemistry, which says nothing about the contrast a reader would actually run,
# within the annotated subset.  Comparing a level against the subset it is
# drawn from removes that artefact and leaves the real cases standing.
#
# *The composition is drawn, not thresholded.*  A count of levels past 90%
# answers "is this contrast blocked" and nothing else; the useful question is
# the one before it -- what does this level's chemistry, kit, run and library
# quality look like next to everything else annotated the same way.  So each
# level gets its full composition over the small-strata variables, its
# concentration in the large-strata ones, and a standardised difference on
# every numeric technical variable, against the same reference row.  The
# threshold survives as a mark on the cells that cross it.

def _compose(sub: pd.Series, levels: list[str]) -> np.ndarray:
    """Share of `sub` in each of `levels`, with the unrecorded share last."""
    out = np.zeros(len(levels) + 1)
    n = len(sub)
    if not n:
        out[-1] = 1.0
        return out
    for i, lv in enumerate(levels):
        out[i] = float((sub == lv).sum()) / n
    out[-1] = float(sub.isna().sum()) / n
    return out


def _confounding(meta: pd.DataFrame, biological: list[str],
                 technical: list[str]) -> dict:
    """Every biological level at n >= MIN_LEVEL_N against every technical variable.

    Each biological variable contributes a reference row (all samples carrying
    it) followed by one row per level.  A categorical technical variable gets
    the level's share in its largest stratum against the same share over the
    subset; a numeric one gets Cohen's d against the rest of the subset.
    """
    cat_bio = [b for b in biological if b not in NUMERIC and b not in REDUNDANT]
    cat_tech = [t for t in technical if t not in NUMERIC and t not in REDUNDANT]
    num_tech = [t for t in technical if t in NUMERIC and t not in REDUNDANT]

    rows: list[dict] = []
    for b in cat_bio:
        obs = meta.index[meta[b].notna()]
        if len(obs) < MIN_LEVEL_N:
            continue
        rows.append({"biological": b, "level": None, "idx": obs,
                     "n": len(obs), "ref": True})
        vc = meta.loc[obs, b].value_counts()
        for lv in vc.index:
            if int(vc[lv]) < MIN_LEVEL_N:
                continue
            rows.append({"biological": b, "level": str(lv),
                         "idx": meta.index[meta[b] == lv],
                         "n": int(vc[lv]), "ref": False})

    recs: list[dict] = []
    for r in rows:
        sub = meta.loc[r["idx"]]
        obs = meta.loc[meta[r["biological"]].notna()]
        rest = obs.drop(index=r["idx"], errors="ignore")
        for t in cat_tech:
            s = sub[t].dropna()
            base = obs[t].dropna()
            if not len(s) or base.nunique() < 2:
                continue
            share = s.value_counts(normalize=True)
            top = str(share.index[0])
            top_share = float(share.iloc[0])
            base_share = float((base == top).mean())
            recs.append({
                "biological": r["biological"], "level": r["level"] or "(all)",
                "is_reference_row": r["ref"], "n": r["n"],
                "technical": t, "measure": "share_in_largest_stratum",
                "stratum": top, "n_strata": int(s.nunique()),
                "value": round(top_share, 4),
                "baseline": round(base_share, 4),
                "n_level_annotated": int(len(s)),
                "subset_n": int(len(base)),
                "confounded": bool(not r["ref"]
                                   and len(s) >= MIN_LEVEL_N
                                   and top_share >= CONFOUND_SHARE
                                   and base_share < CONFOUND_SHARE)})
        for t in num_tech:
            x = sub[t].dropna().to_numpy(float)
            y = rest[t].dropna().to_numpy(float)
            if len(x) < MIN_LEVEL_N or len(y) < MIN_LEVEL_N:
                continue
            sd = np.sqrt(((len(x) - 1) * x.var(ddof=1)
                          + (len(y) - 1) * y.var(ddof=1))
                         / max(len(x) + len(y) - 2, 1))
            d = 0.0 if sd == 0 else float((x.mean() - y.mean()) / sd)
            recs.append({
                "biological": r["biological"], "level": r["level"] or "(all)",
                "is_reference_row": r["ref"], "n": r["n"],
                "technical": t, "measure": "cohens_d_vs_rest_of_subset",
                "stratum": "", "n_strata": 0,
                "value": round(d, 4), "baseline": 0.0,
                "n_level_annotated": int(len(x)),
                "subset_n": int(len(x) + len(y)),
                "confounded": bool(not r["ref"] and abs(d) >= BIG_D)})
    return {"rows": rows, "cells": pd.DataFrame(recs),
            "cat_tech": cat_tech, "num_tech": num_tech}


def table_s13b() -> pd.DataFrame:
    """Table S13b: confounding level by level -- every biological level against
    every technical variable, with the support and the call behind each cell."""
    meta, biological, technical = load_annotation()
    lvdf = _confounding(meta, biological, technical)["cells"]
    _write(lvdf, "tableS13b_confounding_levels.tsv", index=False)
    lv = lvdf[~lvdf.is_reference_row]
    print(f"   tableS13b: {lv.groupby(['biological', 'level']).ngroups} "
          f"biological levels, {len(lvdf)} cells; {int(lv['confounded'].sum())} "
          f"flagged across "
          f"{lvdf[lvdf['confounded']].groupby(['biological', 'level']).ngroups} "
          f"levels")
    return lvdf


@styled
def figure_s11() -> plt.Figure:
    """Fig. S11: where every biological level at n >= 8 sits in every technical
    variable, each against the samples in which its own field is recorded."""
    meta, biological, technical = load_annotation()
    c11 = _confounding(meta, biological, technical)
    rows, lvdf = c11["rows"], c11["cells"]
    cat_tech, num_tech = c11["cat_tech"], c11["num_tech"]
    stackable = [t for t in cat_tech
                 if meta[t].nunique() <= MAX_STRATA_TO_STACK]
    concentr = [t for t in cat_tech if t not in stackable]

    def cell(r: dict, t: str) -> dict | None:
        lv = r["level"] or "(all)"
        m = lvdf[(lvdf.biological == r["biological"]) & (lvdf.level == lv)
                 & (lvdf.is_reference_row == r["ref"]) & (lvdf.technical == t)]
        return None if m.empty else m.iloc[0].to_dict()

    # ---- layout, in inches ------------------------------------------------
    ROW_IN, LBL_IN = 0.152, 1.62
    STACK_IN, CONC_IN, HEAT_COL_IN, GAP_IN = 1.10, 1.02, 0.148, 0.20
    # The foot has to hold two things side by side: the rotated names of the
    # numeric variables under the heat map, and a legend under each composition
    # panel.  Both are sized here in inches, because a legend positioned as a
    # fraction of a body that grows with the number of levels ends up off the
    # page as soon as a level is added.
    LEG_BOT_IN, LEG_LINE_IN = 0.13, 0.125
    n_leg = max([meta[t].nunique() + 1 for t in stackable] or [1])
    TOP_IN = 0.86
    BOT_IN = max(0.98, LEG_BOT_IN + LEG_LINE_IN * (n_leg + 0.6) + 0.08)
    heat_in = HEAT_COL_IN * len(num_tech)
    body_in = (STACK_IN * len(stackable) + CONC_IN * len(concentr) + heat_in
               + GAP_IN * (len(stackable) + len(concentr)))
    fig_w = LBL_IN + body_in + 0.34
    fig_h = TOP_IN + ROW_IN * len(rows) + BOT_IN
    fig = plt.figure(figsize=(fig_w, fig_h))
    body_bot, body_h = BOT_IN / fig_h, (ROW_IN * len(rows)) / fig_h

    def panel(x_in: float, w_in: float):
        a = fig.add_axes((x_in / fig_w, body_bot, w_in / fig_w, body_h))
        a.set_ylim(len(rows) - 0.5, -0.5)
        a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
        return a

    # ---- row labels, grouped by biological variable -----------------------
    lab = panel(0.0, LBL_IN)
    lab.set_xlim(0, 1)
    lab.set_xticks([])
    for i, r in enumerate(rows):
        if r["ref"]:
            lab.text(0.985, i, f"{LABELS.get(r['biological'], r['biological'])}"
                               f"  —  all {r['n']} annotated",
                     ha="right", va="center", fontsize=5.6, color=INK,
                     fontweight="bold")
        else:
            lab.text(0.985, i, f"{r['level']}  ({r['n']})", ha="right",
                     va="center", fontsize=5.4, color=INK_2)
    axes = [lab]
    x = LBL_IN

    # ---- composition over the small-strata technical variables ------------
    for t in stackable:
        a = panel(x, STACK_IN)
        axes.append(a)
        levels = list(meta[t].value_counts().index.astype(str))
        cmap = _cat_colours(levels, QUAL)
        for i, r in enumerate(rows):
            comp = _compose(meta.loc[r["idx"], t].astype("object"), levels)
            left = 0.0
            for j, lv in enumerate(levels + ["not recorded"]):
                w = comp[j]
                if w <= 0:
                    continue
                a.barh(i, w, left=left, height=0.74,
                       color=cmap.get(lv, MISSING) if j < len(levels) else MISSING,
                       edgecolor="white", lw=0.25)
                left += w
            c = cell(r, t)
            if c and c["confounded"]:
                a.add_patch(Rectangle((0, i - 0.44), 1.0, 0.88, fill=False,
                                      edgecolor=INK, lw=0.9, zorder=5))
        a.set_xlim(0, 1)
        a.set_xticks([0, 0.5, 1])
        a.set_xticklabels(["0", "", "1"], fontsize=5)
        a.tick_params(length=1.5, pad=1)
        a.set_title(LABELS.get(t, t), fontsize=6.5, pad=4, color=INK)
        # a legend under each composition panel, so no colour is unnamed
        entries = levels + ["not recorded"]
        la = fig.add_axes((x / fig_w, LEG_BOT_IN / fig_h, STACK_IN / fig_w,
                           (LEG_LINE_IN * len(entries)) / fig_h))
        la.axis("off")
        la.set_xlim(0, 1)
        la.set_ylim(len(entries), 0)
        for k, lv in enumerate(entries):
            la.add_patch(Rectangle((0.01, k + 0.16), 0.075, 0.62,
                                   facecolor=cmap.get(lv, MISSING),
                                   edgecolor="white", lw=0.3))
            la.text(0.105, k + 0.5, lv, fontsize=4.9, color=INK_2,
                    va="center")
        x += STACK_IN + GAP_IN

    # ---- concentration in the many-strata technical variables -------------
    for t in concentr:
        a = panel(x, CONC_IN)
        axes.append(a)
        for i, r in enumerate(rows):
            c = cell(r, t)
            if c is None:
                continue
            v = c["value"]
            a.plot([0, v], [i, i], color=GRID, lw=0.8, zorder=1,
                   solid_capstyle="butt")
            a.scatter([v], [i], s=11, zorder=3,
                      color=INK if c["confounded"] else TYPE_COLOURS["TUMOR"],
                      edgecolor="white", linewidth=0.3)
            a.text(1.03, i, f"{c['n_strata']}", transform=a.get_yaxis_transform(),
                   fontsize=4.6, color=MUTED, va="center", ha="left",
                   clip_on=False)
            if c["confounded"]:
                a.text(v - 0.03, i, c["stratum"], fontsize=4.4, color=INK,
                       va="center", ha="right")
        a.axvline(CONFOUND_SHARE, color=INK_2, lw=0.6, ls=(0, (2, 2)))
        a.set_xlim(0, 1.0)
        a.set_xticks([0, 0.5, 1])
        a.set_xticklabels(["0", "0.5", "1"], fontsize=5)
        a.tick_params(length=1.5, pad=1)
        a.set_title(f"{LABELS.get(t, t)}\nshare on its largest run",
                    fontsize=6.5, pad=13, color=INK)
        a.text(1.03, -0.78, "runs", transform=a.get_yaxis_transform(),
               fontsize=4.6, color=MUTED, va="bottom", ha="left",
               clip_on=False)
        x += CONC_IN + GAP_IN

    # ---- standardised difference on every numeric technical variable ------
    hm = np.full((len(rows), len(num_tech)), np.nan)
    for i, r in enumerate(rows):
        for j, t in enumerate(num_tech):
            c = cell(r, t)
            if c is not None and not r["ref"]:
                hm[i, j] = c["value"]
    a = panel(x, heat_in)
    axes.append(a)
    # NaN is the reference rows, which have nothing to be different from.  It
    # gets the same grey the rest of the study uses for "not annotated", so an
    # empty row cannot be misread as a row of zero differences.
    div = plt.get_cmap("RdBu_r").copy()
    div.set_bad(MISSING)
    im = a.imshow(np.ma.masked_invalid(hm), aspect="auto", cmap=div,
                  vmin=-2, vmax=2, interpolation="nearest",
                  extent=(0, len(num_tech), len(rows) - 0.5, -0.5))
    for i in range(len(rows)):
        for j in range(len(num_tech)):
            if np.isfinite(hm[i, j]) and abs(hm[i, j]) >= BIG_D:
                a.add_patch(Rectangle((j + 0.08, i - 0.42), 0.84, 0.84,
                                      fill=False, edgecolor=INK, lw=0.6))
    a.set_xticks(np.arange(len(num_tech)) + 0.5)
    a.set_xticklabels([LABELS.get(t, t) for t in num_tech], fontsize=5,
                      rotation=90, ha="center")
    a.tick_params(length=1.5, pad=1)
    a.set_title("Numeric technical variables\n"
                "standardised difference from the rest of the same subset",
                fontsize=6.5, pad=4, color=INK)
    cax = fig.add_axes(((x + heat_in + 0.10) / fig_w, body_bot + body_h * 0.55,
                        0.075 / fig_w, body_h * 0.35))
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Cohen's d", fontsize=5)
    cb.ax.tick_params(labelsize=4.5, length=1.5)
    cb.outline.set_visible(False)

    # ---- separators between biological blocks, across every panel ---------
    for i, r in enumerate(rows):
        if r["ref"] and i:
            for aa in axes:
                aa.axhline(i - 0.5, color=INK_2, lw=0.6,
                           xmin=-0.02, xmax=1.02, clip_on=False)

    n_lv = lvdf[~lvdf.is_reference_row].groupby(
        ["biological", "level"]).ngroups
    fig.suptitle(
        "Where each biological level actually sits in each technical variable\n"
        f"{n_lv} levels at n ≥ {MIN_LEVEL_N}, each against the samples in "
        "which its own field is recorded (bold row). Outlined = "
        f"≥ {CONFOUND_SHARE:.0%} in one stratum while that subset is not, "
        f"or |d| ≥ {BIG_D}",
        fontsize=8, x=0.01, y=1 - 0.28 / fig_h, ha="left", va="top")

    _save(fig, "figS11_confounding")
    print(f"   figS11: {n_lv} biological levels x "
          f"{len(cat_tech) + len(num_tech)} technical variables")
    return fig


# ---------------------------------------------------------------- Table S14
def table_s14() -> pd.DataFrame:
    """Table S14: the chemistry-stratified test -- sample type on PC1 and PC4
    in both collections, unconditionally and within each chemistry stratum.

    An association between a component and `sample_type` and an association
    between the same component and `library_type` do not say which of the two
    the component is carrying.  This holds chemistry fixed and re-tests
    sample_type inside each stratum: a signal that survives is not chemistry, a
    signal that collapses was.

    The four recorded chemistries are collapsed to the two poles that separate
    them -- polyA selection against ribo-depletion.  `polya_stranded` is polyA
    with a stranded protocol and `ribominus` is a ribo-depletion kit; the
    stranded/kit distinction is not a pole, and splitting on it would leave
    groups below the size floor.
    """
    meta, _, _ = load_annotation()
    chem = meta.get("library_type")
    if chem is None or chem.isna().all():
        return pd.DataFrame()
    pole = chem.replace({"polya": "polya", "polya_stranded": "polya",
                         "ribozero": "ribozero", "ribominus": "ribozero"})
    rows = []
    for gs, d, stem in (("HALLMARK", HALL_PCA, "HALLMARK___all__"),
                        ("GOBP", ALL_PCA, "GOBP___all__")):
        sc = pd.read_csv(d / f"{stem}_pca_scores.tsv", sep="\t", index_col=0)
        for pc in ("PC1", "PC4"):
            v = sc[pc].reindex(meta.index)
            for label, cat, idx in (
                    ("library_type (run sheet)", meta["library_type"], meta.index),
                    ("chemistry (Brohl, 2 poles)", pole, meta.index),
                    ("sample_type", meta["sample_type"], meta.index),
                    ("sample_type | chemistry=ribozero", meta["sample_type"],
                     meta.index[pole == "ribozero"]),
                    ("sample_type | chemistry=polya", meta["sample_type"],
                     meta.index[pole == "polya"])):
                e, p, n = _eta_squared_raw(cat.reindex(idx), v.reindex(idx))
                rows.append({"collection": gs, "component": pc, "variable": label,
                             "eta_squared": e, "p_value": p, "n": n})
    out = pd.DataFrame(rows)
    out["fdr"] = benjamini_hochberg(out["p_value"].to_numpy())
    _write(out, "tableS14_chemistry_conditional.tsv", index=False)
    for gs in ("HALLMARK", "GOBP"):
        g = out[(out.collection == gs) & (out.component == "PC4")].set_index("variable")
        print(f"   {gs} PC4: sample_type {g.loc['sample_type', 'eta_squared']:.3f} "
              f"-> {g.loc['sample_type | chemistry=ribozero', 'eta_squared']:.3f} "
              f"(ribozero) / "
              f"{g.loc['sample_type | chemistry=polya', 'eta_squared']:.3f} (polya)")
    return out


# ---------------------------------------------------------------- shipping
def copy_annotation_data() -> list[str]:
    """Copy the processed annotation the tables derive from into submission/data/.

    Returns the names copied; a missing source file is reported, not skipped
    silently.
    """
    sdata = SUB / "data"
    sdata.mkdir(parents=True, exist_ok=True)
    copied, missing = [], []
    for name in SHIPPED_DATA:
        src = PROC / name
        if src.exists():
            (sdata / name).write_bytes(src.read_bytes())
            copied.append(name)
        else:
            missing.append(name)
    print(f"   submission/data: {len(copied)} of {len(SHIPPED_DATA)} "
          f"annotation tables copied")
    if missing:
        warnings.warn(f"not in {PROC}, so not shipped: {missing}", stacklevel=2)
    return copied


def run_all() -> None:
    """Write Figs. S4 and S9-S11 and Tables S11-S14, and ship the processed
    annotation.  Needs the pipeline runs and Tables S8-S9
    (:func:`analysis.figures.run_all`) first."""
    meta, biological, technical = load_annotation(verbose=True)
    print(f"   {len(meta)} samples; {len(biological)} biological and "
          f"{len(technical)} technical fields")
    for table in (table_s11, table_s12, table_s13a, table_s13b, table_s14):
        table()
    for draw in (figure_s9, figure_s4, figure_s10, figure_s11):
        plt.close(draw())
    copy_annotation_data()
