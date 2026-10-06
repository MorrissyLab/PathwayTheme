"""Diagnostic figures: scree, metadata-coloured PCA scatter, component
attribution, and theme summarisation.

Every function here draws onto an ``ax`` you supply, and each has a
``render_*`` wrapper that makes its own figure and saves it.  That split is
deliberate: the pipeline calls the wrappers to emit one file per plot, and a
manuscript can call the ``draw_*`` primitives into a GridSpec to compose a
multi-panel figure — without a second implementation of the same plot existing
anywhere.

The layout details here (where a bar label goes so it misses the cumulative
line, when heatmap text switches to white, what a blank cell means) are part of
the plot, not the caller's problem.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from .palettes import class_colour_lut
from .pca_figures import _safe_savefig

__all__ = [
    "draw_workflow", "render_workflow", "WORKFLOW_STAGES",
    "draw_scree", "render_scree",
    "draw_pc_loadings", "render_pc_loadings",
    "draw_pc_pairs", "render_pc_pairs",
    "draw_pca_scatter", "render_pca_scatter",
    "draw_metadata_attribution", "render_metadata_attribution",
    "draw_level_attribution", "render_level_attribution",
    "draw_theme_heatmap", "render_theme_heatmap",
    "draw_theme_reduction",
    "draw_effect_concordance", "render_effect_concordance",
]

_MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*", "<", ">")


def _short(term: str) -> str:
    """``GOBP_CELL_MOTILITY`` -> ``cell motility``."""
    for sep in ("_",):
        if sep in term:
            head = term.split(sep, 1)[0]
            if head.isupper() and len(head) > 1:
                term = term.split(sep, 1)[1]
                break
    return term.replace("_", " ").lower()


# ── the pipeline itself ───────────────────────────────────────────────────
#: One entry per stage: (number, name, what it consumes, what it does, lane).
#: ``lane`` is "trunk" for the shared step, then "unsupervised" / "supervised"
#: for the two branches off the score matrix.  Data-free -- it describes the
#: package, so it lives with the package rather than being redrawn in every
#: paper and doc that needs it.
WORKFLOW_STAGES: tuple[tuple[str, str, str, str, str], ...] = (
    ("1", "Enrichment", "gene × sample matrix, gene sets",
     "ssGSEA / Enrichr / GO-slim", "trunk"),
    ("2", "PCA of pathway space", "term × sample scores",
     "scores, loadings, variance", "unsupervised"),
    ("3", "Component annotation", "loadings, variance",
     "scree + top loadings per PC", "unsupervised"),
    ("4", "Metadata association", "PC scores, sample metadata",
     "Kruskal–Wallis / Spearman, BH", "unsupervised"),
    ("5", "Differential analysis", "term × sample scores, group labels",
     "moderated $t$, Welch or Mann–Whitney", "supervised"),
    ("6", "Theme summarisation", "differential results, term → theme map",
     "GO-slim, split by significance", "supervised"),
)

_INK = "#12304f"
_EDGE = "#4575b4"
_FILL = "#eaf0f7"
_KEY = "#b2182b"          # highlight outline only
_IN = "#3d5a6c"           # the "what it consumes" line
_LANE_FILL = {"unsupervised": "#f4f7fb", "supervised": "#fbf5f3"}


def _glyph_matrix(ax, x, y, w, h, *, rows=6, cols=5, shade=0.30,
                  seed=0) -> None:
    """A small matrix pictogram, so a data object reads as data not as a step."""
    from matplotlib.patches import Rectangle
    rng = np.random.default_rng(seed)
    cw, ch = w / cols, h / rows
    for r in range(rows):
        for c in range(cols):
            v = shade * rng.uniform(0.25, 1.0)
            ax.add_patch(Rectangle(
                (x + c * cw, y + r * ch), cw, ch, transform=ax.transAxes,
                facecolor=str(1.0 - v), edgecolor="white", linewidth=0.25,
                zorder=3))
    ax.add_patch(Rectangle((x, y), w, h, transform=ax.transAxes, fill=False,
                           edgecolor="0.45", linewidth=0.6, zorder=4))


def _arrow(ax, xy_from, xy_to, *, rad: float = 0.0, lw: float = 0.9) -> None:
    from matplotlib.patches import FancyArrowPatch
    ax.add_patch(FancyArrowPatch(
        xy_from, xy_to, transform=ax.transAxes, arrowstyle="-|>",
        mutation_scale=8, color="0.45", linewidth=lw, zorder=1,
        connectionstyle=f"arc3,rad={rad}", shrinkA=1, shrinkB=1))


def _workflow_vertical(ax, by_lane, box, content, lane_labels, aspect,
                       metrics) -> None:
    """The same pipeline running down the panel instead of across it.

    The trunk goes down the middle -- input, enrichment, score matrix -- and
    the two branches down two columns beneath it.  Everything below is in axes
    fractions, laid out by hand for the same reason the horizontal version is:
    the shape carries the meaning, and a layout engine does not know that the
    two lanes are alternatives rather than a sequence.
    """
    from matplotlib.patches import Rectangle

    PAD_X, PAD_Y, LANE_PT, _h = metrics
    W_LANE, W_GLYPH, H_GLYPH = 0.470, 0.070, 0.055
    GAP = 0.050                       # gap between stacked boxes
    # the boxes are as tall as the wordiest stage needs them to be, so a method
    # line that wraps to two finishes inside its border
    BH = max(content(n, s, W_LANE)[2] for _, n, _, s, _ in
             sum(by_lane.values(), [])) + 2 * PAD_Y
    X_L, X_R = 0.014, 0.516           # the two lane columns
    X_MID = 0.5
    Y_IN, Y_S1, Y_SM = 0.945, 0.775, 0.668
    Y_HUB = 0.640
    Y_LANE0 = 0.520                   # top of the first box in each lane
    Y_LANE_TOP = Y_LANE0 + 0.082      # top of the lane's background block

    # ── lane backgrounds, with the lane named on one line so the arrow that
    # lands on the block has somewhere to land that is not text ──────────
    for lane, x0 in (("unsupervised", X_L), ("supervised", X_R)):
        k = len(by_lane[lane])
        if not k:
            continue
        tall = k * BH + (k - 1) * GAP
        ax.add_patch(Rectangle(
            (x0 - 0.020, Y_LANE0 - tall - 0.026), W_LANE + 0.040,
            Y_LANE_TOP - (Y_LANE0 - tall - 0.026), transform=ax.transAxes,
            facecolor=_LANE_FILL[lane], edgecolor="none", zorder=0))
        if lane_labels:
            question = ("what is each group?" if lane == "unsupervised"
                        else "what separates them?")
            ax.text(x0 - 0.006, Y_LANE0 + 0.056, lane.upper(),
                    transform=ax.transAxes, ha="left", va="center",
                    fontsize=LANE_PT, fontweight="bold", color=_EDGE,
                    zorder=1)
            ax.text(x0 - 0.006, Y_LANE0 + 0.028, question,
                    transform=ax.transAxes, ha="left", va="center",
                    fontsize=LANE_PT - 0.4, color="0.45", style="italic",
                    zorder=1)

    # ── trunk, down the middle ───────────────────────────────────────────
    _glyph_matrix(ax, X_MID - W_GLYPH / 2, Y_IN, W_GLYPH, H_GLYPH, seed=1)
    ax.text(X_MID + W_GLYPH / 2 + 0.014, Y_IN + H_GLYPH / 2, "genes × samples",
            transform=ax.transAxes, ha="left", va="center", fontsize=5.2,
            color="0.35")
    _arrow(ax, (X_MID, Y_IN), (X_MID, Y_S1 + BH))

    num, name, _key, sub, _ = by_lane["trunk"][0]
    box(X_MID - W_LANE / 2, Y_S1, num, name, sub, w=W_LANE, h=BH)
    _arrow(ax, (X_MID, Y_S1), (X_MID, Y_SM + H_GLYPH))

    _glyph_matrix(ax, X_MID - W_GLYPH / 2, Y_SM, W_GLYPH, H_GLYPH,
                  shade=0.55, seed=2)
    ax.text(X_MID + W_GLYPH / 2 + 0.014, Y_SM + H_GLYPH / 2, "pathways × samples",
            transform=ax.transAxes, ha="left", va="center", fontsize=5.2,
            color=_INK, fontweight="bold")

    # ── the split: one stem down, two arms out ───────────────────────────
    ax.plot([X_MID, X_MID], [Y_SM, Y_HUB], transform=ax.transAxes,
            color="0.45", lw=1.0, zorder=1, solid_capstyle="round")
    ax.plot([X_MID], [Y_HUB], transform=ax.transAxes, marker="o", ms=3.0,
            color="0.45", zorder=2)
    for x0 in (X_L, X_R):
        _arrow(ax, (X_MID, Y_HUB), (x0 + W_LANE / 2, Y_LANE_TOP),
               rad=0.20 if x0 < X_MID else -0.20)

    # ── the two lanes, each running down ─────────────────────────────────
    for lane, x0 in (("unsupervised", X_L), ("supervised", X_R)):
        for i, (num, name, _key, sub, _) in enumerate(by_lane[lane]):
            y = Y_LANE0 - BH - i * (BH + GAP)
            box(x0, y, num, name, sub, w=W_LANE, h=BH)
            if i:
                _arrow(ax, (x0 + W_LANE / 2, y + BH + GAP),
                       (x0 + W_LANE / 2, y + BH))


def draw_workflow(ax, *, stages=WORKFLOW_STAGES,
                  panel_refs: Optional[dict] = None,
                  highlight: Iterable[str] = (),
                  lane_labels: bool = True,
                  title: Optional[str] = None,
                  orientation: str = "horizontal") -> None:
    """The pipeline as a data flow that branches at the score matrix.

    The shape is the architecture: one enrichment step produces a score matrix,
    from which an unsupervised branch (PCA, component annotation, metadata
    association) and a supervised branch (differential test, theme
    summarisation) run independently off the same cached scores.  Drawing it as
    a single 1-6 chain would imply an order that does not exist.

    ``panel_refs`` maps a stage number to a panel letter, drawn as a tag on that
    stage's box, so a manuscript can point each stage at the panel showing its
    output.  ``highlight`` outlines named stages for a figure about one part of
    the pipeline.

    ``orientation="vertical"`` runs the trunk down the middle and the two
    branches down two columns.  Horizontally the six boxes share the width, so
    each gets a sixth of it and the text inside shrinks with the panel; laid
    out down two columns each box gets nearly half, which is what makes it
    readable when the figure is placed at one column of a journal page.
    """
    if orientation not in ("horizontal", "vertical"):
        raise ValueError("orientation must be 'horizontal' or 'vertical'")
    from matplotlib.patches import FancyBboxPatch, Rectangle, Ellipse

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    hl = set(highlight)
    by_lane = {ln: [s for s in stages if s[4] == ln]
               for ln in ("trunk", "unsupervised", "supervised")}

    # ── geometry ─────────────────────────────────────────────────────────
    # Everything below is derived from the axes' real size in inches rather
    # than from axes fractions with hand-tuned constants.  A schematic drawn in
    # fractions only looks right at the one panel shape it was tuned on: put it
    # in a taller or narrower cell and the type stays the same size while the
    # boxes do not, so a wrapped method line walks out through the bottom
    # border and the number badge sits on top of it.  Sizing the box from the
    # text it has to hold makes the drawing correct at any panel shape.
    bb = ax.get_position()
    fw, fh = ax.figure.get_size_inches()
    ax_w_in = max(bb.width * fw, 1e-6)
    ax_h_in = max(bb.height * fh, 1e-6)
    aspect = ax_w_in / ax_h_in

    NAME_PT, SUB_PT, LANE_PT = 5.8, 5.0, 5.8
    PAD_X, PAD_Y = 0.050 / ax_w_in, 0.042 / ax_h_in

    def _h(pt: float, lead: float = 1.0) -> float:
        """A line of ``pt`` type as a fraction of the axes height."""
        return lead * (pt / 72.0) / ax_h_in

    def _fit(text: str, width: float, pt: float, bold: bool = False):
        """``text`` wrapped to what fits across ``width`` axes fractions."""
        em = (0.62 if bold else 0.60) * pt / 72.0
        return textwrap.wrap(text, max(6, int(width * ax_w_in / em))) or [""]

    R_Y = 0.75 * _h(NAME_PT)          # number badge, tied to the title type
    R_X = R_Y / aspect
    T_H, S_H = _h(NAME_PT, 1.34), _h(SUB_PT, 1.30)
    GAP_TS = 0.34 * S_H               # title block to method block

    def _content(name: str, sub: str, w: float):
        head = _fit(name, w - 2 * PAD_X, NAME_PT, bold=True)
        body = _fit(sub, w - 2 * PAD_X, SUB_PT)
        return head, body, len(head) * T_H + GAP_TS + len(body) * S_H

    def box(x, y, num, name, sub, w=None, h=None):
        w = W_BOX if w is None else w
        h = BH if h is None else h
        keyed = num in hl
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.002,rounding_size=0.016",
            transform=ax.transAxes, facecolor=_FILL,
            edgecolor=_KEY if keyed else _EDGE,
            linewidth=1.8 if keyed else 1.0, zorder=2))
        head, body, _total = _content(name, sub, w)
        # set from the top rather than centred: every box is as tall as the
        # wordiest stage needs, so centring would drop the title of a terse box
        # below the title of a wordy one standing beside it
        top = y + h - PAD_Y
        for i, ln in enumerate(head):
            ax.text(x + PAD_X, top - T_H * (i + 0.5), ln,
                    transform=ax.transAxes, ha="left", va="center",
                    fontsize=NAME_PT, fontweight="bold", color=_INK, zorder=3)
        y0 = top - len(head) * T_H - GAP_TS
        for i, ln in enumerate(body):
            ax.text(x + PAD_X, y0 - S_H * (i + 0.5), ln,
                    transform=ax.transAxes, ha="left", va="center",
                    fontsize=SUB_PT, color="0.34", style="italic", zorder=3)
        # the numbered badge rides the top-left corner rather than sitting
        # inside: the interior is the text's, and a badge in it either eats the
        # width the stage name needs or is drawn over the method line
        ax.add_patch(Ellipse((x, y + h), 2 * R_X, 2 * R_Y,
                             transform=ax.transAxes, facecolor=_EDGE,
                             edgecolor="white", linewidth=0.8, zorder=4))
        ax.text(x, y + h, num, transform=ax.transAxes, ha="center",
                va="center", fontsize=NAME_PT - 1.0, fontweight="bold",
                color="white", zorder=5)
        if panel_refs and num in panel_refs:
            # a tag on the box, not a caption under it: it says "this stage's
            # output is that panel", which is a property of the stage.  Sized
            # in y and corrected in x, so it stays square whatever shape the
            # panel holding the workflow happens to be
            th = 0.052 * min(1.0, h / 0.250)
            tw = th / aspect
            tx, ty = x + w - 0.010 - tw / 2, y + 0.012 + th / 2
            ax.add_patch(FancyBboxPatch(
                (tx - tw / 2, ty - th / 2), tw, th,
                boxstyle="round,pad=0.001,rounding_size=0.008",
                transform=ax.transAxes, facecolor="white",
                edgecolor=_EDGE, linewidth=0.7, zorder=4))
            ax.text(tx, ty, panel_refs[num], transform=ax.transAxes,
                    ha="center", va="center", fontsize=5.6,
                    fontweight="bold", color=_INK, zorder=5)

    if orientation == "vertical":
        _workflow_vertical(ax, by_lane, box, _content, lane_labels, aspect,
                           (PAD_X, PAD_Y, LANE_PT, _h))
        if title:
            ax.set_title(title, fontsize=9)
        return

    # ── the band's layout, solved rather than tuned ──────────────────────
    # Across, the drawing is trunk + three lanes' worth of box, so the box
    # width follows from the width there is; down, it is two lanes each with a
    # label strip, so the lane baselines follow from the tallest box.
    X_IN, W_GLYPH, GAP = 0.004, 0.042, 0.026
    LEAD = X_IN + 2 * W_GLYPH + 0.030 + 0.024 + 0.026 + 0.030
    W_BOX = (0.978 - LEAD - 2 * GAP) / 4
    X_S1 = X_IN + W_GLYPH + 0.030
    X_SM = X_S1 + W_BOX + 0.024
    X_BRANCH = X_SM + W_GLYPH + 0.026
    X_LANE0 = X_BRANCH + 0.030

    # the lane's name sits clear of the badge that rides the corner of
    # the first box under it, so the strip is a badge tall plus a line
    LS = (R_Y + _h(LANE_PT, 1.75)) if lane_labels else 0.014
    BH = max(_content(n, s, W_BOX)[2] for _, n, _, s, _ in stages) + 2 * PAD_Y
    BH = min(BH, (1.0 - 2 * LS - 0.07) / 2)
    slack = 1.0 - 2 * (BH + LS)
    lane_gap = max(0.020, min(0.10, slack * 0.5))
    top_pad = max(0.004, (slack - lane_gap) / 2)
    Y_UP = 1.0 - top_pad - LS - BH
    Y_DN = Y_UP - lane_gap - LS - BH
    Y_MID = (Y_UP + Y_DN) / 2
    Y_HUB = Y_MID + BH / 2

    # ── lane backgrounds, so the branch reads at a glance ────────────────
    for lane, ybase in (("unsupervised", Y_UP), ("supervised", Y_DN)):
        k = len(by_lane[lane])
        if not k:
            continue
        wide = k * W_BOX + (k - 1) * GAP
        pad = 0.5 * PAD_X
        ax.add_patch(Rectangle(
            (X_LANE0 - pad, ybase - PAD_Y), wide + 2 * pad,
            BH + LS + 2 * PAD_Y, transform=ax.transAxes,
            facecolor=_LANE_FILL[lane], edgecolor="none", zorder=0))
        if lane_labels:
            question = ("what is each group?" if lane == "unsupervised"
                        else "what separates them?")
            ly = ybase + BH + R_Y + _h(LANE_PT, 0.90)
            # the lane's name is set on the stage names' own left margin, and
            # the question it answers against the far edge of the lane, so the
            # two cannot meet however long either happens to be
            ax.text(X_LANE0 + PAD_X, ly, lane.upper(), transform=ax.transAxes,
                    ha="left", va="center", fontsize=LANE_PT,
                    fontweight="bold", color=_EDGE, zorder=1)
            ax.text(X_LANE0 + wide - PAD_X, ly, question,
                    transform=ax.transAxes, ha="right", va="center",
                    fontsize=LANE_PT, color="0.45", style="italic", zorder=1)

    # ── trunk: input matrix -> enrichment -> score matrix ────────────────
    G_H = 0.52 * BH
    G_Y = Y_HUB - G_H / 2 + 0.030
    _glyph_matrix(ax, X_IN, G_Y, W_GLYPH, G_H, seed=1)
    ax.text(X_IN + W_GLYPH / 2, G_Y - 0.012, "genes ×\nsamples",
            transform=ax.transAxes, ha="center", va="top", fontsize=SUB_PT,
            color="0.35", linespacing=1.25)
    _arrow(ax, (X_IN + W_GLYPH, Y_HUB), (X_S1, Y_HUB))

    num, name, _key, sub, _ = by_lane["trunk"][0]
    box(X_S1, Y_MID, num, name, sub)
    _arrow(ax, (X_S1 + W_BOX, Y_HUB), (X_SM, Y_HUB))

    _glyph_matrix(ax, X_SM, G_Y, W_GLYPH, G_H, shade=0.55, seed=2)
    ax.text(X_SM + W_GLYPH / 2, G_Y - 0.012, "pathways ×\nsamples",
            transform=ax.transAxes, ha="center", va="top", fontsize=SUB_PT,
            color=_INK, linespacing=1.25, fontweight="bold")

    # ── the split, drawn as one stem into two curved arms ────────────────
    ax.plot([X_SM + W_GLYPH, X_BRANCH], [Y_HUB, Y_HUB],
            transform=ax.transAxes, color="0.45", lw=1.0, zorder=1,
            solid_capstyle="round")
    ax.plot([X_BRANCH], [Y_HUB], transform=ax.transAxes, marker="o", ms=3.0,
            color="0.45", zorder=2)
    for ybase in (Y_UP, Y_DN):
        _arrow(ax, (X_BRANCH, Y_HUB), (X_LANE0, ybase + BH / 2),
               rad=0.22 if ybase > Y_HUB else -0.22)
    

    # ── the two lanes ────────────────────────────────────────────────────
    for lane, ybase in (("unsupervised", Y_UP), ("supervised", Y_DN)):
        for i, (num, name, _key, sub, _) in enumerate(by_lane[lane]):
            x = X_LANE0 + i * (W_BOX + GAP)
            box(x, ybase, num, name, sub)
            if i:
                _arrow(ax, (x - GAP, ybase + BH / 2), (x, ybase + BH / 2))

    if title:
        ax.set_title(title, fontsize=8)


def render_workflow(out_path: str | Path, **kw) -> Path:
    dpi = kw.pop("dpi", 150)
    fig, ax = plt.subplots(figsize=(10.0, 2.2))
    # full bleed: the box widths are axes fractions, and the default margins
    # would make a standalone render narrower than the same drawing composed
    # into a manuscript GridSpec, overflowing the stage names
    fig.subplots_adjust(left=0.012, right=0.988, top=0.985, bottom=0.015)
    draw_workflow(ax, **kw)
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


# ── scree / elbow ─────────────────────────────────────────────────────────
def draw_scree(ax, variance_explained: Sequence[float], *,
               elbow: Optional[int] = None, annotate: int = 3,
               mark: Optional[dict] = None,
               mark_label: str = "",
               title: Optional[str] = None) -> None:
    """Variance per component with its cumulative total.

    ``elbow`` draws a marker after that component; pass the number of components
    you intend to retain.  ``annotate`` labels that many leading bars.

    ``mark`` maps a component number to a short label and colours those bars,
    so the panel can say which components carry a result rather than leaving
    the reader to infer it from the height of a bar.  A scree plot read alone
    says "keep three"; where the components past the elbow turn out to carry
    something, that is the more useful thing for it to say.
    """
    ve = np.asarray(variance_explained, dtype=float) * 100.0
    cum = np.cumsum(ve)
    pcs = np.arange(1, len(ve) + 1)

    marked = dict(mark or {})
    colours = ["#b2182b" if i in marked else "#7f9fc4" for i in pcs]
    ax.bar(pcs, ve, color=colours, width=0.66, zorder=2)
    ax.plot(pcs, cum, "o-", color="#b2182b", ms=3, lw=1.2, zorder=3,
            label="cumulative")
    if elbow:
        ax.axvline(elbow + 0.5, color="0.55", ls=":", lw=1.0, zorder=1)
        ax.text(elbow + 0.6, 92, "elbow", fontsize=6, color="0.35", va="top")
    if annotate:
        # A text block, not per-bar labels: leading bars are adjacent so their
        # labels collide, and an above-bar label runs into the cumulative line,
        # which passes exactly through the top of the first bar.
        n = min(annotate, len(ve))
        ax.text(0.97, 0.74,
                "\n".join(f"PC{i + 1}  {v:.1f}%" for i, v in enumerate(ve[:n]))
                + f"\ncum.  {cum[n - 1]:.1f}%",
                transform=ax.transAxes, fontsize=6.5, ha="right", va="top",
                linespacing=1.5, family="monospace", color="0.2")
    if marked:
        # under the axis, not over the bars: the marked components are the
        # short ones, and a label above a 2% bar is a label in the middle of
        # the panel
        for i, what in marked.items():
            ax.annotate(what, xy=(i, 0), xytext=(i, -0.155),
                        textcoords=("data", "axes fraction"), fontsize=5.4,
                        color="#b2182b", ha="center", va="top", rotation=90,
                        annotation_clip=False)
    ax.set_xticks(pcs)
    ax.set_xticklabels(pcs, fontsize=7)
    ax.set_ylim(0, 100)
    # the marks are set under the axis and run vertically, so the axis label
    # has to clear the longest of them rather than a fixed distance
    pad = 4.0
    if marked:
        pad = 9.0 + 0.56 * 5.4 * max(len(str(v)) for v in marked.values())
    ax.set_xlabel("component", fontsize=8, labelpad=pad)
    ax.set_ylabel("% of variance", fontsize=8)
    ax.tick_params(labelsize=7)
    handles, labels = ax.get_legend_handles_labels()
    if marked and mark_label:
        handles.append(Line2D([], [], marker="s", ls="", color="#b2182b",
                              ms=4, label=mark_label))
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="lower right",
              bbox_to_anchor=(1.0, 0.02 if marked else 0.10))
    if title:
        ax.set_title(title, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)


def draw_pc_loadings(ax, loadings: pd.Series, *,
                     top_n: int = 10,
                     variance_explained: Optional[float] = None,
                     component: str = "PC1",
                     wrap: int = 46,
                     label_fontsize: float = 5.6,
                     title: Optional[str] = None) -> None:
    """The strongest signed loadings of one component, as a diverging bar plot.

    ``loadings`` is that component's loading per pathway (the row of the
    loadings matrix).  The top ``top_n`` in each direction are drawn, positive
    above negative, so a component is read as the programmes that define it and
    the two ends are named separately.  A component whose loadings are one-signed
    -- which happens on a leading component of a score matrix -- simply has no
    bars on that side, and the panel says so rather than inventing an opposite
    end by rank.
    """
    lv = pd.Series(loadings).dropna().astype(float)
    up = lv[lv > 0].sort_values(ascending=False).head(top_n)
    dn = lv[lv < 0].sort_values().head(top_n)
    # drawn bottom-up: most negative at the bottom, most positive at the top
    sel = pd.concat([dn.sort_values(ascending=False), up.sort_values()])
    if sel.empty:
        ax.set_axis_off()
        return

    y = np.arange(len(sel))
    colours = ["#b2182b" if v > 0 else "#2166ac" for v in sel.values]
    ax.barh(y, sel.values, color=colours, height=0.72, zorder=2)
    ax.axvline(0.0, color="0.25", lw=0.8, zorder=3)
    ax.set_yticks(y)
    # a long term name has to wrap, and at this size a wrapped label runs into
    # its neighbour unless the wrapped lines are pulled together
    ax.set_yticklabels(["\n".join(textwrap.wrap(_short(t), wrap))
                        for t in sel.index], fontsize=label_fontsize,
                       linespacing=0.92)
    ax.set_ylim(-0.7, len(sel) - 0.3)
    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("loading", fontsize=7)
    if title is None:
        title = component
        if variance_explained is not None:
            title += f" ({variance_explained * 100:.1f}%)"
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="0.9", lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def render_pc_loadings(loadings_df: pd.DataFrame, out_path: str | Path, *,
                       components: Sequence[str] = ("PC1", "PC2", "PC3", "PC4"),
                       variance_explained: Optional[Sequence[float]] = None,
                       top_n: int = 10, dpi: int = 150, **kw) -> Path:
    """One ``draw_pc_loadings`` panel per component, stacked in a column."""
    comps = [c for c in components if c in loadings_df.index]
    ve = (dict(zip(loadings_df.index, np.asarray(variance_explained, float)))
          if variance_explained is not None else {})
    fig, axes = plt.subplots(len(comps), 1,
                             figsize=(5.0, 2.55 * top_n / 10 * len(comps)))
    for ax, c in zip(np.atleast_1d(axes), comps):
        draw_pc_loadings(ax, loadings_df.loc[c], component=c, top_n=top_n,
                         variance_explained=ve.get(c), **kw)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def draw_pc_pairs(fig, scores: pd.DataFrame,
                  variance_explained: Sequence[float], *,
                  n_components: int = 10,
                  colour_by: Optional[pd.Series] = None,
                  colours: Optional[dict] = None,
                  subplot_spec=None,
                  point_size: float = 1.6,
                  legend: bool = True,
                  legend_loc: tuple = (0.80, 0.80),
                  title: Optional[str] = None):
    """Every pair of the leading components in one triangular grid.

    Only the lower triangle is drawn, with a per-group density on the diagonal:
    the upper triangle of a pairs plot is the same panels transposed, and
    printing both doubles the ink without adding a comparison.

    This is the exhaustive visual scan that the attribution grid replaces --
    ``n_components`` of 10 is 45 panels to read by eye, against one row per
    component of a corrected test.
    """
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    ve = np.asarray(variance_explained, dtype=float)
    comps = [f"PC{i + 1}" for i in range(n_components) if f"PC{i + 1}" in scores]
    k = len(comps)

    cser = (colour_by.reindex(scores.index) if colour_by is not None
            else pd.Series("all", index=scores.index))
    groups = [g for g in pd.unique(cser.dropna())]
    lut = dict(colours) if colours else class_colour_lut(sorted(map(str, groups)))
    lut = {g: lut.get(str(g), lut.get(g, "#4575b4")) for g in groups}

    kws = dict(wspace=0.10, hspace=0.10)
    gs = (GridSpecFromSubplotSpec(k, k, subplot_spec=subplot_spec, **kws)
          if subplot_spec is not None else GridSpec(k, k, figure=fig, **kws))

    axes = {}
    for r in range(k):
        for c in range(k):
            if c > r:
                continue
            ax = fig.add_subplot(gs[r, c])
            axes[(r, c)] = ax
            ax.tick_params(labelsize=4.6, length=1.6, pad=1.0)
            for sp in ax.spines.values():
                sp.set_linewidth(0.5)
                sp.set_color("0.6")
            if r == c:
                # a density per group, so the diagonal carries the marginal the
                # scatter panels flatten
                v_all = scores[comps[r]].astype(float)
                lo, hi = float(v_all.min()), float(v_all.max())
                grid = np.linspace(lo, hi, 200)
                for g in groups:
                    v = v_all[cser == g].dropna()
                    if len(v) < 5:
                        continue
                    from scipy.stats import gaussian_kde
                    d = gaussian_kde(v.values)(grid)
                    ax.fill_between(grid, d, color=lut[g], alpha=0.45, lw=0.0)
                    ax.plot(grid, d, color=lut[g], lw=0.6)
                ax.set_yticks([])
                ax.set_xlim(lo, hi)
            else:
                for g in groups:
                    sel = cser == g
                    if not sel.any():
                        continue
                    ax.scatter(scores.loc[sel, comps[c]],
                               scores.loc[sel, comps[r]], c=lut[g],
                               s=point_size, alpha=0.65, linewidths=0.0,
                               zorder=2)
            if c != 0 or r == 0:
                ax.set_yticklabels([])
            if r != k - 1:
                ax.set_xticklabels([])
            if c == 0:
                # the top-left cell is PC1's own density, and it is the only
                # place PC1 can be named on the vertical axis
                ax.set_ylabel(f"{comps[r]}\n{ve[r] * 100:.1f}%", fontsize=5.4)
            if r == k - 1:
                ax.set_xlabel(f"{comps[c]}\n{ve[c] * 100:.1f}%", fontsize=5.4)

    if legend and colour_by is not None:
        handles = [Line2D([], [], marker="o", ls="", color=lut[g], ms=5,
                          label=str(g)) for g in groups]
        fig.legend(handles=handles, fontsize=7, frameon=False,
                   loc="center", bbox_to_anchor=legend_loc,
                   title=str(colour_by.name or ""), title_fontsize=7.5)
    if title:
        fig.suptitle(title, fontsize=9)
    return axes


def render_pc_pairs(scores, variance_explained, out_path: str | Path, *,
                    n_components: int = 10, dpi: int = 200, **kw) -> Path:
    side = 0.78 * n_components + 1.4
    fig = plt.figure(figsize=(side, side))
    draw_pc_pairs(fig, scores, variance_explained,
                  n_components=n_components, **kw)
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def render_scree(variance_explained, out_path: str | Path, *,
                 elbow: Optional[int] = None, title: Optional[str] = None,
                 dpi: int = 150) -> Path:
    fig, ax = plt.subplots(figsize=(4.8, 3.3))
    draw_scree(ax, variance_explained, elbow=elbow, title=title)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


# ── PCA scatter coloured by metadata ──────────────────────────────────────
def draw_pca_scatter(ax, scores: pd.DataFrame,
                     variance_explained: Sequence[float], *,
                     colour_by: Optional[pd.Series] = None,
                     marker_by: Optional[pd.Series] = None,
                     components: tuple[str, str] = ("PC1", "PC2"),
                     colours: Optional[dict] = None,
                     markers: Optional[dict] = None,
                     group_means: bool = False,
                     min_group_for_mean: int = 3,
                     point_size: float = 13.0,
                     legend: bool = True,
                     legend_loc: str = "lower left",
                     title: Optional[str] = None) -> None:
    """Two components of ``scores``, coloured and optionally shaped by metadata.

    ``components`` selects which pair to draw, so the biological and the
    technical axis of one decomposition can be shown side by side.

    ``group_means=True`` draws a dashed line at each colour group's mean on the
    first component and labels it, which is how an ordering claim about that
    component is shown rather than asserted.  Groups below
    ``min_group_for_mean`` are omitted from that -- a line at a single sample is
    not a group position.
    """
    x, y = components
    ve = np.asarray(variance_explained, dtype=float)
    pc_index = {f"PC{i + 1}": i for i in range(len(ve))}

    cser = (colour_by.reindex(scores.index) if colour_by is not None
            else pd.Series("all", index=scores.index))
    groups = [g for g in pd.unique(cser.dropna())]
    lut = dict(colours) if colours else class_colour_lut(sorted(map(str, groups)))
    lut = {g: lut.get(str(g), lut.get(g, "#4575b4")) for g in groups}

    if marker_by is not None:
        mser = marker_by.reindex(scores.index)
        mlevels = sorted(pd.unique(mser.dropna()), key=str)
        mmap = dict(markers) if markers else {}
        for i, m in enumerate(mlevels):
            mmap.setdefault(m, _MARKERS[i % len(_MARKERS)])
    else:
        mser, mlevels, mmap = None, [], {}

    for g in groups:
        for m in (mlevels or [None]):
            sel = cser == g
            if m is not None:
                sel = sel & (mser == m)
            if not sel.any():
                continue
            ax.scatter(scores.loc[sel, x], scores.loc[sel, y], c=lut[g],
                       marker=mmap.get(m, "o"), s=point_size, alpha=0.75,
                       linewidths=0.25, edgecolors="white", zorder=2)

    if group_means and len(groups) > 1:
        # A single observation has no mean worth drawing, and a group of two
        # states an ordering the data cannot support -- skip below the floor
        # rather than draw a line the reader will read as a group position.
        means = {}
        for g in groups:
            v = scores.loc[cser == g, x].dropna()
            if len(v) >= min_group_for_mean:
                means[g] = float(v.mean())
        if means:
            # headroom first, so the labels sit inside the axes not clipped
            lo, hi = ax.get_ylim()
            span = hi - lo
            ax.set_ylim(lo, hi + 0.16 * span)
            # Two nearby groups put their labels on top of each other, which is
            # exactly the case an ordering claim cares about.  Stagger by rank
            # along the axis so neighbours land on different rows.
            order = sorted(means, key=means.get)
            # Two rows separate two neighbours, but four groups on one axis put
            # the first and third back on the same row -- which is where the
            # chemistry means collided.  Give every group its own row when
            # there are more than three.
            n_rows = 1 if len(order) <= 2 else (2 if len(order) <= 3
                                                else len(order))
            rows = tuple(0.985 - 0.055 * i for i in range(n_rows))
            for i, g in enumerate(order):
                mu = means[g]
                ax.axvline(mu, color=lut[g], lw=1.1, ls="--", alpha=0.85,
                           zorder=1)
                # each label sits on its own group's dashed line, and with
                # several groups it also sits across its neighbours' -- a plain
                # white backing keeps the digits readable through both
                ax.annotate(f"{mu:+.1f}", xy=(mu, rows[i % len(rows)]),
                            xycoords=("data", "axes fraction"), color=lut[g],
                            fontsize=6.5, fontweight="bold", ha="center",
                            va="top", zorder=4,
                            bbox=dict(boxstyle="square,pad=0.12", fc="white",
                                      ec="none", alpha=0.85))

    ax.set_xlabel(f"{x} ({ve[pc_index[x]] * 100:.1f}%)", fontsize=8)
    ax.set_ylabel(f"{y} ({ve[pc_index[y]] * 100:.1f}%)", fontsize=8)
    ax.tick_params(labelsize=7)
    if title:
        ax.set_title(title, fontsize=8)
    if legend:
        chandles = [Line2D([], [], marker="o", ls="", color=lut[g], ms=5,
                           label=str(g)) for g in groups]
        mhandles = [Line2D([], [], marker=mmap[m], ls="", color="0.45", ms=4,
                           label=str(m)) for m in mlevels]
        if chandles and mhandles:
            # Matplotlib fills a multi-column legend column-major, so unequal
            # lists spill markers into the colour column -- where a shape reads
            # as a colour group.  Pad both to the same length with blanks.
            n = max(len(chandles), len(mhandles))
            blank = lambda: Line2D([], [], ls="", marker="", label=" ")
            chandles += [blank() for _ in range(n - len(chandles))]
            mhandles += [blank() for _ in range(n - len(mhandles))]
            ax.legend(handles=chandles + mhandles, fontsize=5.5, frameon=False,
                      ncol=2, loc=legend_loc, handletextpad=0.3,
                      columnspacing=0.8)
        elif chandles or mhandles:
            ax.legend(handles=chandles or mhandles, fontsize=5.5,
                      frameon=False, loc=legend_loc, handletextpad=0.3)
    ax.spines[["top", "right"]].set_visible(False)


def render_pca_scatter(scores, variance_explained, out_path: str | Path,
                       **kw) -> Path:
    dpi = kw.pop("dpi", 150)
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    draw_pca_scatter(ax, scores, variance_explained, **kw)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


# ── component x metadata attribution ──────────────────────────────────────
def draw_metadata_attribution(ax, assoc: pd.DataFrame, *,
                              variance_explained: Optional[Sequence[float]] = None,
                              technical: Iterable[str] = (),
                              n_components: Optional[int] = None,
                              column_order: Optional[Sequence[str]] = None,
                              vmax: float = 0.9,
                              show_values: bool = True,
                              colorbar: bool = True,
                              colorbar_horizontal: bool = False,
                              title: Optional[str] = None) -> None:
    """Components x metadata variables, cell = effect size.

    ``assoc`` is the long table from
    :func:`pathwaytheme.pca.metadata.associate_metadata`.  An asterisk marks
    ``significant``; names in ``technical`` are labelled in red so a technical
    source is visually separable from a biological one.
    """
    if assoc.empty:
        ax.set_axis_off()
        return
    wide = assoc.pivot_table(index="component", columns="variable",
                             values="effect")
    sig = assoc.pivot_table(index="component", columns="variable",
                            values="significant", aggfunc="first")
    pcs = sorted(wide.index, key=lambda s: int(str(s).lstrip("PC") or 0))
    if n_components:
        pcs = pcs[:n_components]
    cols = ([c for c in column_order if c in wide.columns] if column_order
            else list(wide.columns))
    wide = wide.loc[pcs, cols]
    sig = sig.reindex(index=pcs, columns=cols).fillna(False)

    im = ax.imshow(wide.to_numpy(dtype=float), cmap="magma_r", vmin=0,
                   vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("_", " ") for c in cols], rotation=45,
                       ha="right", fontsize=6.5)
    tech = set(technical)
    for lbl, name in zip(ax.get_xticklabels(), cols):
        if name in tech:
            lbl.set_color("#b2182b")
    if variance_explained is not None:
        ve = np.asarray(variance_explained, dtype=float)
        ylabels = [f"{p} ({ve[i] * 100:.0f}%)" if i < len(ve) else str(p)
                   for i, p in enumerate(pcs)]
    else:
        ylabels = [str(p) for p in pcs]
    ax.set_yticks(range(len(pcs)))
    ax.set_yticklabels(ylabels, fontsize=6.5)
    if show_values:
        for i in range(len(pcs)):
            for j in range(len(cols)):
                v = wide.iat[i, j]
                if pd.isna(v):
                    continue
                star = "*" if bool(sig.iat[i, j]) else ""
                ax.text(j, i, f"{star}{v:.2f}", ha="center", va="center",
                        fontsize=5.2, color="white" if v > 0.55 * vmax else "0.15")
    if title:
        ax.set_title(title, fontsize=8)
    if colorbar:
        # horizontal when the panel has a neighbour: a vertical bar lands on the
        # next panel's y labels
        if colorbar_horizontal:
            cb = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.85,
                              pad=0.30, aspect=28)
        else:
            cb = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
        cb.set_label("effect size", fontsize=7)
        cb.ax.tick_params(labelsize=6)


def draw_level_attribution(ax, levels: pd.DataFrame, *,
                           variable: Optional[str] = None,
                           variance_explained: Optional[Sequence[float]] = None,
                           n_components: Optional[int] = None,
                           level_order: Optional[Sequence[str]] = None,
                           vmax: float = 1.0,
                           show_values: bool = True,
                           show_n: bool = True,
                           colorbar: bool = True,
                           colorbar_horizontal: bool = False,
                           title: Optional[str] = None) -> None:
    """Components x individual levels, cell = signed one-versus-rest effect.

    ``levels`` is the long table from
    :func:`pathwaytheme.pca.metadata.associate_levels`.  The scale is diverging
    because the sign is the reading: red means samples of that level sit at the
    positive end of the component, blue the negative end, white means the level
    is not separated from the rest.  An asterisk marks ``significant``.

    Where :func:`draw_metadata_attribution` answers "which components carry the
    diagnosis", this answers "which component is the CCSK axis".
    """
    from ..pca.metadata import level_matrix

    wide = level_matrix(levels, variable=variable, signed=True)
    if wide.empty:
        ax.set_axis_off()
        return
    df = levels if variable is None else levels[levels["variable"] == variable]
    sigl = df.copy()
    sigl["_col"] = (sigl["level"] if variable is not None
                    else sigl["variable"] + ": " + sigl["level"])
    sig = sigl.pivot_table(index="component", columns="_col",
                           values="significant", aggfunc="first")

    pcs = list(wide.index)
    if n_components:
        pcs = pcs[:n_components]
    cols = ([c for c in level_order if c in wide.columns] if level_order
            else list(wide.columns))
    wide = wide.loc[pcs, cols]
    sig = sig.reindex(index=pcs, columns=cols).fillna(False)

    im = ax.imshow(wide.to_numpy(dtype=float), cmap="RdBu_r", vmin=-vmax,
                   vmax=vmax, aspect="auto")
    # a level of 13 samples separates perfectly far more easily than one of
    # 200, so the group size belongs on the axis, not in a caption
    n_of = (sigl.groupby("_col")["n_level"].first().to_dict() if show_n
            and "n_level" in sigl.columns else {})
    xlabels = [c.replace("_", " ") + (f"  (n={int(n_of[c])})" if c in n_of else "")
               for c in cols]
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=6.5)
    if variance_explained is not None:
        ve = np.asarray(variance_explained, dtype=float)
        ylabels = [f"{p} ({ve[i] * 100:.0f}%)" if i < len(ve) else str(p)
                   for i, p in enumerate(pcs)]
    else:
        ylabels = [str(p) for p in pcs]
    ax.set_yticks(range(len(pcs)))
    ax.set_yticklabels(ylabels, fontsize=6.5)
    if show_values:
        for i in range(len(pcs)):
            for j in range(len(cols)):
                v = wide.iat[i, j]
                if pd.isna(v):
                    continue
                star = "*" if bool(sig.iat[i, j]) else ""
                ax.text(j, i, f"{star}{v:+.2f}", ha="center", va="center",
                        fontsize=5.0,
                        color="white" if abs(v) > 0.62 * vmax else "0.15")
    if title:
        ax.set_title(title, fontsize=8)
    if colorbar:
        if colorbar_horizontal:
            cb = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.85,
                              pad=0.30, aspect=28)
        else:
            cb = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
        cb.set_label("rank-biserial, level vs rest", fontsize=7)
        cb.ax.tick_params(labelsize=6)


def render_level_attribution(levels: pd.DataFrame, out_path: str | Path,
                             **kw) -> Path:
    dpi = kw.pop("dpi", 150)
    var = kw.get("variable")
    sub = levels if var is None or levels.empty else levels[levels["variable"] == var]
    n_col = max(sub["level"].nunique() if not sub.empty else 1, 1)
    if var is None and not levels.empty:
        n_col = levels.groupby("variable")["level"].nunique().sum()
    n_pc = kw.get("n_components") or (levels["component"].nunique()
                                      if not levels.empty else 1)
    fig, ax = plt.subplots(figsize=(2.0 + 0.60 * n_col, 1.4 + 0.42 * n_pc))
    draw_level_attribution(ax, levels, **kw)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def render_metadata_attribution(assoc: pd.DataFrame, out_path: str | Path,
                                **kw) -> Path:
    dpi = kw.pop("dpi", 150)
    n_var = assoc["variable"].nunique() if not assoc.empty else 1
    n_pc = kw.get("n_components") or (assoc["component"].nunique()
                                      if not assoc.empty else 1)
    fig, ax = plt.subplots(figsize=(1.6 + 0.66 * n_var, 1.4 + 0.42 * n_pc))
    draw_metadata_attribution(ax, assoc, **kw)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


# ── theme summarisation heatmap ─────────────────────────────────────────────────
def draw_theme_reduction(ax, term_effects: Sequence[float],
                         theme_effects: pd.Series, *,
                         theme_sizes: Optional[pd.Series] = None,
                         annotate: int = 2,
                         label_chars: int = 26,
                         title: Optional[str] = None) -> None:
    """One contrast at both grains on one effect axis: terms above, themes below.

    The summarisation stage is a claim about *compression* — thousands of terms
    become tens of themes — and neither grain drawn alone makes that claim.  A
    panel of theme values shows the destination without the journey, and a panel
    of the strongest individual terms shows the thing being replaced.  Sharing
    one axis puts the two populations against each other, which is the only
    arrangement in which the reduction is the thing you see.

    ``term_effects`` are the individual terms that cleared significance;
    ``theme_effects`` is one mean per theme, and ``theme_sizes`` how many terms
    each stands for, drawn as the marker area so a theme cannot look important
    on the strength of one member.
    """
    x = np.asarray(term_effects, dtype=float)
    x = x[np.isfinite(x)]
    th = theme_effects.dropna()
    if not len(x) or not len(th):
        ax.set_axis_off()
        return

    ax.scatter(x, np.ones(len(x)), marker="|", s=90, linewidths=0.5,
               color="0.62", alpha=0.30, zorder=1)

    if theme_sizes is not None:
        n = theme_sizes.reindex(th.index).fillna(1.0).to_numpy(dtype=float)
        size = 12.0 + 90.0 * np.sqrt(n / max(n.max(), 1.0))
    else:
        size = np.full(len(th), 34.0)
    ax.scatter(th.to_numpy(), np.zeros(len(th)), s=size,
               color=["#b2182b" if v > 0 else "#2166ac" for v in th],
               edgecolor="white", linewidth=0.4, zorder=3)

    # the extremes carry the reading, and only they are named: a label per
    # theme would be sixty-eight labels on a panel two inches wide
    if annotate:
        # the named themes sit at the two ends, so each end's labels are near
        # neighbours: they are stacked at alternating depths rather than set at
        # one, which is what made them overprint each other
        k = min(annotate, len(th) // 2)
        lo, hi = float(np.min(th)), float(np.max(th))
        for side in (th.nsmallest(k), th.nlargest(k)):
            for i, (name, v) in enumerate(side.items()):
                ax.annotate("\n".join(textwrap.wrap(str(name), label_chars)),
                            xy=(v, 0), xytext=(v, -0.52 - 0.62 * (i % 2)),
                            textcoords="data", fontsize=5.8, color="0.25",
                            ha="left" if v > (lo + hi) / 2 else "right",
                            va="top",
                            arrowprops=dict(arrowstyle="-", color="0.7",
                                            lw=0.6, shrinkA=0, shrinkB=3))
    ax.axvline(0, color="0.45", lw=0.8, zorder=0)
    # the deeper of the two label rows runs to three wrapped lines, and the
    # axis has to leave room for it rather than let it cross the spine
    ax.set_ylim(-2.55 if annotate else -0.5, 1.6)
    ax.set_yticks([1, 0])
    ax.set_yticklabels([f"{len(x):,} terms", f"{len(th)} themes"], fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("effect", fontsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    if title:
        ax.set_title(title, fontsize=8)


def draw_theme_heatmap(ax, category_summary: pd.DataFrame, *,
                       significance: str = "significant",
                       top_n: int = 20,
                       drop_unmapped: bool = True,
                       comparisons: Optional[Sequence[str]] = None,
                       show_counts: bool = True,
                       row_order: str = "magnitude",
                       label_chars: int = 26,
                       fontsize: float = 6.5,
                       colorbar: bool = True,
                       colorbar_horizontal: bool = False,
                       title: Optional[str] = None) -> None:
    """Themes x comparisons from the category summarisation.

    Reads the summariser's own output, so every cell that has a significant term
    is filled.  A blank cell therefore means *no significant term in that
    comparison* — which matters, because a blank in a heatmap otherwise reads as
    an omission.
    """
    c = category_summary
    if "significance" in c.columns:
        c = c[c["significance"] == significance]
    if drop_unmapped:
        c = c[c["category"] != "unmapped"]
    if c.empty:
        ax.set_axis_off()
        return
    wide = c.pivot_table(index="category", columns="comparison",
                         values="mean_effect", aggfunc="first")
    counts = c.pivot_table(index="category", columns="comparison",
                           values="n_pathways", aggfunc="first")
    if comparisons:
        keep = [x for x in comparisons if x in wide.columns]
        wide, counts = wide[keep], counts[keep]
    # "magnitude" ranks by the strongest effect a theme reaches anywhere, which
    # is the right pick when the question is which themes matter.  It also
    # interleaves the two directions, so a reader scanning for "what goes up"
    # reads colour row by row; "signed" keeps the same selection but sorts it on
    # the value, which puts the two directions in blocks.
    top = wide.abs().max(axis=1).sort_values(ascending=False).index[:top_n]
    if row_order == "signed":
        key = wide.mean(axis=1).sort_values(ascending=False)
        order = [c for c in key.index if c in set(top)]
    elif row_order == "magnitude":
        order = list(top)
    else:
        raise ValueError("row_order must be 'magnitude' or 'signed'")
    wide, counts = wide.loc[order], counts.loc[order]

    lim = float(np.nanmax(np.abs(wide.to_numpy(dtype=float)))) or 1.0
    im = ax.imshow(wide.to_numpy(dtype=float), cmap="RdBu_r", vmin=-lim,
                   vmax=lim, aspect="auto")
    ax.set_xticks(range(wide.shape[1]))
    ax.set_xticklabels([str(x).replace("_vs_", "\nvs ") for x in wide.columns],
                       fontsize=fontsize)
    ax.set_yticks(range(wide.shape[0]))
    # theme names run long; wrapped rather than left to reach into whatever
    # panel sits to the left of this one
    ax.set_yticklabels(["\n".join(textwrap.wrap(str(t), label_chars))
                        for t in wide.index], fontsize=fontsize)
    if show_counts:
        for i in range(wide.shape[0]):
            for j in range(wide.shape[1]):
                n, v = counts.iat[i, j], wide.iat[i, j]
                if pd.isna(n):
                    continue
                dark = pd.notna(v) and abs(float(v)) > 0.62 * lim
                ax.text(j, i, f"{int(n)}", ha="center", va="center",
                        fontsize=fontsize - 1.1,
                        color="white" if dark else "0.15")
    if title:
        ax.set_title(title, fontsize=8)
    if colorbar:
        if colorbar_horizontal:
            cb = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.85,
                              pad=0.17, aspect=28)
        else:
            cb = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
        cb.set_label("mean effect", fontsize=7)
        cb.ax.tick_params(labelsize=6)


def draw_effect_concordance(ax, x: pd.Series, y: pd.Series, *,
                            xlabel: str = "x", ylabel: str = "y",
                            title: Optional[str] = None,
                            annotate_rho: bool = True,
                            point_size: float = 18.0,
                            unit: str = "shared",
                            poles: Optional[tuple] = None,
                            colour_by_agreement: bool = False,
                            annotate: int | Sequence[str] = 0,
                            annotate_disagreements: bool = False,
                            label_chars: int = 30) -> float:
    """Effect sizes from two analyses against each other, on a shared square.

    Use to ask whether a choice changed a conclusion: the same terms or themes
    scored under two contrasts, two collections, or two settings.  Returns the
    Spearman correlation so the caller can report it.  Axes are symmetric about
    zero with the identity line drawn, because the question is agreement in sign
    and magnitude, not fit.

    A scatter plot of two effect columns is only readable if the reader knows
    what one point is, which way is which, and which points are the ones to
    look at, so the optional arguments put those on the panel rather than in a
    caption:

    ``unit``
        what one point is -- ``"GO-slim theme"``, ``"pathway"``.  It is named
        beside the count in the subtitle.
    ``poles``
        ``(negative, positive)``, what the two ends of the *shared* scale mean
        -- ``("higher in cell lines", "higher in tumours")``.  Both axes carry
        the same quantity, so the direction is written onto both axis labels
        and the two quadrants in which the analyses agree are shaded.
    ``colour_by_agreement``
        colour the points that disagree in sign differently from those that
        agree, so a disagreement is visible rather than counted.
    ``annotate``
        how many of the strongest points to name, split between the two ends,
        or the index labels to name explicitly.  ``annotate_disagreements``
        additionally names every point whose sign differs between the two
        analyses -- in a panel about agreement those are the points a reader
        goes looking for, and they are otherwise anonymous dots.
    """
    j = pd.concat([pd.Series(x).rename("x"), pd.Series(y).rename("y")],
                  axis=1).dropna()
    if j.empty:
        ax.set_axis_off()
        return float("nan")
    rho = float(j["x"].corr(j["y"], method="spearman"))
    lim = float(np.abs(j.to_numpy(dtype=float)).max()) * 1.18 or 1.0
    agree = np.sign(j["x"]) == np.sign(j["y"])

    if poles is not None:
        # the two quadrants in which the analyses agree.  Shading them is what
        # makes "the sign never changes" a thing to see rather than to be told
        from matplotlib.patches import Rectangle
        for sx, sy in ((-1, -1), (1, 1)):
            ax.add_patch(Rectangle((0 if sx > 0 else -lim,
                                    0 if sy > 0 else -lim), lim, lim,
                                   facecolor="#eef2f7", edgecolor="none",
                                   zorder=0))
    ax.axhline(0, color="0.8", lw=0.7, zorder=1)
    ax.axvline(0, color="0.8", lw=0.7, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], ls="--", color="0.55", lw=0.9, zorder=2)
    # the identity line needs saying: without it the reader has no reason to
    # read distance from the diagonal as anything
    if poles is None:
        # only where nothing else explains the line.  With quadrant labels the
        # corners are taken and the points hug the diagonal, so this lands on
        # top of one or the other; there, the caller's caption says it instead
        ax.text(lim * 0.55, lim * 0.55, "equal effect", fontsize=6,
                color="0.45", rotation=45, rotation_mode="anchor",
                ha="left", va="top", zorder=3)

    colours = (np.where(agree, "#4575b4", "#d73027") if colour_by_agreement
               else "#4575b4")
    ax.scatter(j["x"], j["y"], s=point_size, c=colours, alpha=0.85,
               edgecolors="white", linewidths=0.3, zorder=4)

    # ── the points worth naming ──────────────────────────────────────────
    # An anonymous cloud invites the reader to take the correlation on trust.
    # Naming the extremes and every disagreement turns it into a panel that can
    # be read: those are the points an argument about agreement rests on.
    mid = (j["x"] + j["y"]) / 2
    if isinstance(annotate, int):
        up, down = (annotate + 1) // 2, annotate // 2
        chosen = list(mid.nlargest(up).index) + list(mid.nsmallest(down).index)
    else:
        chosen = [k for k in annotate if k in j.index]
    if annotate_disagreements:
        chosen += [k for k in j.index[~agree] if k not in chosen]

    # The points cluster where they are strongest, so a label at a fixed offset
    # from its point lands on a neighbour or outside the axes.  The two
    # off-diagonal triangles are empty by construction -- the cloud hugs the
    # diagonal -- so labels are stacked down the upper-left one and up the
    # lower-right one, each with a leader to its point, and the leaders cross
    # nothing.
    up, down = [], []
    for key in chosen:
        px, py = float(j.loc[key, "x"]), float(j.loc[key, "y"])
        (up if py > px else down).append((px, py, key))
    def _wrap(key):
        # wrapped to two lines rather than truncated: a theme cut off at
        # "extracellular matrix..." is not a name a reader can use
        lines = textwrap.wrap(str(key), label_chars) or [str(key)]
        if len(lines) > 2:
            lines = [lines[0], textwrap.shorten(" ".join(lines[1:]),
                                                label_chars, placeholder="…")]
        return "\n".join(lines), len(lines)

    for items, side in ((up, 1), (down, -1)):
        items.sort(key=lambda t: t[1], reverse=side > 0)
        x_text = -0.94 * lim if side > 0 else 0.94 * lim
        cursor = 0.94
        # the pitch has to clear a two-line label, not a one-line one: the
        # panel is square and its height is set by the narrower of the cell's
        # two sides, so a stack tuned on a wide panel collides on a tall one
        PITCH, LINE = 0.118, 0.046
        for px, py, key in items:
            text, n_lines = _wrap(key)
            cursor -= LINE * (n_lines - 1)           # centre the taller ones
            ax.annotate(
                text, xy=(px, py), xytext=(x_text, side * cursor * lim),
                textcoords="data", fontsize=5.2, color="0.25",
                linespacing=1.2, ha="left" if side > 0 else "right",
                va="center", zorder=6,
                arrowprops=dict(arrowstyle="-", lw=0.4, color="0.55",
                                shrinkA=1.0, shrinkB=2.0))
            cursor -= LINE * (n_lines - 1) + PITCH

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    if poles is not None:
        # both axes carry the same quantity, so an arrow cue on each of them
        # sets the same two phrases four times over one panel.  Said once
        # instead, in the two corners where the analyses agree -- which is
        # where the shading already is, so the label explains the shading
        # rather than repeating the axis.
        neg, pos = poles
        for txt, px, py, ha, va in ((pos, 0.975, 0.975, "right", "top"),
                                    (neg, 0.025, 0.025, "left", "bottom")):
            ax.text(px, py, txt, transform=ax.transAxes, ha=ha, va=va,
                    fontsize=6.4, color="0.40", style="italic", zorder=5)
    head = title or ""
    if annotate_rho:
        # "68 themes, Spearman 0.983, 3 disagree in sign" needed a footnote to
        # read: 68 of what, and disagree about what.  Said in words instead --
        # what is being counted, and what the exceptions actually do.
        n_bad = int((~agree).sum())
        one = unit[:-1] if unit.endswith("s") else unit
        # the count of shared points is a property of the collection, not of
        # this panel: down a row of panels it is the same number every time and
        # belongs in the legend.  What varies panel to panel is the agreement.
        sub = (f"Spearman ρ = {rho:.3f} · "
               + (f"all {len(j)} point the same way" if not n_bad else
                  f"{n_bad} {unit if n_bad > 1 else one} "
                  f"point{'' if n_bad > 1 else 's'} the opposite way"))
        head = f"{head}\n{sub}" if head else sub
    if head:
        ax.set_title(head, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    return rho


def render_effect_concordance(x, y, out_path: str | Path, **kw) -> Path:
    dpi = kw.pop("dpi", 150)
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    draw_effect_concordance(ax, x, y, **kw)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def render_theme_heatmap(category_summary: pd.DataFrame, out_path: str | Path,
                         **kw) -> Path:
    dpi = kw.pop("dpi", 150)
    top_n = kw.get("top_n", 20)
    n_comp = category_summary["comparison"].nunique() or 1
    fig, ax = plt.subplots(figsize=(2.2 + 1.35 * n_comp, 1.7 + 0.32 * top_n))
    draw_theme_heatmap(ax, category_summary, **kw)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return Path(out_path)
