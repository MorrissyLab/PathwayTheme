"""Figures for the differential pathway stage: one volcano + one top-pathway
bar panel per comparison."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..contracts import DiffResult
from .palettes import truncate
from .diagnostics import _short
from .pca_figures import _safe_savefig

_NEGLOG_CAP = 320.0   # -log10(p) cap so p==0 rows stay on-plot


def _neglog10(p: pd.Series) -> np.ndarray:
    p = p.to_numpy(dtype=float)
    with np.errstate(divide="ignore"):
        v = -np.log10(p)
    v[~np.isfinite(v)] = _NEGLOG_CAP
    return np.clip(v, 0.0, _NEGLOG_CAP)


def draw_volcano(ax, df: pd.DataFrame, *, comparison: str | None = None,
                 alpha: float = 0.05, n_labels: int = 3, label_chars: int = 24,
                 short_labels: bool = False, legend: bool = True,
                 fontsize: float = 8.0) -> None:
    """Volcano onto a supplied ``ax``.

    ``n_labels`` terms are labelled *per side*, staggered with leader lines, and
    the x-axis is widened to make room.  The strongest effects in a pathway
    contrast are usually near-duplicate terms at almost identical coordinates,
    so labelling many of them stacks text rather than adding information.
    """
    d = df.dropna(subset=["effect", "fdr"])
    if d.empty:
        ax.set_axis_off()
        return
    y = _neglog10(d["fdr"])
    eff = d["effect"].to_numpy(dtype=float)
    sig = d["fdr"].to_numpy(dtype=float) < alpha
    ax.scatter(eff[~sig], y[~sig], s=5, c="#bbbbbb", alpha=0.6, linewidths=0,
               label=f"FDR ≥ {alpha}")
    for sel, colour, lab in ((sig & (eff >= 0), "#d62728", "up"),
                             (sig & (eff < 0), "#1f77b4", "down")):
        ax.scatter(eff[sel], y[sel], s=6, c=colour, alpha=0.8, linewidths=0,
                   label=lab)
    ax.axhline(-np.log10(alpha), color="#888", lw=0.7, ls="--")
    ax.text(0.0, -np.log10(alpha), f" FDR {alpha}", fontsize=6, color="0.35",
            va="bottom", ha="left", transform=ax.get_yaxis_transform())

    if n_labels:
        lo, hi = ax.get_xlim()
        ax.set_xlim(lo - 0.42 * (hi - lo), hi + 0.42 * (hi - lo))
        dd = d.assign(_y=y)[sig]
        fmt = _short if short_labels else (lambda t: truncate(t, 60))
        for side, sub in (("+", dd.nlargest(n_labels, "effect")),
                          ("-", dd.nsmallest(n_labels, "effect"))):
            for k, (_, r) in enumerate(sub.iterrows()):
                ax.annotate(textwrap.fill(fmt(r["pathway"]), label_chars),
                            (r["effect"], r["_y"]), fontsize=4.6,
                            xytext=(6 if side == "+" else -6, 12 - 24 * k),
                            textcoords="offset points",
                            ha="left" if side == "+" else "right",
                            va="center", color="0.2",
                            arrowprops=dict(arrowstyle="-", lw=0.4,
                                            color="0.55", shrinkA=0, shrinkB=1))
    n_sig, n_up = int(sig.sum()), int((sig & (eff >= 0)).sum())
    ax.set_xlabel("effect (mean_case − mean_reference)", fontsize=fontsize)
    ax.set_ylabel("$-\\log_{10}$ FDR", fontsize=fontsize)
    if comparison:
        ax.set_title(f"{comparison}\n{n_sig} of {len(d)} significant "
                     f"({n_up} up, {n_sig - n_up} down)", fontsize=fontsize + 1)
    if legend:
        ax.legend(fontsize=6, frameon=False, loc="lower right", markerscale=1.8)
    ax.tick_params(labelsize=fontsize - 1)
    ax.spines[["top", "right"]].set_visible(False)


def draw_top_bars(ax, df: pd.DataFrame, *, comparison: str | None = None,
                  top_n: int = 25, per_side: int | None = None,
                  label_chars: int = 46, short_labels: bool = False,
                  wrap: bool = True, fontsize: float = 8.0) -> None:
    """Strongest effects as a horizontal bar panel, onto a supplied ``ax``.

    Default selects ``top_n`` by FDR.  ``per_side`` instead takes that many by
    largest positive and largest negative effect, which is the balanced view a
    figure panel usually wants.

    ``wrap=False`` shortens a long term to ``label_chars`` instead of folding
    it onto a second line.  In a composed panel the bars are as far apart as
    the panel is tall divided by their number, and a two-line label is taller
    than that gap, so wrapping is what makes neighbouring labels collide.
    """
    d = df.dropna(subset=["effect"]).copy()
    if d.empty:
        ax.set_axis_off()
        return
    if per_side:
        d = pd.concat([d.nlargest(per_side, "effect"),
                       d.nsmallest(per_side, "effect")])
    else:
        d = d.sort_values("fdr", na_position="last").head(top_n)
    d = d.sort_values("effect")
    fmt = _short if short_labels else (lambda t: truncate(t, 60))
    colors = ["#d62728" if e >= 0 else "#1f77b4" for e in d["effect"]]
    ax.barh(range(len(d)), d["effect"], color=colors, height=0.7)
    ax.set_yticks(range(len(d)))
    lab = ((lambda s: textwrap.fill(s, label_chars)) if wrap
           else (lambda s: truncate(s, label_chars)))
    ax.set_yticklabels([lab(fmt(p)) for p in d["pathway"]],
                       fontsize=max(5.0, fontsize - 2.0))
    ax.axvline(0, color="#333", lw=0.7)
    ax.set_xlabel("effect (mean_case − mean_reference)", fontsize=fontsize)
    if comparison:
        ax.set_title(comparison, fontsize=fontsize + 1)
    ax.tick_params(axis="x", labelsize=fontsize - 1)
    ax.spines[["top", "right"]].set_visible(False)


def _draw_volcano(df: pd.DataFrame, comparison: str, out_path: Path, *,
                  alpha: float = 0.05, dpi: int = 120) -> None:
    if df.dropna(subset=["effect", "fdr"]).empty:
        return
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    draw_volcano(ax, df, comparison=comparison, alpha=alpha)
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def _draw_top_bars(df: pd.DataFrame, comparison: str, out_path: Path, *,
                   top_n: int = 25, dpi: int = 120) -> None:
    d = df.dropna(subset=["effect"])
    if d.empty:
        return
    n = min(top_n, len(d))
    fig, ax = plt.subplots(figsize=(9.0, max(3.0, 0.32 * n + 1.0)))
    draw_top_bars(ax, df, comparison=f"{comparison} — top {n} by FDR",
                  top_n=top_n)
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def render_diff_figures(result: DiffResult, out_dir: str | Path, prefix: str, *,
                        top_n: int = 25, dpi: int = 120) -> list[Path]:
    """Volcano + top-pathway bars for each comparison in the DiffResult."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for comp in result.comparisons:
        sub = result.table[result.table["comparison"] == comp]
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in comp)
        vp = out_dir / f"{prefix}_diff_{safe}_volcano.pdf"
        bp = out_dir / f"{prefix}_diff_{safe}_top.pdf"
        _draw_volcano(sub, comp, vp, dpi=dpi)
        _draw_top_bars(sub, comp, bp, top_n=top_n, dpi=dpi)
        for p in (vp, bp):
            if p.exists():
                written.append(p)
    return written
