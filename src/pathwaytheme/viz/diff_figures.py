"""Figures for the differential pathway stage: one volcano + one top-pathway
bar panel per comparison."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..contracts import DiffResult
from .palettes import truncate
from .pca_figures import _safe_savefig

_NEGLOG_CAP = 320.0   # -log10(p) cap so p==0 rows stay on-plot


def _neglog10(p: pd.Series) -> np.ndarray:
    p = p.to_numpy(dtype=float)
    with np.errstate(divide="ignore"):
        v = -np.log10(p)
    v[~np.isfinite(v)] = _NEGLOG_CAP
    return np.clip(v, 0.0, _NEGLOG_CAP)


def _draw_volcano(df: pd.DataFrame, comparison: str, out_path: Path, *,
                  alpha: float = 0.05, dpi: int = 120) -> None:
    d = df.dropna(subset=["effect", "fdr"])
    if d.empty:
        return
    y = _neglog10(d["fdr"])
    sig = d["fdr"].to_numpy() < alpha
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    ax.scatter(d["effect"][~sig], y[~sig], s=10, c="#bbbbbb", alpha=0.6,
               linewidths=0, label=f"FDR ≥ {alpha}")
    up = sig & (d["effect"].to_numpy() >= 0)
    dn = sig & (d["effect"].to_numpy() < 0)
    ax.scatter(d["effect"][up], y[up], s=14, c="#d62728", alpha=0.85,
               linewidths=0, label="up")
    ax.scatter(d["effect"][dn], y[dn], s=14, c="#1f77b4", alpha=0.85,
               linewidths=0, label="down")
    ax.axhline(-np.log10(alpha), color="#888", lw=0.7, ls="--")
    ax.axvline(0, color="#bbb", lw=0.6)
    # label the strongest few
    for _, r in d[d["fdr"] < alpha].sort_values("fdr").head(8).iterrows():
        ax.annotate(truncate(r["pathway"], 40), (r["effect"], -np.log10(max(r["fdr"], 1e-320))),
                    fontsize=6, alpha=0.9)
    ax.set_xlabel("effect (mean_case − mean_reference)")
    ax.set_ylabel("−log10(FDR)")
    ax.set_title(comparison, fontsize=10, loc="left")
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def _draw_top_bars(df: pd.DataFrame, comparison: str, out_path: Path, *,
                   top_n: int = 25, dpi: int = 120) -> None:
    d = df.dropna(subset=["effect"]).copy()
    if d.empty:
        return
    d = d.sort_values("fdr", na_position="last").head(top_n)
    d = d.sort_values("effect")
    colors = ["#d62728" if e >= 0 else "#1f77b4" for e in d["effect"]]
    fig, ax = plt.subplots(figsize=(9.0, max(3.0, 0.32 * len(d) + 1.0)))
    ax.barh([truncate(p, 60) for p in d["pathway"]], d["effect"], color=colors)
    ax.axvline(0, color="#333", lw=0.7)
    ax.set_xlabel("effect (mean_case − mean_reference)")
    ax.set_title(f"{comparison} — top {len(d)} by FDR", fontsize=10, loc="left")
    ax.tick_params(axis="y", labelsize=6)
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
