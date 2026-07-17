"""Full-matrix sanity heatmap — a QC view of the whole score matrix.

Z-scores each pathway across samples (ddof=0) and draws a clustered heatmap,
optionally with a colour strip for one metadata label.  Mirrors the "all GO:BP
pathways" sanity heatmap in the reference ASPS workflow.  Very large matrices
are capped to the most-variable pathways (logged) so the figure stays legible.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

from ..contracts import ScoreMatrix
from .palettes import R_BWR, class_colour_lut, truncate
from .pca_figures import _safe_savefig


def render_sanity_heatmap(sm: ScoreMatrix, out_path: str | Path, *,
                          label_col: Optional[str] = None,
                          max_pathways: int = 200, dpi: int = 120,
                          title: str = "Pathway score matrix (z-scored)") -> Optional[Path]:
    """Write a z-scored pathway x sample clustered heatmap.  Returns the path."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    X = sm.data.fillna(0.0)
    sds = X.std(axis=1, ddof=0)
    X = X.loc[sds > 1e-12]
    if X.shape[0] < 2 or X.shape[1] < 2:
        return None

    if max_pathways and X.shape[0] > max_pathways:
        keep = X.std(axis=1, ddof=0).sort_values(ascending=False).head(max_pathways).index
        X = X.loc[keep]
        title = f"{title} — top {max_pathways} most-variable of {sds[sds > 1e-12].size}"

    Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1, ddof=0), axis=0)
    Z.index = [truncate(t, 60) for t in Z.index]

    col_colors = None
    lut = None
    if label_col and sm.metadata.has(label_col):
        labs = sm.metadata.get(label_col).reindex(sm.samples).astype(str)
        lut = class_colour_lut(sorted(set(labs)))
        col_colors = [lut.get(v, "#dddddd") for v in labs]

    h = min(40.0, max(4.0, 0.16 * Z.shape[0] + 2.0))
    w = min(30.0, max(6.0, 0.18 * Z.shape[1] + 3.0))
    g = sns.clustermap(
        Z, cmap=R_BWR, center=0, vmin=-2, vmax=2,
        col_colors=col_colors,
        row_cluster=True, col_cluster=True,
        xticklabels=Z.columns, yticklabels=(Z.shape[0] <= 120),
        figsize=(w, h), dendrogram_ratio=(0.06, 0.08),
        cbar_kws={"label": "z-score"},
    )
    g.ax_heatmap.tick_params(axis="x", labelsize=6, rotation=90)
    g.ax_heatmap.tick_params(axis="y", labelsize=5)
    g.ax_heatmap.set_title(title, fontsize=10, pad=8, loc="left")
    if lut is not None:
        handles = [Patch(facecolor=lut[c], label=c or "(none)") for c in sorted(lut)]
        g.ax_heatmap.legend(handles=handles, title=label_col,
                            bbox_to_anchor=(1.15, 1.0), loc="upper left",
                            fontsize=7, title_fontsize=8, frameon=False)
    _safe_savefig(g.fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(g.fig)
    return out_path
