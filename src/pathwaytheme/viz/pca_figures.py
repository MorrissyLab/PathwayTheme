"""PCA figure generation — the 11 PDFs, ported from _pca_rollout.py.

Blocks (matching pca_analysis.ipynb):
  1. cluster x PC scores heatmap
  2. PC1-PC2 biplot + pairwise grid  (x2 colourings: primary + secondary label)
  3. per-observation pathway signatures
  4. per-PC loading bars + pathway x PC heatmap + clustered heatmap
  5. union of top +/- loadings per PC  (clustermap + signature view)

All functions consume a :class:`PCAResult`.  Colouring uses the result's
``labels`` (primary) and ``sublabels`` (secondary).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

from ..config import PCAConfig, VizConfig
from ..contracts import PCAResult
from .palettes import R_BWR, truncate, class_colour_lut, celltype_colour_lut

try:
    from adjustText import adjust_text as _adjust_text
except ImportError:
    _adjust_text = None


# ── save helper (retries on transient Windows PermissionError) ────────────
def _safe_savefig(fig, path, **kw) -> bool:
    for _ in range(8):
        try:
            fig.savefig(path, **kw)
            return True
        except PermissionError:
            time.sleep(2)
    print(f"  SAVE FAILED: {path}", flush=True)
    return False


def _short_label(lbl: str) -> str:
    """`cl3 (tumor)` -> `cl3`; otherwise the label unchanged."""
    s = str(lbl)
    return s.split(" (", 1)[0] if " (" in s else s


# ─────────────────────────────────────────────────────────────────────────
# Block 1 — cluster x PC scores heatmap
# ─────────────────────────────────────────────────────────────────────────
def _draw_scores_heatmap(res, title, out_path, dpi):
    scores_df = res.scores
    var_expl = res.variance_explained
    n_obs, K = scores_df.shape
    class_lut = class_colour_lut(sorted(set(res.labels)))
    col_colors = [class_lut.get(res.labels.iloc[i], "#dddddd") for i in range(n_obs)]

    body_w = 0.35 * K + 1.5
    body_h = 0.32 * n_obs + 1.5
    vmax = max(1.0, float(np.percentile(np.abs(scores_df.values), 98)))

    g = sns.clustermap(
        scores_df, cmap=R_BWR, vmin=-vmax, vmax=vmax, center=0,
        row_cluster=(n_obs > 1), col_cluster=False,
        row_colors=col_colors,
        xticklabels=[f"{p}\n({var_expl[i]*100:.0f}%)" for i, p in enumerate(res.pc_names)],
        yticklabels=scores_df.index,
        figsize=(max(7.0, body_w + 4.0), max(5.0, body_h + 2.0)),
        dendrogram_ratio=(0.10, 0.04),
        cbar_pos=(0.02, 0.83, 0.015, 0.12),
        cbar_kws={"label": "PC score"},
    )
    g.ax_heatmap.tick_params(axis="x", labelsize=8, rotation=0)
    g.ax_heatmap.tick_params(axis="y", labelsize=7, rotation=0)
    g.ax_heatmap.set_title(title, fontsize=10, pad=8, loc="left")
    handles = [Patch(facecolor=class_lut[c], label=c or "(none)") for c in sorted(class_lut)]
    g.ax_heatmap.legend(handles=handles, title=res.labels.name or "group",
                        bbox_to_anchor=(1.18, 1.0), loc="upper left",
                        fontsize=7, title_fontsize=8, frameon=False)
    _safe_savefig(g.fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(g.fig)


# ─────────────────────────────────────────────────────────────────────────
# Block 2 — biplot + pairwise grid
# ─────────────────────────────────────────────────────────────────────────
def _draw_biplot(scores_df, loadings_df, var_expl, point_labels, group_labels,
                 sizes, lut, legend_title, title_str, out_path, dpi):
    point_colors = [lut.get(g, "#dddddd") for g in group_labels]
    n = np.log1p(sizes)
    pt_sizes = 40 + 320 * (n / max(1.0, n.max()))

    fig, ax = plt.subplots(figsize=(12.0, 9.0))
    ax.axhline(0, color="#bbb", lw=0.6, zorder=0)
    ax.axvline(0, color="#bbb", lw=0.6, zorder=0)
    x, y = scores_df["PC1"].values, scores_df["PC2"].values
    ax.scatter(x, y, s=pt_sizes, c=point_colors, edgecolor="#222",
               linewidth=0.7, alpha=0.9, zorder=3)

    label_objs = [ax.text(xi, yi, pl, fontsize=7, color="#222",
                          fontweight="bold", zorder=5)
                  for pl, xi, yi in zip(point_labels, x, y)]

    lp = loadings_df.loc[["PC1", "PC2"]].T
    mag = (lp["PC1"]**2 + lp["PC2"]**2).pow(0.5)
    top10 = mag.sort_values(ascending=False).head(10).index
    span_x = (x.max() - x.min()) * 0.5 or 1.0
    span_y = (y.max() - y.min()) * 0.5 or 1.0
    arrow_scale = 0.7 * min(span_x, span_y) / lp.loc[top10].abs().values.max()

    arrow_label_objs = []
    for p in top10:
        ax_x = lp.loc[p, "PC1"] * arrow_scale
        ax_y = lp.loc[p, "PC2"] * arrow_scale
        ax.annotate("", xy=(ax_x, ax_y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.8, alpha=0.7),
                    zorder=2)
        arrow_label_objs.append(
            ax.text(ax_x * 1.08, ax_y * 1.08, truncate(p, 55), fontsize=6,
                    color="#444", style="italic", zorder=5,
                    ha="left" if ax_x >= 0 else "right",
                    va="bottom" if ax_y >= 0 else "top"))

    if _adjust_text is not None:
        _adjust_text(label_objs + arrow_label_objs, ax=ax,
                     arrowprops=dict(arrowstyle="-", color="#bbb", lw=0.4, alpha=0.7),
                     expand=(1.1, 1.2), time_limit=2.0, only_move={"text": "xy"})

    xmin, xmax = ax.get_xlim(); ymin, ymax = ax.get_ylim()
    ax.set_xlim(xmin - 0.06 * (xmax - xmin), xmax + 0.06 * (xmax - xmin))
    ax.set_ylim(ymin - 0.06 * (ymax - ymin), ymax + 0.06 * (ymax - ymin))
    ax.set_xlabel(f"PC1 ({var_expl[0]*100:.1f}% variance)", fontsize=10)
    ax.set_ylabel(f"PC2 ({var_expl[1]*100:.1f}% variance)", fontsize=10)
    ax.set_title(title_str, fontsize=11, loc="left")

    present = sorted({g for g in group_labels if g})
    handles = [Patch(facecolor=lut.get(c, "#dddddd"), label=c) for c in present]
    if any(not g for g in group_labels):
        handles.append(Patch(facecolor="#dddddd", label="(none)"))
    ax.legend(handles=handles, title=legend_title, bbox_to_anchor=(1.02, 1.0),
              loc="upper left", fontsize=8, title_fontsize=9, frameon=False)
    fig.tight_layout()
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def _draw_pairwise(scores_df, var_expl, pc_names, point_labels, group_labels,
                   sizes, lut, legend_title, title_str, out_path, dpi):
    K = scores_df.shape[1]
    point_colors = [lut.get(g, "#dddddd") for g in group_labels]
    n = np.log1p(sizes)
    pt_sizes = 30 + 200 * (n / max(1.0, n.max()))
    n_show = min(4, K)
    show_pcs = pc_names[:n_show]
    fig, axes = plt.subplots(n_show, n_show, figsize=(n_show * 3.6, n_show * 3.4),
                             squeeze=False)
    for i, pci in enumerate(show_pcs):
        for j, pcj in enumerate(show_pcs):
            ax = axes[i][j]
            if i == j:
                ax.text(0.5, 0.55, pci, ha="center", va="center", fontsize=18,
                        fontweight="bold", transform=ax.transAxes)
                ax.text(0.5, 0.30, f"{var_expl[i]*100:.1f}% var", ha="center",
                        va="center", fontsize=9, color="#555", transform=ax.transAxes)
                ax.set_xticks([]); ax.set_yticks([])
                for s in ax.spines.values():
                    s.set_visible(False)
                continue
            xi, yi = scores_df[pcj].values, scores_df[pci].values
            ax.axhline(0, color="#ddd", lw=0.5, zorder=0)
            ax.axvline(0, color="#ddd", lw=0.5, zorder=0)
            ax.scatter(xi, yi, s=pt_sizes * 0.6, c=point_colors, edgecolor="#333",
                       linewidth=0.4, alpha=0.85, zorder=3)
            for pl, x_v, y_v in zip(point_labels, xi, yi):
                ax.text(x_v, y_v, pl, fontsize=4.5, color="#333", zorder=4,
                        ha="left", va="bottom")
            ax.tick_params(axis="both", labelsize=6)
            if i == n_show - 1:
                ax.set_xlabel(pcj, fontsize=8)
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel(pci, fontsize=8)
            else:
                ax.set_yticklabels([])
    fig.suptitle(title_str, fontsize=11, y=0.995)
    present = sorted({g for g in group_labels if g})
    handles = [Patch(facecolor=lut.get(c, "#dddddd"), label=c) for c in present]
    if any(not g for g in group_labels):
        handles.append(Patch(facecolor="#dddddd", label="(none)"))
    fig.legend(handles=handles, title=legend_title, bbox_to_anchor=(1.005, 0.97),
               loc="upper left", fontsize=8, title_fontsize=9, frameon=False)
    fig.tight_layout(rect=[0, 0, 0.92, 0.97])
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
# Block 3 — per-observation signature panels
# ─────────────────────────────────────────────────────────────────────────
def _per_obs_from_signatures(sig_df: pd.DataFrame):
    """Rebuild {obs: (top_pcs, signed_series)} from the long signature table."""
    per_obs: dict = {}
    for obs, grp in sig_df.groupby("cluster", sort=False):
        top_pcs = grp["top_pcs"].iloc[0].split(",")
        ser = pd.Series(grp["signature_score"].values, index=grp["pathway"].values)
        per_obs[obs] = (top_pcs, ser)
    return per_obs


def _draw_signature_panels(res, title, out_path, dpi):
    per_obs = _per_obs_from_signatures(res.signatures)
    n_obs = len(per_obs)
    n_cols = 4
    n_rows = int(np.ceil(n_obs / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5.0, n_rows * 4.2),
                             squeeze=False)
    for idx, (c_label, (top_pcs, top_paths)) in enumerate(per_obs.items()):
        ax = axes[idx // n_cols][idx % n_cols]
        tp = top_paths.sort_values(ascending=True)
        colors = ["#d62728" if v > 0 else "#1f77b4" for v in tp.values]
        ypos = np.arange(len(tp))
        ax.barh(ypos, tp.values, color=colors, edgecolor="none")
        ax.set_yticks(ypos)
        ax.set_yticklabels([truncate(p, 70) for p in tp.index], fontsize=6)
        ax.tick_params(axis="x", labelsize=7)
        ax.axvline(0, color="#333", lw=0.5)
        ax.set_title(f"{c_label}\nsignature PCs: {', '.join(top_pcs)}",
                     fontsize=8, loc="left")
        ax.set_xlabel("signature score (red=up, blue=down)", fontsize=7)
    for j in range(n_obs, n_rows * n_cols):
        axes[j // n_cols][j % n_cols].set_visible(False)
    fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    _safe_savefig(fig, out_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
# Block 4 — per-PC loadings bars + pathway x PC heatmaps
# ─────────────────────────────────────────────────────────────────────────
def _draw_pc_loadings(res, tag_title, out_dir, prefix, dpi, top_n=25):
    loadings_df, var_expl, pc_names = res.loadings, res.variance_explained, res.pc_names
    K = loadings_df.shape[0]

    n_cols = min(3, K)
    n_rows = int(np.ceil(K / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6.0, n_rows * 4.5),
                             squeeze=False)
    for idx, pc in enumerate(pc_names):
        ax = axes[idx // n_cols][idx % n_cols]
        lv = loadings_df.loc[pc]
        pos = lv.nlargest(top_n)
        neg = lv.nsmallest(top_n).iloc[::-1]
        combined = pd.concat([neg, pos]).sort_values(ascending=True)
        colors = ["#d62728" if v > 0 else "#1f77b4" for v in combined.values]
        ypos = np.arange(len(combined))
        ax.barh(ypos, combined.values, color=colors, edgecolor="none")
        ax.set_yticks(ypos)
        ax.set_yticklabels([truncate(p, 65) for p in combined.index], fontsize=5.5)
        ax.tick_params(axis="x", labelsize=7)
        ax.axvline(0, color="#333", lw=0.5)
        ax.set_title(f"{pc} ({var_expl[idx]*100:.1f}% variance)", fontsize=9, loc="left")
        ax.set_xlabel("loading (red = positive axis end, blue = negative)", fontsize=7)
    for j in range(K, n_rows * n_cols):
        axes[j // n_cols][j % n_cols].set_visible(False)
    fig.suptitle(f"{tag_title} — top pathway loadings per principal component",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    _safe_savefig(fig, out_dir / f"{prefix}_pca_pc_loadings.pdf", bbox_inches="tight", dpi=dpi)
    plt.close(fig)

    # pathway x PC heatmap (top 80 by max |loading|)
    max_abs = loadings_df.abs().max(axis=0)
    top_paths = max_abs.sort_values(ascending=False).head(80).index
    sub = loadings_df.loc[:, top_paths].T
    vmax = float(np.percentile(np.abs(sub.values), 98))
    fig, ax = plt.subplots(figsize=(max(7.0, 0.6 * K + 2.0),
                                    max(8.0, 0.18 * len(top_paths) + 1.5)))
    im = ax.imshow(sub.values, aspect="auto", cmap=R_BWR, vmin=-vmax, vmax=vmax,
                   interpolation="nearest")
    im.set_rasterized(True)
    ax.set_yticks(np.arange(len(top_paths)))
    ax.set_yticklabels([truncate(p, 70) for p in top_paths], fontsize=6)
    ax.set_xticks(np.arange(K))
    ax.set_xticklabels([f"{p}\n({var_expl[i]*100:.0f}%)" for i, p in enumerate(pc_names)],
                       fontsize=7)
    ax.tick_params(axis="x", which="both", length=0)
    ax.set_title(f"{tag_title}\nTop {len(top_paths)} pathways x top {K} PCs "
                 f"- loadings heatmap (ordered by max |loading|)", fontsize=10, loc="left", pad=6)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("loading", fontsize=8); cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    _safe_savefig(fig, out_dir / f"{prefix}_pca_pathway_pc_heatmap.pdf",
                  bbox_inches="tight", dpi=dpi)
    plt.close(fig)

    # clustered pathway x PC heatmap
    cg = sns.clustermap(
        sub, cmap=R_BWR, vmin=-vmax, vmax=vmax, center=0,
        row_cluster=True, col_cluster=False, method="average", metric="euclidean",
        xticklabels=[f"{p} ({var_expl[i]*100:.0f}%)" for i, p in enumerate(pc_names)],
        yticklabels=[truncate(p, 70) for p in sub.index],
        figsize=(max(8.0, 0.7 * K + 3.0), max(9.0, 0.18 * len(top_paths) + 2.0)),
        dendrogram_ratio=(0.12, 0.04), cbar_pos=(0.02, 0.82, 0.015, 0.12),
        cbar_kws={"label": "loading"},
    )
    cg.ax_heatmap.tick_params(axis="x", labelsize=7, rotation=0)
    cg.ax_heatmap.tick_params(axis="y", labelsize=6, rotation=0)
    for c in cg.ax_heatmap.collections:
        c.set_rasterized(True)
    cg.ax_heatmap.set_title(f"{tag_title}\nTop {len(top_paths)} pathways x top {K} PCs "
                            f"- loadings clustered by pathway", fontsize=10, loc="left", pad=8)
    _safe_savefig(cg.fig, out_dir / f"{prefix}_pca_pathway_pc_clustermap.pdf",
                  bbox_inches="tight", dpi=dpi)
    plt.close(cg.fig)


# ─────────────────────────────────────────────────────────────────────────
# Block 5 — union of top +/- loadings per PC
# ─────────────────────────────────────────────────────────────────────────
def _draw_pc_loading_union(res, tag_title, out_dir, prefix, dpi, top_per_pc_per_sign=10):
    loadings_df, var_expl, pc_names = res.loadings, res.variance_explained, res.pc_names
    K = loadings_df.shape[0]

    union_paths = []
    for pc in pc_names:
        lv = loadings_df.loc[pc]
        union_paths.extend(lv.nlargest(top_per_pc_per_sign).index)
        union_paths.extend(lv.nsmallest(top_per_pc_per_sign).index)
    union_paths = list(dict.fromkeys(union_paths))
    sub_u = loadings_df.loc[:, union_paths].T
    vmax_u = float(np.percentile(np.abs(sub_u.values), 98))

    x_labels = [f"{p}\n{var_expl[i]*100:.0f}% var" for i, p in enumerate(pc_names)]
    fig_w = max(10.0, 0.8 * K + 5.0)
    fig_h = max(11.0, 0.18 * len(union_paths) + 3.0)

    # Plot A — clustered union
    cg = sns.clustermap(
        sub_u, cmap=R_BWR, vmin=-vmax_u, vmax=vmax_u, center=0,
        row_cluster=True, col_cluster=False, method="average", metric="euclidean",
        xticklabels=x_labels, yticklabels=[truncate(p, 55) for p in sub_u.index],
        figsize=(fig_w, fig_h), dendrogram_ratio=(0.12, 0.04),
        cbar_pos=(0.02, 0.84, 0.018, 0.10), cbar_kws={"label": "loading"},
    )
    for tick in cg.ax_heatmap.get_xticklabels():
        tick.set_rotation(35); tick.set_ha("right"); tick.set_va("top"); tick.set_fontsize(8)
    for tick in cg.ax_heatmap.get_yticklabels():
        tick.set_rotation(0); tick.set_fontsize(6)
    for c in cg.ax_heatmap.collections:
        c.set_rasterized(True)
    cg.fig.suptitle(f"{tag_title} - union of top +/-{top_per_pc_per_sign} loadings per PC "
                    f"({len(union_paths)} pathways) - clustered by pathway",
                    fontsize=11, y=0.995, x=0.5, ha="center")
    _safe_savefig(cg.fig, out_dir / f"{prefix}_pca_pc_loading_union_clustermap.pdf",
                  bbox_inches="tight", dpi=dpi)
    plt.close(cg.fig)

    # Plot B — signature view (rows grouped by home PC + sign)
    home_idx = sub_u.abs().values.argmax(axis=1)
    signed_home = np.array([sub_u.iloc[r, home_idx[r]] for r in range(len(sub_u))])
    order = pd.DataFrame({"pathway": sub_u.index, "home_idx": home_idx,
                          "signed_home": signed_home})
    order["sign_group"] = (order["signed_home"] < 0).astype(int)
    order["rank_in_block"] = -order["signed_home"].abs()
    order = order.sort_values(["home_idx", "sign_group", "rank_in_block"])
    sub_ord = sub_u.loc[order["pathway"].tolist()]

    tab10 = plt.get_cmap("tab10").colors
    pc_color = {p: tab10[i % 10] for i, p in enumerate(pc_names)}
    row_colors = [pc_color[pc_names[i]] for i in order["home_idx"]]

    cg = sns.clustermap(
        sub_ord, cmap=R_BWR, vmin=-vmax_u, vmax=vmax_u, center=0,
        row_cluster=False, col_cluster=False, row_colors=row_colors,
        xticklabels=x_labels, yticklabels=[truncate(p, 55) for p in sub_ord.index],
        figsize=(fig_w, fig_h), dendrogram_ratio=(0.02, 0.04),
        cbar_pos=(0.02, 0.84, 0.018, 0.10), cbar_kws={"label": "loading"},
    )
    for tick in cg.ax_heatmap.get_xticklabels():
        tick.set_rotation(35); tick.set_ha("right"); tick.set_va("top"); tick.set_fontsize(8)
    for tick in cg.ax_heatmap.get_yticklabels():
        tick.set_rotation(0); tick.set_fontsize(6)
    for c in cg.ax_heatmap.collections:
        c.set_rasterized(True)
    handles = [Patch(facecolor=pc_color[p], label=p) for p in pc_names]
    cg.fig.legend(handles=handles, title="top PC (max |loading|)", loc="upper left",
                  bbox_to_anchor=(0.005, 0.70), fontsize=7, title_fontsize=8,
                  frameon=False, ncol=1, labelspacing=0.3)
    cg.fig.suptitle(f"{tag_title} - union of top +/-{top_per_pc_per_sign} loadings per PC "
                    f"({len(union_paths)} pathways) - grouped by top PC + sign",
                    fontsize=11, y=0.995, x=0.5, ha="center")
    _safe_savefig(cg.fig, out_dir / f"{prefix}_pca_pc_loading_union_signature.pdf",
                  bbox_inches="tight", dpi=dpi)
    plt.close(cg.fig)


# ─────────────────────────────────────────────────────────────────────────
# Public entry point
# ─────────────────────────────────────────────────────────────────────────
def render_pca_figures(res: PCAResult, out_dir: str | Path, prefix: str,
                       title: str, viz: VizConfig | None = None) -> list[Path]:
    """Render all 11 PCA PDFs for one PCAResult.  Returns the written paths."""
    viz = viz or VizConfig()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dpi = viz.dpi

    scores_df, loadings_df = res.scores, res.loadings
    var_expl, pc_names = res.variance_explained, res.pc_names
    group_labels = list(res.labels.values)
    sub_labels = list(res.sublabels.values) if len(res.sublabels) else group_labels
    point_labels = [_short_label(l) for l in scores_df.index]
    sizes = res.sizes if res.sizes is not None else np.ones(scores_df.shape[0])

    class_lut = class_colour_lut(sorted(set(group_labels)))
    ct_lut = celltype_colour_lut(set(sub_labels))
    primary_name = res.labels.name or "group"
    secondary_name = res.sublabels.name or "secondary"

    # Block 1
    _draw_scores_heatmap(res, f"{title}\nPCA cluster scores (top {res.K} PCs)",
                         out_dir / f"{prefix}_pca_cluster_scores_heatmap.pdf", dpi)
    # Block 2 — biplots + pairwise, primary colouring
    _draw_biplot(scores_df, loadings_df, var_expl, point_labels, group_labels, sizes,
                 class_lut, primary_name, f"{title}\nPCA biplot (colour = {primary_name})",
                 out_dir / f"{prefix}_pca_scatter_pc1_pc2_by_class.pdf", dpi)
    _draw_pairwise(scores_df, var_expl, pc_names, point_labels, group_labels, sizes,
                   class_lut, primary_name, f"{title} - pairwise PCA scatter (colour = {primary_name})",
                   out_dir / f"{prefix}_pca_scatter_pairwise_by_class.pdf", dpi)
    # Block 2 — secondary colouring (only meaningful if sublabels differ)
    _draw_biplot(scores_df, loadings_df, var_expl, point_labels, sub_labels, sizes,
                 ct_lut, secondary_name, f"{title}\nPCA biplot (colour = {secondary_name})",
                 out_dir / f"{prefix}_pca_scatter_pc1_pc2_by_celltype.pdf", dpi)
    _draw_pairwise(scores_df, var_expl, pc_names, point_labels, sub_labels, sizes,
                   ct_lut, secondary_name, f"{title} - pairwise PCA scatter (colour = {secondary_name})",
                   out_dir / f"{prefix}_pca_scatter_pairwise_by_celltype.pdf", dpi)
    # Block 3
    _draw_signature_panels(res, f"{title} - per-cluster pathway signatures (PCA)",
                           out_dir / f"{prefix}_pca_per_cluster_signatures.pdf", dpi)
    # Block 4 (+ block 5 called within, matching the rollout structure)
    _draw_pc_loadings(res, title, out_dir, prefix, dpi)
    # Block 5
    _draw_pc_loading_union(res, title, out_dir, prefix, dpi)

    return sorted(out_dir.glob(f"{prefix}_pca_*.pdf"))
