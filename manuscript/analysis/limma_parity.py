"""Fig. S12 and Table S16 -- the differential test against limma itself.

Python has no limma, so ``pathwaytheme.diff.stats`` reimplements limma's
empirical-Bayes moderated t (Smyth 2004) in pure Python: ``fit_fdist`` and
``trigamma_inverse`` are ports of ``limma::fitFDist`` and
``limma::trigammaInverse``.  Every differential effect in this manuscript comes
out of that port, so the port is checked against the thing it ports.

It is one port shared by every run, so what is being checked is the
implementation and not a particular contrast.  Four runs are compared -- both
collections, and both contrast forms (one-vs-rest and pairwise).

What is compared is not a re-run.  The Python side is read from the tables the
four runs below wrote -- the same files Tables S5 and S6 are built from --
and limma is run on the same score matrix, the same samples and the same model
(``~ 0 + group``, contrast case - reference).  The sample counts are asserted
against the ones the pipeline recorded before limma is run at all, so the two
sides cannot be compared on different contrasts.

Two levels are reported, because an error in the prior could cancel downstream:

  * the fitted hyperparameters d0 and s0^2, with ``fit_fdist`` given limma's
    own residual variances so that only the estimator differs;
  * the per-pathway effect, moderated t and p, over every pathway of every
    contrast.

R with the limma package is required; the ``Rscript`` executable is taken from
the environment variable ``RSCRIPT``, else from ``PATH``.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pathwaytheme.diff.stats import fit_fdist

from .paths import PIPE, PROC, SCORES, SFIG, STAB

R_DIR = Path(__file__).resolve().parent / "R"

# (run directory, gene set, the column the run grouped on, the run's row filter)
RUNS = [
    ("hallmark_sampletype", "HALLMARK", "sample_type", None),
    ("gobp_sampletype", "GOBP", "sample_type", None),
    ("gobp_tumor_normal", "GOBP", "sample_type", None),
    ("gobp_diagnosis", "GOBP", "diagnosis", ("sample_type", ["TUMOR"])),
]
INK, ACCENT, GRID = "#22262b", "#c0392b", "#c8ccd2"
EPS = float(np.finfo(float).eps)


def rscript() -> str:
    """The Rscript executable: ``$RSCRIPT``, else the one on ``PATH``."""
    exe = shutil.which(os.environ.get("RSCRIPT") or "Rscript")
    if not exe:
        raise RuntimeError(
            "Rscript not found. Install R with the limma package and put "
            "Rscript on PATH, or set the environment variable RSCRIPT to the "
            "full path of the Rscript executable.")
    return exe


def limma_version() -> str:
    """The limma and R versions the comparison runs against."""
    proc = subprocess.run(
        [rscript(), "-e", 'suppressPackageStartupMessages(library(limma));'
                          'cat(as.character(packageVersion("limma")),'
                          'R.version.string)'], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"limma is not installed:\n{proc.stderr[-800:]}")
    return proc.stdout.strip()


def _contrasts(run: str, geneset: str) -> pd.DataFrame:
    """The differential table the run wrote, as it ships."""
    path = PIPE / run / geneset / f"{geneset}_differential_moderated_t.tsv"
    df = pd.read_csv(path, sep="\t")
    if set(df["method"]) != {"moderated_t"}:
        raise RuntimeError(f"{path.name} is not a moderated-t table")
    return df


def _samples(run_filter) -> pd.DataFrame:
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)
    if run_filter is not None:
        col, keep = run_filter
        meta = meta[meta[col].isin(keep)]
    return meta


def _run_limma(scores: pd.DataFrame, groups: pd.Series, case: str,
               ref: str) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        spath, gpath, out = tmp / "s.tsv", tmp / "g.tsv", tmp / "o.tsv"
        scores.to_csv(spath, sep="\t")
        pd.DataFrame({"sample": groups.index,
                      "group": groups.to_numpy()}).to_csv(gpath, sep="\t",
                                                          index=False)
        proc = subprocess.run(
            [rscript(), str(R_DIR / "limma_reference.R"), str(spath),
             str(gpath), case, ref, str(out)], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"limma failed for {case} vs {ref}:\n"
                               f"{proc.stderr[-1500:]}")
        return pd.read_csv(out, sep="\t").set_index("pathway")


def compare_with_limma() -> pd.DataFrame:
    """Run limma on every contrast of the four runs and compare, per contrast.

    The returned frame is Table S16 plus a ``_t`` column holding the two
    per-pathway moderated-t vectors, which Fig. S12 plots.
    """
    print(f"limma {limma_version()}")
    rows = []
    for run, geneset, target_col, run_filter in RUNS:
        meta = _samples(run_filter)
        sm = pd.read_csv(SCORES / f"{geneset}_ssgsea_scores.tsv", sep="\t",
                         index_col=0)
        cols = [c for c in sm.columns if c in meta.index]
        sm, labels = sm[cols], meta.loc[cols, target_col].astype(str)
        table = _contrasts(run, geneset)
        for comparison, got in table.groupby("comparison", sort=True):
            case, _, ref = str(comparison).partition("_vs_")
            got = got.set_index("pathway")
            in_case = labels == case
            in_ref = (~in_case) if ref == "rest" else (labels == ref)
            # the pipeline recorded the sizes it ran at; if the groups rebuilt
            # here are not those groups, the comparison below is meaningless
            n_case, n_ref = int(in_case.sum()), int(in_ref.sum())
            want = (int(got["n_case"].iloc[0]), int(got["n_reference"].iloc[0]))
            if (n_case, n_ref) != want:
                raise RuntimeError(f"{run} {comparison}: rebuilt groups are "
                                   f"{n_case}+{n_ref}, the run recorded {want}")
            keep = labels.index[in_case | in_ref]
            g = pd.Series(np.where(labels[keep] == case, case, "reference"),
                          index=keep)
            lim = _run_limma(sm[keep], g, case, "reference")
            lim = lim.reindex(got.index)
            if lim["effect"].isna().any():
                raise RuntimeError(
                    f"{run} {comparison}: limma returned no row for "
                    f"{int(lim['effect'].isna().sum())} pathways")

            d0_py, s02_py = fit_fdist(lim["s2_ordinary"].to_numpy(float),
                                      lim["df_residual"].to_numpy(float))
            d0_r = float(lim["df_prior"].iloc[0])
            s02_r = float(lim["s2_prior"].iloc[0])
            e_py, e_r = got["effect"].to_numpy(float), lim["effect"].to_numpy(float)
            t_py, t_r = got["statistic"].to_numpy(float), lim["t_moderated"].to_numpy(float)
            p_py, p_r = got["p_value"].to_numpy(float), lim["p_moderated"].to_numpy(float)
            rows.append({
                "run": run, "geneset": geneset, "contrast": str(comparison),
                "form": "one-vs-rest" if ref == "rest" else "pairwise",
                "n_case": n_case, "n_reference": n_ref,
                "n_pathways": len(got),
                "d0_limma": d0_r, "d0_pathwaytheme": d0_py,
                "d0_rel_err": abs(d0_py - d0_r) / max(abs(d0_r), 1e-12),
                "s0sq_limma": s02_r, "s0sq_pathwaytheme": s02_py,
                "s0sq_rel_err": abs(s02_py - s02_r) / max(abs(s02_r), 1e-12),
                "effect_max_abs_diff": float(np.max(np.abs(e_py - e_r))),
                "t_max_abs_diff": float(np.max(np.abs(t_py - t_r))),
                "t_pearson": float(np.corrcoef(t_py, t_r)[0, 1]),
                "p_max_abs_diff": float(np.max(np.abs(p_py - p_r))),
                "n_disagreeing_calls": int(np.sum((p_py < 0.05) != (p_r < 0.05))),
                "_t": (t_py, t_r),
            })
            print(f"  {run:22s} {comparison:24s} {n_case:4d}+{n_ref:4d}  "
                  f"max|dt| = {rows[-1]['t_max_abs_diff']:.2e}")
    return pd.DataFrame(rows)


def figure_s12(parity: pd.DataFrame | None = None) -> plt.Figure:
    """Fig. S12: PathwayTheme's moderated t against limma::eBayes.

    ``parity`` is the output of :func:`compare_with_limma`; it is computed if
    not given.
    """
    res = compare_with_limma() if parity is None else parity
    t_py = np.concatenate([r[0] for r in res["_t"]])
    t_r = np.concatenate([r[1] for r in res["_t"]])
    # the right panel's row labels are long and hang to its left, so the gap
    # between the two panels is where they live
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.2),
                             gridspec_kw={"width_ratios": [1.0, 1.15],
                                          "wspace": 0.95})

    ax = axes[0]
    lim = float(np.max(np.abs(np.r_[t_py, t_r]))) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color=GRID, lw=1.0, zorder=1)
    ax.scatter(t_r, t_py, s=1.4, color=ACCENT, alpha=0.35, lw=0, zorder=2)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.set_xlabel("moderated t, limma::eBayes (R)", fontsize=8)
    ax.set_ylabel("moderated t, PathwayTheme (Python)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title("a  Every pathway of every contrast",
                 fontsize=8.5, loc="left", pad=8)

    ax = axes[1]
    order = res.sort_values(["geneset", "form", "contrast"]).reset_index(drop=True)
    y = np.arange(len(order))
    series = [("effect_max_abs_diff", "effect", "#4575b4", "o"),
              ("t_max_abs_diff", "moderated t", ACCENT, "s"),
              ("p_max_abs_diff", "p value", "#1a9850", "^")]
    floor = EPS / 4                       # so an exact zero still has a place
    for col, name, colour, marker in series:
        v = np.maximum(order[col].to_numpy(float), floor)
        ax.scatter(v, y, s=13, color=colour, marker=marker, label=name,
                   lw=0, zorder=3)
    ax.axvline(EPS, color=INK, lw=0.8, ls=(0, (3, 2)), zorder=2,
               label="double precision")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{'GO:BP' if r.geneset == 'GOBP' else 'hallmark'}   "
                        f"{r.contrast.replace('_vs_', ' vs ')}  "
                        f"({r.n_case}+{r.n_reference})"
                        for r in order.itertuples()], fontsize=6.0)
    ax.set_ylim(len(order) - 0.5, -0.8)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("largest absolute difference, over that contrast's pathways",
                  fontsize=7.5)
    ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(fontsize=6.8, frameon=False, ncol=4, handletextpad=0.2,
              columnspacing=1.2, loc="lower left", bbox_to_anchor=(0.0, 1.005))
    ax.set_title("b  Per contrast, the worst pathway", fontsize=8.5, loc="left",
                 pad=18)
    for sp in ("top", "right"):
        axes[0].spines[sp].set_visible(False)
        ax.spines[sp].set_visible(False)

    fig.suptitle("Fig. S12  PathwayTheme against limma::eBayes",
                 fontsize=10, y=1.005)
    SFIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(SFIG / f"figS12_limma_parity.{ext}", bbox_inches="tight",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)
    return fig


def table_s16(parity: pd.DataFrame | None = None) -> pd.DataFrame:
    """Table S16: PathwayTheme against limma, one row per contrast.

    ``parity`` is the output of :func:`compare_with_limma`; it is computed if
    not given.
    """
    res = compare_with_limma() if parity is None else parity
    out = res.drop(columns=["_t"])
    STAB.mkdir(parents=True, exist_ok=True)
    out.to_csv(STAB / "tableS16_limma_parity.tsv", sep="\t", index=False)
    print(f"{len(out)} contrasts, {int(out.n_pathways.sum()):,} pathway tests")
    print(f"  worst |d effect|   = {out.effect_max_abs_diff.max():.3e}")
    print(f"  worst |d t|        = {out.t_max_abs_diff.max():.3e}")
    print(f"  worst |d p|        = {out.p_max_abs_diff.max():.3e}")
    print(f"  worst d0 rel err   = {out.d0_rel_err.max():.3e}")
    print(f"  worst s0^2 rel err = {out.s0sq_rel_err.max():.3e}")
    print(f"  contrasts with any p<0.05 call differing: "
          f"{int((out.n_disagreeing_calls > 0).sum())}")
    return out


def run_all() -> None:
    """Fig. S12 and Table S16, from one limma pass."""
    parity = compare_with_limma()
    figure_s12(parity)
    table_s16(parity)
