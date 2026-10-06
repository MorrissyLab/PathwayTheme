"""Reshape the pipeline's outputs into the summary tables the paper ships.

Reads the per-run outputs under ``results/pipeline`` and writes, into
``results/tables``:

    contrast_counts.tsv      per contrast: terms tested, significant, up, down
                             (shipped as Table S5)
    theme_summary_full.tsv   per contrast: every GO-slim theme, significant and
                             not, with its mean effect and up/down split
                             (shipped as Table S6)
    pc_loadings_named.tsv    per run: the top positive and negative pathways
                             per component (shipped as Table S7)

Nothing here recomputes a statistic; it only reshapes what the pipeline wrote,
so every number in these tables can be traced to a pipeline file.
:mod:`analysis.figures` copies them into ``submission/tables``, so this module
runs first.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .paths import PIPE, TABLES

ALPHA = 0.05

# run directory -> the gene-set collection it scored
RUNS = {
    "hallmark_sampletype": "HALLMARK",
    "gobp_sampletype": "GOBP",
    "gobp_diagnosis": "GOBP",
    "gobp_tumor_normal": "GOBP",
    "gobp_matched": "GOBP",
}


def available_runs() -> dict[str, str]:
    """The runs of :data:`RUNS` that exist under ``results/pipeline``."""
    runs = {r: g for r, g in RUNS.items() if (PIPE / r).exists()}
    if not runs:
        raise FileNotFoundError(f"no pipeline runs under {PIPE}; run the pipeline first")
    return runs


def find(run: str, pattern: str) -> Path | None:
    """The first file under a run directory matching ``pattern``, if any."""
    hits = list((PIPE / run).rglob(pattern))
    return hits[0] if hits else None


def short_label(term: str) -> str:
    """A pathway name without its collection prefix, in lower case words."""
    for pre in ("GOBP_", "HALLMARK_"):
        if term.startswith(pre):
            term = term[len(pre):]
    return term.replace("_", " ").lower()


def contrast_counts(runs: dict[str, str] | None = None) -> pd.DataFrame:
    """Table S5 source: terms tested, significant (FDR < 0.05), up and down per contrast.

    Writes ``results/tables/contrast_counts.tsv``.
    """
    runs = runs or available_runs()
    rows = []
    for run, geneset in runs.items():
        p = find(run, f"{geneset}_differential_*.tsv")
        if p is None:
            continue
        d = pd.read_csv(p, sep="\t")
        for contrast, g in d.groupby("comparison"):
            sig = g[g["fdr"] < ALPHA]
            rows.append({
                "run": run, "geneset": geneset, "contrast": contrast,
                "n_case": int(g["n_case"].iloc[0]),
                "n_reference": int(g["n_reference"].iloc[0]),
                "n_tested": len(g), "n_significant": len(sig),
                "n_up": int((sig["effect"] > 0).sum()),
                "n_down": int((sig["effect"] < 0).sum()),
            })
    out = pd.DataFrame(rows)
    TABLES.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLES / "contrast_counts.tsv", sep="\t", index=False)
    return out


def theme_summary_full(runs: dict[str, str] | None = None) -> pd.DataFrame:
    """Table S6 source: every GO-slim theme of every contrast, both significance strata.

    The non-significant stratum is kept rather than dropped: the summarisation
    reports the members that did not reach significance alongside those that
    did, and a table cut to the significant themes could not show that.

    Writes ``results/tables/theme_summary_full.tsv``.
    """
    runs = runs or available_runs()
    full = []
    for run, geneset in runs.items():
        p = find(run, f"{geneset}_category_summary.tsv")
        if p is None:
            continue
        full.append(pd.read_csv(p, sep="\t").assign(run=run, geneset=geneset))
    every = pd.concat(full, ignore_index=True) if full else pd.DataFrame()
    if not every.empty:
        every = every[["run", "geneset", "comparison", "significance",
                       "category", "mean_effect", "n_pathways", "n_up",
                       "n_down"]].rename(columns={"comparison": "contrast",
                                                  "category": "theme",
                                                  "n_pathways": "n_terms"})
        every["mean_effect"] = every["mean_effect"].round(4)
    TABLES.mkdir(parents=True, exist_ok=True)
    every.to_csv(TABLES / "theme_summary_full.tsv", sep="\t", index=False)
    return every


def pc_loadings_named(runs: dict[str, str] | None = None,
                      per_side: int = 8) -> pd.DataFrame:
    """Table S7 source: the strongest positive and negative pathways per component.

    Taken from the pipeline's top-pathways-per-component table, ``per_side``
    at each end.  Writes ``results/tables/pc_loadings_named.tsv``.
    """
    runs = runs or available_runs()
    rows = []
    for run, geneset in runs.items():
        p = find(run, f"*{geneset}*_pca_top_pathways_per_pc.tsv")
        if p is None:
            continue
        t = pd.read_csv(p, sep="\t")
        for pc, g in t.groupby("PC", sort=False):
            up = g[g["loading"] > 0].nlargest(per_side, "loading")
            dn = g[g["loading"] < 0].nsmallest(per_side, "loading")
            for direction, sub in (("positive", up), ("negative", dn)):
                for rank, row in enumerate(sub.itertuples(), 1):
                    rows.append({
                        "run": run, "geneset": geneset, "component": pc,
                        "variance_explained": round(float(row.variance_explained), 5),
                        "axis": direction, "rank": rank,
                        "pathway": row.pathway, "label": short_label(row.pathway),
                        "loading": round(float(row.loading), 5),
                    })
    out = pd.DataFrame(rows)
    TABLES.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLES / "pc_loadings_named.tsv", sep="\t", index=False)
    return out


def run_all() -> None:
    """Write the three summary tables behind Tables S5-S7."""
    runs = available_runs()
    contrast_counts(runs)
    theme_summary_full(runs)
    pc_loadings_named(runs)
