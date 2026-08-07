"""Generate shared fixtures and dump every Python stage output for R to check.

Run from the repo root with the project venv:

    .venv/Scripts/python.exe r/parity/run_python.py

Writes ``r/parity/data/`` (inputs both implementations read) and
``r/parity/python/`` (reference outputs).  ``r/parity/compare.R`` then runs the
same stages in R and reports the agreement.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import pathwaytheme as pt                                    # noqa: E402
from pathwaytheme.enrichment.coverage import geneset_coverage  # noqa: E402
from pathwaytheme.diff.stats import fit_fdist, trigamma_inverse  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "python"


def build_fixtures() -> None:
    """Deterministic inputs: a dense matrix, a tie-heavy matrix, and a gmt."""
    rng = np.random.default_rng(20260802)
    genes = [f"G{i:04d}" for i in range(600)]
    samples = [f"S{i:02d}" for i in range(16)]

    # --- gene sets: 40 terms of 8-60 genes, plus two that barely overlap ---
    lines = []
    for t in range(40):
        k = int(rng.integers(8, 60))
        members = sorted(rng.choice(genes, size=k, replace=False).tolist())
        lines.append(f"TERM_{t:02d}\tdesc\t" + "\t".join(members))
    lines.append("TERM_SPARSE\tdesc\t" + "\t".join(
        [genes[0], genes[1]] + [f"ABSENT{i}" for i in range(20)]))
    lines.append("TERM_ABSENT\tdesc\t" + "\t".join(
        [f"NOWHERE{i}" for i in range(15)]))
    (DATA / "genesets.gmt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # --- continuous matrix with group structure ---
    x = rng.normal(0, 1, size=(len(genes), len(samples)))
    x[:120, :8] += 1.6
    expr = pd.DataFrame(x, index=genes, columns=samples)
    expr.to_csv(DATA / "expr.tsv", sep="\t")

    # --- tie-heavy matrix: ~45% exact zeros, to pin down rank tie handling ---
    y = rng.normal(0, 1, size=(len(genes), len(samples)))
    y[rng.random(y.shape) < 0.45] = 0.0
    y[:120, :8] += 1.0
    pd.DataFrame(y, index=genes, columns=samples).to_csv(
        DATA / "expr_ties.tsv", sep="\t")

    # --- metadata ---
    meta = pd.DataFrame({
        "sample_id": ["A"] * 8 + ["B"] * 8,
        "cluster_id": [str(i) for i in range(8)] * 2,
        "grp": ["hi"] * 8 + ["lo"] * 8,
        "arm": (["ctrl"] * 4 + ["treat"] * 4) * 2,
        "cell_type": ["tumor", "macrophage", "T-cell", "fibroblast"] * 4,
        "n_cells": rng.integers(20, 500, size=len(samples)),
    }, index=samples)
    meta.to_csv(DATA / "meta.tsv", sep="\t")

    # --- term -> category map ---
    cats = pd.DataFrame({
        "pathway": [f"TERM_{t:02d}" for t in range(40)],
        "category": [["metabolism", "signalling", "immune"][t % 3]
                     for t in range(40)],
    })
    cats.to_csv(DATA / "categories.tsv", sep="\t", index=False)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    build_fixtures()

    gmt = str(DATA / "genesets.gmt")

    # ---- stage 1-2: input + ssGSEA (both the clean and the tie-heavy matrix)
    for tag, path in (("", "expr.tsv"), ("_ties", "expr_ties.tsv")):
        fm = pt.load_matrix(str(DATA / path), metadata=str(DATA / "meta.tsv"))
        sm = pt.enrich(fm, backend="ssgsea", gmt=gmt, geneset="SYN", threads=1)
        sm.data.to_csv(OUT / f"ssgsea{tag}.tsv", sep="\t")

    fm = pt.load_matrix(str(DATA / "expr.tsv"), metadata=str(DATA / "meta.tsv"))
    sm = pt.enrich(fm, backend="ssgsea", gmt=gmt, geneset="SYN", threads=1)

    # ---- gene-set coverage ----
    from pathwaytheme.io.genesets import parse_gmt
    geneset_coverage(fm.features, parse_gmt(gmt), 5, 5000).to_csv(
        OUT / "coverage.tsv", sep="\t", index=False)

    # ---- over-representation backend ----
    ora = pt.enrich(fm, backend="enrichr", gmt=gmt, geneset="SYN", top_n=120)
    ora.data.to_csv(OUT / "ora.tsv", sep="\t")
    goslim = pt.enrich(fm, backend="goslim", gmt=gmt, geneset="SYN", top_n=120)
    goslim.data.to_csv(OUT / "goslim.tsv", sep="\t")

    # ---- stage 3-4: grouping + PCA ----
    groups = pt.group(sm, mode="target", target_col="grp")
    res = pt.pca(sm, groups, size_col="n_cells")[0]
    pd.DataFrame({"PC": res.pc_names,
                  "variance_explained": res.variance_explained}).to_csv(
        OUT / "pca_variance.tsv", sep="\t", index=False)
    res.scores.to_csv(OUT / "pca_scores.tsv", sep="\t")
    res.loadings.to_csv(OUT / "pca_loadings.tsv", sep="\t")
    res.signatures.to_csv(OUT / "pca_signatures.tsv", sep="\t", index=False)

    # scoped PCA (one per sample_id) exercises the scope machinery
    groups_scoped = pt.group(sm, mode="target", target_col="grp",
                             scope_col="sample_id")
    scoped = pt.pca(sm, groups_scoped)
    pd.concat([pd.DataFrame({"scope": r.scope, "PC": r.pc_names,
                             "variance_explained": r.variance_explained})
               for r in scoped]).to_csv(OUT / "pca_scoped_variance.tsv",
                                        sep="\t", index=False)

    # ---- stage 5: differential, all three methods + contrast forms ----
    for method in ("welch", "mannwhitney", "moderated_t"):
        d = pt.diff(sm, groups, method=method)
        d.table.to_csv(OUT / f"diff_{method}.tsv", sep="\t", index=False)

    d_ref = pt.diff(sm, method="welch", group_col="arm", reference="ctrl")
    d_ref.table.to_csv(OUT / "diff_reference.tsv", sep="\t", index=False)
    d_con = pt.diff(sm, method="moderated_t", group_col="cell_type",
                    contrasts=[["tumor", "T-cell"]])
    d_con.table.to_csv(OUT / "diff_contrast.tsv", sep="\t", index=False)

    # ---- categories ----
    d = pt.diff(sm, groups, method="moderated_t")
    summary = pt.summarize_categories(d, str(DATA / "categories.tsv"),
                                      significance_col="fdr")
    summary.to_csv(OUT / "categories.tsv", sep="\t", index=False)

    # ---- limma helper functions, probed directly ----
    rng = np.random.default_rng(7)
    var = rng.chisquare(4, size=200) / 4
    d0, s0 = fit_fdist(var, np.full(200, 4.0))
    pd.DataFrame({
        "quantity": ["d0", "s0_sq"] + [f"trigamma_inverse({v})" for v in
                                       (0.01, 0.1, 0.5, 1.0, 5.0)],
        "value": [d0, s0] + [trigamma_inverse(v) for v in
                             (0.01, 0.1, 0.5, 1.0, 5.0)],
    }).to_csv(OUT / "limma_helpers.tsv", sep="\t", index=False)
    pd.DataFrame({"var": var}).to_csv(OUT / "limma_var.tsv", sep="\t",
                                      index=False)

    print(f"wrote python reference outputs to {OUT}")


if __name__ == "__main__":
    main()
