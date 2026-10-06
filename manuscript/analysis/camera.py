"""Supplementary Methods §12 -- how much of the term count survives set correlation.

PathwayTheme scores each gene set per sample and then runs a moderated t on the
score matrix, one row per set, with Benjamini-Hochberg across sets.  That
treats the sets as independent features, and over this cohort's GO:BP scores
they are not: :func:`term_correlation` measures how correlated the score rows
are and how many independent tests they amount to (Li & Ji).  The correlation
is not set overlap -- the same function reports the gene-set Jaccard over the
same terms -- so it is co-regulation of the underlying programmes, which no
membership-based correction would touch.

`limma::camera` is the competitive test that does model it, inflating each
set's variance by 1 + (m - 1) * rho.  :func:`camera_check` runs it on the
gene-level matrix with the same gene sets, the same samples and the same
contrasts, so the significant-term counts can be read against a number that
has paid for the correlation.

The two tests answer different questions (self-contained on the scores against
competitive on the genes), so nothing in the manuscript is recomputed from
camera; it bounds how the reported counts should be read.

R with the limma package is required for :func:`camera_check`; the
``Rscript`` executable is taken from the environment variable ``RSCRIPT``,
else from ``PATH``.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .limma_parity import rscript
from .paths import GENESETS, MYCN, PIPE, PROC, RESULTS, SCORES

R_DIR = Path(__file__).resolve().parent / "R"
OUT = RESULTS / "camera_check"

# (label, run directory holding our side, the contrast, the metadata column,
#  the sample type the contrast is restricted to)
CASES = [
    ("TUMOR vs NORMAL", "gobp_tumor_normal", ("TUMOR", "NORMAL"),
     "sample_type", None),
    ("CELL_LINE vs NORMAL", "gobp_tumor_normal", ("CELL_LINE", "NORMAL"),
     "sample_type", None),
    ("CELL_LINE vs TUMOR", "gobp_tumor_normal", ("CELL_LINE", "TUMOR"),
     "sample_type", None),
    ("MYCN.A vs MYCN.NA", None, ("MYCN.A", "MYCN.NA"), "mycn_status", "TUMOR"),
]

# term_correlation: how many score rows, drawn how
N_TERMS, SEED = 900, 0


def _ours(run: str | None, case: str, ref: str) -> pd.DataFrame:
    """PathwayTheme's score-level result for one contrast."""
    if run is None:
        p = (MYCN / "results" / "gobp_nbl_mycn" / "GOBP"
             / "GOBP_differential_moderated_t.tsv")
    else:
        p = PIPE / run / "GOBP" / "GOBP_differential_moderated_t.tsv"
    if not p.exists():
        raise FileNotFoundError(f"{p} is missing")
    d = pd.read_csv(p, sep="\t")
    d = d[d["comparison"] == f"{case}_vs_{ref}"]
    if d.empty:
        raise RuntimeError(f"{p.name} has no comparison {case}_vs_{ref}")
    return d.set_index("pathway")


def _groups(column: str, subset: str | None) -> pd.Series:
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0)
    if column == "mycn_status":
        call = pd.read_csv(MYCN / "mycn_status.tsv", sep="\t", index_col=0)
        s = call["mycn_status"].reindex(meta.index)
        if subset:
            s = s.where(meta["sample_type"] == subset)
        return s
    return meta[column]


def _run_camera(labels: pd.Series, case: str, ref: str) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        gpath, out = tmp / "g.tsv", tmp / "o.tsv"
        keep = labels[labels.isin([case, ref])]
        pd.DataFrame({"sample": keep.index,
                      "group": keep.to_numpy()}).to_csv(gpath, sep="\t",
                                                        index=False)
        proc = subprocess.run(
            [rscript(), str(R_DIR / "camera_reference.R"),
             str(PROC / "expression_log2.tsv"), str(gpath),
             str(GENESETS / "GOBP.gmt"), case, ref, str(out)],
            capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"camera failed for {case} vs {ref}:\n"
                               f"{proc.stderr[-2000:]}")
        print("   " + proc.stdout.strip())
        return pd.read_csv(out, sep="\t").set_index("pathway")


def camera_check() -> pd.DataFrame:
    """Supplementary Methods §12: limma::camera against PathwayTheme's term calls.

    Writes ``results/camera_check/camera_summary.tsv`` -- per contrast, how
    many of the terms called significant here camera also calls, and how often
    the two agree on direction -- and one ``camera_<case>_vs_<ref>.tsv`` per
    contrast with both sides term by term.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, run, (case, ref), column, subset in CASES:
        print(label)
        cam = _run_camera(_groups(column, subset), case, ref)
        us = _ours(run, case, ref)
        shared = us.index.intersection(cam.index)
        u, c = us.loc[shared], cam.loc[shared]

        us_sig = u["fdr"] < 0.05
        cam_sig = c["FDR"] < 0.05
        # camera reports a direction word, ours a signed effect
        same_way = (c["Direction"] == np.where(u["effect"] > 0, "Up", "Down"))
        both = us_sig & cam_sig
        rows.append({
            "contrast": label, "n_sets": len(shared),
            "n_sig_pathwaytheme": int(us_sig.sum()),
            "n_sig_camera": int(cam_sig.sum()),
            "n_sig_both": int(both.sum()),
            "pct_of_ours_kept": round(100 * float(both.sum())
                                      / max(int(us_sig.sum()), 1), 1),
            "direction_agreement_on_shared_sig":
                round(float(same_way[both].mean()), 4) if int(both.sum()) else np.nan,
            "spearman_rank": round(float(
                pd.Series(-np.log10(u["p_value"].clip(1e-300))).corr(
                    pd.Series(-np.log10(c["PValue"].clip(1e-300))),
                    method="spearman")), 4)})
        pd.concat([u[["effect", "statistic", "p_value", "fdr"]],
                   c[["NGenes", "Direction", "PValue", "FDR"]]
                   .add_prefix("camera_")], axis=1).to_csv(
            OUT / f"camera_{case}_vs_{ref}.tsv", sep="\t")

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "camera_summary.tsv", sep="\t", index=False)
    return res


def _li_ji_meff(corr: np.ndarray) -> float:
    """Li & Ji (2005) effective number of independent tests.

    M_eff = sum_i [ I(|l_i| >= 1) + (|l_i| - floor(|l_i|)) ] over the
    eigenvalues l_i of the correlation matrix.
    """
    lam = np.abs(np.linalg.eigvalsh(corr))
    return float(np.sum((lam >= 1) + (lam - np.floor(lam))))


def term_correlation() -> pd.DataFrame:
    """Supplementary Methods §12: correlation among GO:BP score rows and Li-Ji M_eff.

    Draws 900 of the scored GO:BP terms (every term that passed the size
    filter) at random without replacement, correlates their ssGSEA score rows
    (Pearson) across all samples of the cached score matrix, and reports the
    mean pairwise r, the mean |r|, and the Li-Ji effective number of
    independent tests among them.  For the same terms it reports the mean
    pairwise Jaccard of their gene sets (members as listed in the GMT) and the
    fraction of pairs sharing at least one gene, which separates co-regulation
    from overlap.  Writes ``results/camera_check/term_correlation.tsv``.
    """
    sm = pd.read_csv(SCORES / "GOBP_ssgsea_scores.tsv", sep="\t", index_col=0)
    rows = np.random.default_rng(SEED).choice(len(sm), N_TERMS, replace=False)
    names = sm.index[rows]
    corr = np.corrcoef(sm.iloc[rows].to_numpy(float))
    iu = np.triu_indices(N_TERMS, 1)

    members = {}
    with open(GENESETS / "GOBP.gmt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            members[parts[0]] = set(parts[2:])
    sets = [members[n] for n in names]
    jac, shared = [], 0
    for i in range(N_TERMS):
        for j in range(i + 1, N_TERMS):
            inter = len(sets[i] & sets[j])
            jac.append(inter / len(sets[i] | sets[j]))
            shared += inter > 0

    out = pd.DataFrame([{
        "n_terms_scored": len(sm), "n_samples": sm.shape[1],
        "n_terms_drawn": N_TERMS, "seed": SEED,
        "mean_r": round(float(corr[iu].mean()), 4),
        "mean_abs_r": round(float(np.abs(corr[iu]).mean()), 4),
        "li_ji_meff": round(_li_ji_meff(corr), 4),
        "mean_jaccard": round(float(np.mean(jac)), 4),
        "frac_pairs_sharing_a_gene": round(shared / len(jac), 4),
    }])
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "term_correlation.tsv", sep="\t", index=False)
    return out


def run_all() -> None:
    """The term correlation, then the camera comparison (the slow step)."""
    term_correlation()
    camera_check()
