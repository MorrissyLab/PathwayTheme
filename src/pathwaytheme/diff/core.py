"""Differential pathway analysis over a ScoreMatrix + grouping.

Runs a two-group test (see :mod:`pathwaytheme.diff.stats`) for one or more
comparisons and returns a single long :class:`~pathwaytheme.contracts.DiffResult`
table.  This is the Python analogue of the ASPS R workflow's limma /
``wilcox_pathway_stats`` step: compare pathway scores between groups, FDR-correct,
rank.

Comparison selection
--------------------
* explicit ``contrasts=[[case, ref], ...]``  — run exactly those, OR
* a ``reference`` group                       — every other group vs reference, OR
* neither                                      — one-vs-rest for every group.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from ..config import DiffConfig
from ..contracts import DiffResult, Grouping, ScoreMatrix
from .stats import benjamini_hochberg, run_test


def _labels_for(sm: ScoreMatrix, config: DiffConfig,
                grouping: Optional[Grouping]) -> pd.Series:
    """Resolve the per-sample group label series used for comparisons."""
    if config.group_col:
        if not sm.metadata.has(config.group_col):
            raise KeyError(f"diff.group_col {config.group_col!r} not in metadata")
        labels = sm.metadata.get(config.group_col).reindex(sm.samples).astype(str)
        labels.name = config.group_col
        return labels
    if grouping is not None:
        return grouping.labels.reindex(sm.samples).astype(str)
    raise ValueError("differential analysis needs diff.group_col or a grouping")


def _plan_comparisons(config: DiffConfig, groups: list[str]) -> list[tuple[str, object]]:
    """Return [(case, reference_or_REST), ...].  Reference ``None`` means one-vs-rest."""
    if config.contrasts:
        return [(str(c), str(r)) for c, r in config.contrasts]
    if config.reference is not None:
        ref = str(config.reference)
        return [(g, ref) for g in groups if g != ref]
    return [(g, None) for g in groups]


def run_diff(sm: ScoreMatrix, config: DiffConfig,
             grouping: Optional[Grouping] = None) -> DiffResult:
    """Compare pathway scores between groups; return a stacked DiffResult table."""
    labels = _labels_for(sm, config, grouping)
    groups = sorted(pd.unique(labels.dropna()))
    plan = _plan_comparisons(config, groups)

    data = sm.data                          # pathways x samples
    pathways = data.index.to_numpy()
    frames: list[pd.DataFrame] = []

    for case, ref in plan:
        case_cols = labels.index[labels == case]
        if ref is None:
            ref_cols = labels.index[labels != case]
            ref_name = "rest"
        else:
            ref_cols = labels.index[labels == ref]
            ref_name = str(ref)
        case_cols = [c for c in case_cols if c in data.columns]
        ref_cols = [c for c in ref_cols if c in data.columns]
        if len(case_cols) < config.min_group_size or len(ref_cols) < config.min_group_size:
            continue

        case_mat = data.loc[:, case_cols].to_numpy(dtype=float)
        ref_mat = data.loc[:, ref_cols].to_numpy(dtype=float)
        effect, stat, pval = run_test(config.method, case_mat, ref_mat)
        fdr = benjamini_hochberg(pval)

        frame = pd.DataFrame({
            "comparison": f"{case}_vs_{ref_name}",
            "pathway": pathways,
            "method": config.method,
            "n_case": len(case_cols),
            "n_reference": len(ref_cols),
            "mean_case": np.nanmean(case_mat, axis=1),
            "mean_reference": np.nanmean(ref_mat, axis=1),
            "effect": effect,
            "statistic": stat,
            "p_value": pval,
            "fdr": fdr,
        })
        frame["direction"] = np.where(frame["effect"] >= 0, "up", "down")
        frame = frame.sort_values(["fdr", "p_value"], na_position="last")
        frames.append(frame)

    table = (pd.concat(frames, ignore_index=True) if frames
             else pd.DataFrame(columns=[
                 "comparison", "pathway", "method", "n_case", "n_reference",
                 "mean_case", "mean_reference", "effect", "statistic",
                 "p_value", "fdr", "direction"]))
    return DiffResult(table=table, method=config.method,
                      group_col=str(labels.name or "group"))
