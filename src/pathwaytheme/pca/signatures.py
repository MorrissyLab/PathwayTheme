"""Per-observation pathway signatures from a fitted PCA.

Ports pca_analysis.ipynb section 4 / _pca_rollout.run_pca_for_sample:
for each observation, pick its top-N PCs by |score|, then

    signature(obs, pathway) = sum_{k in top PCs} score[obs, k] * loading[k, pathway]

Rank pathways by |signature|, keep the top M (sign preserved: + = up, - = down).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_signatures(scores_df: pd.DataFrame, loadings_df: pd.DataFrame,
                       top_pcs_per_obs: int, top_pathways: int):
    """Return ``(per_obs, signatures_long)``.

    per_obs : {obs_label: (top_pc_names, top_paths_series)}
    signatures_long : DataFrame [cluster, top_pcs, rank, pathway, signature_score, direction]
    """
    per_obs: dict = {}
    records: list[dict] = []
    for obs_label, s_vec in scores_df.iterrows():
        top_pcs = (s_vec.abs().sort_values(ascending=False)
                   .head(top_pcs_per_obs).index.tolist())
        sig = np.zeros(loadings_df.shape[1])
        for k in top_pcs:
            sig += s_vec[k] * loadings_df.loc[k].values
        sig_ser = pd.Series(sig, index=loadings_df.columns)
        top_idx = sig_ser.abs().sort_values(ascending=False).head(top_pathways).index
        top_paths = sig_ser.loc[top_idx]
        per_obs[obs_label] = (top_pcs, top_paths)
        for rank, (p, v) in enumerate(top_paths.items(), 1):
            records.append({
                "cluster": obs_label,
                "top_pcs": ",".join(top_pcs),
                "rank": rank,
                "pathway": p,
                "signature_score": float(v),
                "direction": "up" if v > 0 else "down",
            })
    return per_obs, pd.DataFrame(records)


def top_pathways_per_pc(loadings_df: pd.DataFrame, var_expl: np.ndarray,
                        top_n: int) -> pd.DataFrame:
    """Long table of the top-|loading| pathways per PC (ports the TSV logic)."""
    pc_names = list(loadings_df.index)
    rows: list[dict] = []
    for i, k in enumerate(pc_names):
        lv = loadings_df.loc[k]
        top = lv.abs().sort_values(ascending=False).head(top_n).index
        for rank, p in enumerate(top, 1):
            rows.append({
                "PC": k,
                "variance_explained": float(var_expl[i]),
                "rank": rank,
                "pathway": p,
                "loading": float(lv[p]),
                "direction": "up" if lv[p] > 0 else "down",
            })
    return pd.DataFrame(rows)
