"""Core PCA computation over a term x sample score matrix.

Ports the math from pca_analysis.ipynb section 3 / _pca_rollout.run_pca_for_sample:

  1. z-score every term (feature) across observations (ddof=0); drop constants.
  2. transpose -> observations x terms; observations are the PCA points.
  3. fit ``PCA(K=min(n_obs-1, K_MAX), random_state=42)``.

Observations are the columns of the score matrix within one *scope*.  Labels /
sub-labels / sizes are optional per-observation annotations used only for
figures + the signature table (they never affect the numeric result).
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from ..config import PCAConfig
from ..contracts import Grouping, PCAResult, ScoreMatrix
from .signatures import compute_signatures


def run_pca(score_df: pd.DataFrame,
            labels: pd.Series,
            config: PCAConfig,
            *,
            scope: str = "__all__",
            display_labels: Optional[pd.Series] = None,
            sublabels: Optional[pd.Series] = None,
            sizes: Optional[pd.Series] = None) -> Optional[PCAResult]:
    """Run PCA on one scope's ``terms x observations`` numeric matrix.

    Returns ``None`` if the scope has too few observations / variable features.
    """
    X = score_df.fillna(0.0)
    n_features, n_obs = X.shape
    if n_obs < config.min_observations or n_features < config.min_features:
        return None

    obs_ids = list(X.columns)
    # display labels for observations (default: the raw column id)
    if display_labels is not None:
        disp = [str(display_labels.get(o, o)) for o in obs_ids]
    else:
        disp = [str(o) for o in obs_ids]

    # ── z-score per feature across observations; drop constant features ──
    sds = X.std(axis=1, ddof=0)
    keep = sds > 1e-12
    Xn = X.loc[keep].sub(X.loc[keep].mean(axis=1), axis=0).div(sds[keep], axis=0)
    if Xn.shape[0] < config.min_features:
        return None

    K = min(n_obs - 1, config.k_max)
    pca = PCA(n_components=K, random_state=config.random_state)
    scores = pca.fit_transform(Xn.values.T)         # [n_obs x K]
    loadings = pca.components_                       # [K x n_features]
    var_expl = pca.explained_variance_ratio_        # [K]

    pc_names = [f"PC{k + 1}" for k in range(K)]
    scores_df = pd.DataFrame(scores, index=disp, columns=pc_names)
    loadings_df = pd.DataFrame(loadings, index=pc_names, columns=Xn.index)

    per_obs, sig_long = compute_signatures(
        scores_df, loadings_df,
        top_pcs_per_obs=config.top_pcs_per_cluster,
        top_pathways=config.top_pathways_per_cluster,
    )

    # align annotation series to the display-labelled observations
    def _series_on_disp(src: Optional[pd.Series], default="") -> pd.Series:
        vals = [(src.get(o, default) if src is not None else default) for o in obs_ids]
        name = src.name if src is not None else None
        return pd.Series(vals, index=disp, name=name)

    labels_disp = _series_on_disp(labels)
    sub_disp = _series_on_disp(sublabels)
    sizes_arr = (np.array([float(sizes.get(o, 1.0)) for o in obs_ids], dtype=float)
                 if sizes is not None else None)

    return PCAResult(
        scope=scope,
        scores=scores_df,
        loadings=loadings_df,
        variance_explained=var_expl,
        signatures=sig_long,
        labels=labels_disp,
        sublabels=sub_disp,
        sizes=sizes_arr,
    )


def run_pca_scopes(sm: ScoreMatrix, grouping: Grouping, config: PCAConfig,
                   *, display_col: Optional[str] = None,
                   sublabel_col: Optional[str] = None,
                   size_col: Optional[str] = None) -> list[PCAResult]:
    """Run one PCA per scope defined by ``grouping``.

    Optional metadata columns provide per-observation display labels, secondary
    colour labels, and dot sizes (looked up from ``sm.metadata``).
    """
    meta = sm.metadata.table
    disp = meta[display_col] if (display_col and display_col in meta) else None
    subl = meta[sublabel_col] if (sublabel_col and sublabel_col in meta) else None
    size = meta[size_col] if (size_col and size_col in meta) else None

    results: list[PCAResult] = []
    for scope, cols in grouping.scopes().items():
        cols = [c for c in cols if c in sm.data.columns]
        if len(cols) < config.min_observations:
            continue
        sub = sm.data.loc[:, cols]
        res = run_pca(
            sub, grouping.labels.reindex(cols), config,
            scope=scope,
            display_labels=disp, sublabels=subl, sizes=size,
        )
        if res is not None:
            results.append(res)
    return results
