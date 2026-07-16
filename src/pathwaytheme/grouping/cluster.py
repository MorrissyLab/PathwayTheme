"""Auto-clustering of samples in term space (used when no labels exist).

Observations are the score-matrix columns (samples); features are the terms.
Features are z-scored across samples, then KMeans / Agglomerative clustering is
run.  ``k`` is chosen by silhouette score over ``[k_min, k_max]`` unless a fixed
``n_clusters`` is given.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import GroupingConfig


def auto_cluster(score_data: pd.DataFrame, config: GroupingConfig) -> pd.Series:
    """Return a Series mapping each sample (column) -> cluster label ('c0', 'c1', ...)."""
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans, AgglomerativeClustering
    from sklearn.metrics import silhouette_score

    samples = list(score_data.columns)
    n = len(samples)
    if n < 2:
        return pd.Series(["c0"] * n, index=samples)

    X = score_data.fillna(0.0).T.values                # samples x terms
    X = StandardScaler().fit_transform(X)

    def _fit(k: int) -> np.ndarray:
        if config.method == "hierarchical":
            return AgglomerativeClustering(n_clusters=k).fit_predict(X)
        return KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(X)

    if config.n_clusters:
        labels = _fit(min(config.n_clusters, n))
    else:
        k_lo = max(2, config.k_min)
        k_hi = min(config.k_max, n - 1)
        best_score, labels = -1.0, None
        for k in range(k_lo, k_hi + 1):
            lab = _fit(k)
            if len(set(lab)) < 2:
                continue
            s = silhouette_score(X, lab)
            if s > best_score:
                best_score, labels = s, lab
        if labels is None:                 # too few samples to search a range
            labels = _fit(min(2, n))
    return pd.Series([f"c{int(l)}" for l in labels], index=samples, name="auto_cluster")
