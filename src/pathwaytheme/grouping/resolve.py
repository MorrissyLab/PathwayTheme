"""Resolve a comparison Grouping from one of three sources.

modes
-----
* ``target``   : colour/compare by a named metadata column (e.g. classification).
* ``existing`` : use cluster labels already present in the metadata column.
* ``auto``     : cluster the score matrix when no labels exist.

An optional ``scope_col`` partitions samples into independent PCA runs (e.g.
one PCA per biological sample).  Without it, all samples form a single scope.
"""

from __future__ import annotations

import pandas as pd

from ..config import GroupingConfig
from ..contracts import Grouping, ScoreMatrix
from .cluster import auto_cluster


def resolve_grouping(sm: ScoreMatrix, config: GroupingConfig) -> Grouping:
    meta = sm.metadata.table
    samples = list(sm.data.columns)

    mode = config.mode.lower()
    if mode in ("target", "existing"):
        if not config.target_col:
            raise ValueError(f"grouping.mode={mode!r} requires grouping.target_col")
        if config.target_col not in meta.columns:
            raise KeyError(f"target_col {config.target_col!r} not in metadata "
                           f"columns {list(meta.columns)}")
        labels = meta.loc[samples, config.target_col].astype(str)
        labels.name = config.target_col
    elif mode == "auto":
        labels = auto_cluster(sm.data, config)
        labels = labels.reindex(samples)
    else:
        raise ValueError(f"Unknown grouping.mode {config.mode!r}; "
                         f"expected 'target' | 'existing' | 'auto'")

    scope = None
    if config.scope_col:
        if config.scope_col not in meta.columns:
            raise KeyError(f"scope_col {config.scope_col!r} not in metadata")
        scope = meta.loc[samples, config.scope_col].astype(str)

    return Grouping(labels=labels, scope=scope,
                    label_name=str(labels.name or "group"),
                    scope_name=config.scope_col)
