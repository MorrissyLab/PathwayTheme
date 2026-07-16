"""Data contracts that flow between pipeline stages.

The whole package is glued together by a small number of deliberately boring
pandas structures.  Keeping these thin (mostly ``@dataclass`` wrappers around a
``DataFrame``) means every stage is independently testable and swappable.

Conventions
-----------
* A **sample** is an observation column: it may be a bulk sample, a
  pseudobulk ``sample|cluster`` column, or any other unit the user chose.
* ``FeatureMatrix``  : rows = features (genes), cols = samples.
* ``ScoreMatrix``    : rows = terms (pathways / GO terms), cols = samples.
  Every enrichment backend returns this identical shape.
* ``SampleMetadata`` : indexed by sample id, one row per column of the matrices.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class SampleMetadata:
    """Per-sample attributes, indexed by sample id (matrix column labels)."""

    table: pd.DataFrame

    def __post_init__(self) -> None:
        self.table = self.table.copy()
        self.table.index = self.table.index.astype(str)

    @property
    def samples(self) -> list[str]:
        return list(self.table.index)

    def align_to(self, columns) -> "SampleMetadata":
        """Return metadata restricted & reordered to ``columns``."""
        cols = [str(c) for c in columns]
        missing = [c for c in cols if c not in self.table.index]
        if missing:
            raise KeyError(f"{len(missing)} sample(s) missing from metadata, "
                           f"e.g. {missing[:3]}")
        return SampleMetadata(self.table.loc[cols])

    def get(self, column: str) -> pd.Series:
        if column not in self.table.columns:
            raise KeyError(f"metadata column {column!r} not found; available: "
                           f"{list(self.table.columns)}")
        return self.table[column]

    def has(self, column: str) -> bool:
        return column in self.table.columns


@dataclass
class FeatureMatrix:
    """features (genes) x samples.  The input to enrichment."""

    data: pd.DataFrame
    metadata: SampleMetadata

    def __post_init__(self) -> None:
        self.data = self.data.copy()
        self.data.columns = self.data.columns.astype(str)
        self.data.index = self.data.index.astype(str)

    @property
    def features(self) -> pd.Index:
        return self.data.index

    @property
    def samples(self) -> list[str]:
        return list(self.data.columns)

    @property
    def shape(self):
        return self.data.shape


@dataclass
class ScoreMatrix:
    """terms (pathways) x samples.  The unified output of every backend."""

    data: pd.DataFrame
    metadata: SampleMetadata
    backend: str = "unknown"
    geneset: str = "unknown"

    def __post_init__(self) -> None:
        self.data = self.data.copy()
        self.data.columns = self.data.columns.astype(str)
        self.data.index = self.data.index.astype(str)

    @property
    def terms(self) -> pd.Index:
        return self.data.index

    @property
    def samples(self) -> list[str]:
        return list(self.data.columns)

    @property
    def shape(self):
        return self.data.shape

    def subset_samples(self, columns) -> "ScoreMatrix":
        cols = [str(c) for c in columns]
        return ScoreMatrix(self.data.loc[:, cols], self.metadata.align_to(cols),
                           backend=self.backend, geneset=self.geneset)


@dataclass
class Grouping:
    """Maps each sample -> a comparison label, and (optionally) a scope.

    ``labels`` drives figure colouring and comparison.  ``scope`` optionally
    partitions samples into independent PCA runs (e.g. one PCA per biological
    sample); when ``scope`` is None the whole matrix is a single scope.
    """

    labels: pd.Series
    scope: Optional[pd.Series] = None
    label_name: str = "group"
    scope_name: Optional[str] = None

    def __post_init__(self) -> None:
        self.labels = self.labels.copy()
        self.labels.index = self.labels.index.astype(str)
        if self.scope is not None:
            self.scope = self.scope.copy()
            self.scope.index = self.scope.index.astype(str)

    def scopes(self) -> dict[str, list[str]]:
        """Return {scope_value: [sample ids]}.  Single '__all__' scope if none."""
        if self.scope is None:
            return {"__all__": list(self.labels.index)}
        out: dict[str, list[str]] = {}
        for sample, sc in self.scope.items():
            out.setdefault(str(sc), []).append(str(sample))
        return out


@dataclass
class PCAResult:
    """Output of one PCA run over a single scope."""

    scope: str
    scores: pd.DataFrame          # observations x PCs
    loadings: pd.DataFrame        # PCs x features
    variance_explained: np.ndarray  # (K,)
    signatures: pd.DataFrame      # long: cluster, top_pcs, rank, pathway, signature_score, direction
    labels: pd.Series             # observation -> comparison label (aligned to scores.index)
    sublabels: pd.Series = field(default_factory=lambda: pd.Series(dtype=object))
    sizes: Optional[np.ndarray] = None   # per-observation weight (e.g. n_cells) for dot sizing

    @property
    def pc_names(self) -> list[str]:
        return list(self.scores.columns)

    @property
    def K(self) -> int:
        return self.scores.shape[1]
