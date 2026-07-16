"""Matrix adapter: read a pre-quantified features x samples table.

This is the most generic "any omic" entry point.  The matrix may be gene
expression, protein abundance, methylation, etc. — the pipeline only requires
that rows are features (matchable to gene-set members) and columns are samples.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..config import InputConfig
from ..contracts import FeatureMatrix, SampleMetadata
from .base import InputAdapter


def _read_table(path: str | Path, sep: str) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in (".parquet", ".pq"):
        return pd.read_parquet(path)
    if path.suffix.lower() in (".csv",):
        return pd.read_csv(path, index_col=0)
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path, index_col=0)
    return pd.read_csv(path, sep=sep, index_col=0)


class MatrixAdapter(InputAdapter):
    def load(self, config: InputConfig) -> FeatureMatrix:
        if not config.matrix_path:
            raise ValueError("InputConfig.matrix_path is required for kind='matrix'")
        data = _read_table(config.matrix_path, config.sep)
        data = data.apply(pd.to_numeric, errors="coerce")
        data.columns = data.columns.astype(str)
        data.index = data.index.astype(str)

        if config.metadata_path:
            meta = _read_table(config.metadata_path, config.sep)
            meta.index = meta.index.astype(str)
            # keep only samples that exist in the matrix, preserve matrix order
            common = [c for c in data.columns if c in meta.index]
            if len(common) < data.shape[1]:
                missing = [c for c in data.columns if c not in meta.index]
                raise KeyError(f"{len(missing)} matrix columns absent from metadata, "
                               f"e.g. {missing[:3]}")
            meta = meta.loc[data.columns]
        else:
            # minimal metadata: sample id as its own attribute
            meta = pd.DataFrame({"sample": data.columns}, index=data.columns)

        return FeatureMatrix(data, SampleMetadata(meta))
