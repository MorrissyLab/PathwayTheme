"""h5ad adapter: aggregate single-cell .h5ad files into per-cluster pseudobulk.

Ports the core pseudobulk logic from ssgsea_analysis_fixed.ipynb (cell 3):
for every ``.h5ad`` and every cluster within it, take the mean (or sum)
expression across cells -> one pseudobulk column ``<sample>|cl<cluster>``.

Domain-specific WorkingID normalisation from the notebook is intentionally
NOT included here — this adapter stays generic.  A project that needs it can
post-process the returned metadata.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import InputConfig
from ..contracts import FeatureMatrix, SampleMetadata
from .base import InputAdapter


class H5adAdapter(InputAdapter):
    def load(self, config: InputConfig) -> FeatureMatrix:
        try:
            import anndata as ad
        except ImportError as e:  # pragma: no cover
            raise ImportError("h5ad input requires anndata (`pip install "
                              "pathwaytheme[h5ad]`)") from e
        from scipy import sparse

        paths = self._resolve_paths(config)
        if not paths:
            raise FileNotFoundError("No .h5ad files found for the h5ad adapter")

        cluster_col, sample_col = config.cluster_col, config.sample_col
        agg = config.aggregate.lower()

        vecs: list[pd.Series] = []
        meta_rows: list[dict] = []
        for fp in paths:
            adata = ad.read_h5ad(fp)
            if cluster_col not in adata.obs.columns:
                continue
            obs = adata.obs.copy()
            obs[cluster_col] = obs[cluster_col].astype(str)
            if sample_col in obs.columns:
                sample_id = (pd.Series(obs[sample_col].astype(str))
                             .dropna().unique().tolist()[0])
            else:
                sample_id = fp.stem

            X = adata.X
            genes = pd.Index(adata.var_names.astype(str), name="gene")
            for cl in sorted(obs[cluster_col].unique(),
                             key=lambda x: (len(str(x)), str(x))):
                idx = np.where(obs[cluster_col].values == cl)[0]
                if idx.size == 0:
                    continue
                block = X[idx, :]
                arr = (np.asarray(block.todense()) if sparse.issparse(block)
                       else np.asarray(block))
                vec = arr.mean(axis=0) if agg == "mean" else arr.sum(axis=0)
                col = f"{sample_id}|cl{cl}"
                vecs.append(pd.Series(np.asarray(vec).ravel(), index=genes, name=col))
                meta_rows.append({
                    "sample_id": str(sample_id),
                    "cluster_id": str(cl),
                    "sample_cluster_id": col,
                    "modality": "snRNA_cluster",
                    "source_file": fp.name,
                    "n_cells": int(idx.size),
                })

        expr = pd.concat(vecs, axis=1).fillna(0.0).groupby(level=0).mean()
        meta = (pd.DataFrame(meta_rows)
                .drop_duplicates(subset=["sample_cluster_id"])
                .set_index("sample_cluster_id"))
        meta = meta.loc[expr.columns]
        return FeatureMatrix(expr, SampleMetadata(meta))

    @staticmethod
    def _resolve_paths(config: InputConfig) -> list[Path]:
        if config.h5ad_paths:
            return [Path(p) for p in config.h5ad_paths]
        if config.h5ad_dir:
            return sorted(Path(config.h5ad_dir).glob("*.h5ad"))
        return []
