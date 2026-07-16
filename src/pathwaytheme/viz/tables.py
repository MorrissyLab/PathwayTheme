"""Write the 3 PCA summary TSVs (ports pca_analysis.ipynb section 12)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..contracts import PCAResult
from ..config import PCAConfig
from ..pca.signatures import top_pathways_per_pc


def write_pca_tables(result: PCAResult, out_dir: str | Path, prefix: str,
                     config: PCAConfig) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    # 1. variance explained
    p = out_dir / f"{prefix}_pca_variance_explained.tsv"
    pd.DataFrame({
        "PC": result.pc_names,
        "variance_explained": result.variance_explained,
        "cumulative_variance": np.cumsum(result.variance_explained),
    }).to_csv(p, sep="\t", index=False)
    written.append(p)

    # 2. top pathways per PC
    p = out_dir / f"{prefix}_pca_top_pathways_per_pc.tsv"
    top_pathways_per_pc(result.loadings, result.variance_explained,
                        config.top_pathways_per_pc).to_csv(p, sep="\t", index=False)
    written.append(p)

    # 3. per-cluster signatures
    p = out_dir / f"{prefix}_pca_per_cluster_signatures.tsv"
    result.signatures.to_csv(p, sep="\t", index=False)
    written.append(p)

    return written
