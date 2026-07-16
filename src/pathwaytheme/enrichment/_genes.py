"""Per-sample gene-list selection shared by list-based backends (EnrichR/GoSlim).

Continuous omic data has no intrinsic "gene list", so for each sample we derive
one from its column: the top-N features by value, or all features above a
threshold.  This is the input to over-representation style enrichment.
"""

from __future__ import annotations

import pandas as pd

from ..config import EnrichmentConfig


def select_genes(col: pd.Series, config: EnrichmentConfig) -> list[str]:
    """Pick the gene list for one sample column."""
    s = col.dropna()
    if config.gene_selection == "threshold":
        if config.threshold is None:
            raise ValueError("gene_selection='threshold' requires enrichment.threshold")
        return list(s.index[s >= config.threshold].astype(str))
    # default: top_n by value
    n = min(config.top_n, len(s))
    return list(s.sort_values(ascending=False).head(n).index.astype(str))
