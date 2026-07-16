"""EnrichR backend: per-sample over-representation -> term x sample matrix.

For each sample we take a gene list (see :mod:`._genes`) and run enrichment:

* **offline** (default, deterministic, testable): ``gseapy.enrich`` against a
  local ``.gmt`` (``enrichment.gmt_path``) using the full feature set as the
  statistical background.
* **online**: ``gseapy.enrichr`` against a named Enrichr library
  (``enrichment.enrichr_library``) — requires internet.

The per-sample results are pivoted into ``terms x samples`` using either
``-log10(Adjusted P-value)`` (default) or the Enrichr ``Combined Score``.
Terms absent for a sample are 0 (not enriched).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import EnrichmentConfig
from ..contracts import FeatureMatrix
from ..io.genesets import parse_gmt
from ._genes import select_genes
from .base import EnrichmentBackend, register


def _score_column(res2d: pd.DataFrame, how: str) -> pd.Series:
    """Turn one enrich result table into a term -> score Series."""
    df = res2d.copy()
    term_col = "Term" if "Term" in df.columns else df.columns[0]
    if how == "combined_score":
        col = "Combined Score"
        vals = pd.to_numeric(df[col], errors="coerce")
    else:  # -log10_padj
        col = "Adjusted P-value"
        p = pd.to_numeric(df[col], errors="coerce").clip(lower=1e-300)
        vals = -np.log10(p)
    return pd.Series(vals.values, index=df[term_col].astype(str)).groupby(level=0).max()


@register
class EnrichrBackend(EnrichmentBackend):
    name = "enrichr"

    def run(self, fm: FeatureMatrix, config: EnrichmentConfig) -> pd.DataFrame:
        import gseapy as gp

        offline = bool(config.gmt_path)
        if offline:
            gene_sets = parse_gmt(config.gmt_path)
        elif config.enrichr_library:
            gene_sets = config.enrichr_library
        else:
            raise ValueError("EnrichR needs enrichment.gmt_path (offline) or "
                             "enrichment.enrichr_library (online)")

        background = list(fm.features.astype(str)) if offline else None
        per_sample: dict[str, pd.Series] = {}
        for sample in fm.samples:
            genes = select_genes(fm.data[sample], config)
            if not genes:
                continue
            try:
                if offline:
                    res = gp.enrich(gene_list=genes, gene_sets=gene_sets,
                                    background=background, outdir=None, verbose=False)
                else:
                    res = gp.enrichr(gene_list=genes, gene_sets=gene_sets,
                                     organism=config.organism, outdir=None, no_plot=True)
                r2 = res.res2d
                if r2 is None or len(r2) == 0:
                    continue
                per_sample[sample] = _score_column(r2, config.enrichr_score)
            except Exception as e:  # keep the run going; report per-sample failures
                print(f"  [enrichr] {sample}: {e}", flush=True)

        if not per_sample:
            raise RuntimeError("EnrichR produced no results for any sample")
        mat = pd.DataFrame(per_sample).reindex(columns=[c for c in fm.samples
                                                        if c in per_sample])
        return mat.fillna(0.0)
