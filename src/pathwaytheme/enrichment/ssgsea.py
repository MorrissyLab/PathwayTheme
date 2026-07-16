"""ssGSEA backend (continuous per-sample NES matrix).

Ports ``run_ssgsea`` from ssgsea_analysis_fixed.ipynb: gseapy.ssgsea with
rank normalisation and weight 0.25, pivoted to ``terms x samples``.  This is
the backend the PCA pipeline was originally designed around.
"""

from __future__ import annotations

import pandas as pd

from ..config import EnrichmentConfig
from ..contracts import FeatureMatrix
from ..io.genesets import parse_gmt, harmonize_case
from .base import EnrichmentBackend, register


@register
class SsgseaBackend(EnrichmentBackend):
    name = "ssgsea"

    def run(self, fm: FeatureMatrix, config: EnrichmentConfig) -> pd.DataFrame:
        import gseapy as gp

        if not config.gmt_path:
            raise ValueError("ssGSEA requires enrichment.gmt_path (a .gmt file)")
        gs = parse_gmt(config.gmt_path)

        expr = fm.data.copy()
        # harmonise gene-symbol case to maximise overlap (ipynb logic)
        transform, gs = harmonize_case(expr.index, gs)
        expr.index = [transform(g) for g in expr.index]
        expr = expr.groupby(expr.index).mean()

        res = gp.ssgsea(
            data=expr, gene_sets=gs,
            sample_norm_method="rank", correl_norm_type="rank",
            outdir=None,
            min_size=config.min_gene_set_size, max_size=config.max_gene_set_size,
            threads=config.threads, weight=config.weight,
            ascending=False, permutation_num=0, no_plot=True, verbose=False,
        )
        r = res.res2d.copy()
        score_col = "NES" if "NES" in r.columns else ("ES" if "ES" in r.columns else None)
        if score_col is None:
            raise RuntimeError("gseapy.ssgsea returned no NES/ES column")
        wide = r.pivot(index="Term", columns="Name", values=score_col)
        wide = wide.apply(pd.to_numeric, errors="coerce")
        # restrict to gene-set terms actually requested, preserve sample order
        return wide.reindex(columns=[c for c in fm.samples if c in wide.columns])
