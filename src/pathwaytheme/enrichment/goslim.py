"""GoSlim backend: coarse GO-slim category scores per sample.

Two modes, both producing a ``slim-term x sample`` matrix:

* **gmt mode** (default, offline, testable): a GO-slim ``.gmt`` where each term
  is a slim category and members are its genes.  Per sample, the score for a
  slim term is the fraction (or count) of the sample's selected genes that are
  members of that term.
* **obo mode** (goatools): if ``obo_path`` + ``slim_obo_path`` + ``gene2go_path``
  are given, genes are mapped to their GO annotations and rolled up to slim
  ancestors via ``goatools.mapslim``; the score is the fraction/count of the
  sample's genes hitting each slim term.
"""

from __future__ import annotations

import pandas as pd

from ..config import EnrichmentConfig
from ..contracts import FeatureMatrix
from ..io.genesets import parse_gmt
from ._genes import select_genes
from .base import EnrichmentBackend, register


def _score_from_membership(genes: list[str], term_members: dict[str, set[str]],
                           how: str) -> pd.Series:
    gset = set(genes)
    n = max(1, len(gset))
    out = {}
    for term, members in term_members.items():
        hit = len(gset & members)
        if hit:
            out[term] = (hit / n) if how == "fraction" else float(hit)
    return pd.Series(out, dtype=float)


def _build_obo_membership(config: EnrichmentConfig, universe: set[str]) -> dict[str, set[str]]:
    """Map genes -> GO-slim terms using goatools; return {slim_term: {genes}}."""
    from goatools.obo_parser import GODag
    from goatools.mapslim import mapslim
    from goatools.anno.idtogos_reader import IdToGosReader

    godag = GODag(config.obo_path)
    slimdag = GODag(config.slim_obo_path)
    # associations: gene_symbol -> set(GO ids).  gene2go_path is a 2-col
    # "id<tab>GO:...;GO:..." file (IdToGos format) keyed by gene symbol.
    assoc = IdToGosReader(config.gene2go_path, godag=godag).get_id2gos()

    term_members: dict[str, set[str]] = {}
    for gene, gos in assoc.items():
        if gene not in universe:
            continue
        slims: set[str] = set()
        for go in gos:
            if go not in godag:
                continue
            direct, _all = mapslim(go, godag, slimdag)
            slims |= direct
        for s in slims:
            name = slimdag[s].name if s in slimdag else s
            term_members.setdefault(name, set()).add(gene)
    return term_members


@register
class GoslimBackend(EnrichmentBackend):
    name = "goslim"

    def run(self, fm: FeatureMatrix, config: EnrichmentConfig) -> pd.DataFrame:
        universe = set(fm.features.astype(str))
        if config.obo_path and config.slim_obo_path and config.gene2go_path:
            term_members = _build_obo_membership(config, universe)
        elif config.gmt_path:
            gs = parse_gmt(config.gmt_path)
            term_members = {t: set(m) for t, m in gs.items()}
        else:
            raise ValueError("GoSlim needs either enrichment.gmt_path (a GO-slim "
                             "gmt) or obo_path+slim_obo_path+gene2go_path")

        per_sample: dict[str, pd.Series] = {}
        for sample in fm.samples:
            genes = select_genes(fm.data[sample], config)
            if not genes:
                continue
            ser = _score_from_membership(genes, term_members, config.goslim_score)
            if len(ser):
                per_sample[sample] = ser

        if not per_sample:
            raise RuntimeError("GoSlim produced no results for any sample")
        mat = pd.DataFrame(per_sample).reindex(columns=[c for c in fm.samples
                                                        if c in per_sample])
        return mat.fillna(0.0)
