"""Gene-set (.gmt) parsing + gene-universe helpers.

Ported from ssgsea_analysis_fixed.ipynb (``parse_gmt`` / ``gene_universe``).
"""

from __future__ import annotations

from pathlib import Path


def parse_gmt(path: str | Path) -> dict[str, list[str]]:
    """Parse a GMT file -> {term: [genes]}.

    GMT format: ``term <tab> description <tab> gene1 <tab> gene2 ...``
    """
    path = Path(path)
    gs: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            sp = line.rstrip("\n").split("\t")
            if len(sp) >= 3:
                gs[sp[0]] = [x for x in sp[2:] if x]
    if not gs:
        raise ValueError(f"No gene sets parsed from {path}")
    return gs


def gene_universe(gs: dict[str, list[str]]) -> set[str]:
    """Union of all genes across every gene set."""
    out: set[str] = set()
    for v in gs.values():
        out.update(v)
    return out


def harmonize_case(expr_index, gs: dict[str, list[str]]):
    """Decide whether upper-casing improves gene overlap (ports the ipynb logic).

    Returns ``(index_transform, gs_out)`` where ``index_transform`` is a callable
    applied to the expression index and ``gs_out`` is the (possibly upper-cased)
    gene-set dict.  If upper-casing does not improve overlap, both are identity.
    """
    universe = gene_universe(gs)
    raw_ov = len(set(map(str, expr_index)) & universe)
    up_ov = len({str(g).upper() for g in expr_index} & {g.upper() for g in universe})
    if up_ov > raw_ov:
        gs_out = {k: [x.upper() for x in v] for k, v in gs.items()}
        return (lambda s: str(s).upper()), gs_out
    return (lambda s: str(s)), gs
