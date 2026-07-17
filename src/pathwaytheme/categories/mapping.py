"""Load a term -> broad-category mapping used for downstream roll-up.

Two sources are supported:

* :func:`load_category_map` — a plain TSV/CSV with a term column and a category
  column (the general, backend-agnostic path).
* :func:`goslim_category_map` — map GO **IDs** to their GO-slim ancestor terms
  via goatools (optional; only when the pathway names are GO IDs).  This mirrors
  the GO:BP -> GO-slim summarization used downstream in the reference ASPS
  workflow, as opposed to the upstream GoSlim *enrichment* backend.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


def load_category_map(path: str | Path, *, key_col: str = "pathway",
                      category_col: str = "category",
                      sep: Optional[str] = None) -> pd.Series:
    """Read a term -> category table into a Series indexed by term."""
    path = Path(path)
    if sep is None:
        sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    for col in (key_col, category_col):
        if col not in df.columns:
            raise KeyError(f"column {col!r} not in {path.name}; "
                           f"found {list(df.columns)}")
    mapping = df.set_index(df[key_col].astype(str))[category_col].astype(str)
    return mapping[~mapping.index.duplicated(keep="first")]


def goslim_category_map(go_ids, obo_path: str | Path, slim_obo_path: str | Path,
                        ) -> pd.Series:
    """Map GO IDs to a representative GO-slim term (requires goatools).

    Returns a Series ``go_id -> slim term name``.  Terms with no slim ancestor
    are omitted (callers fall back to the unmapped label).
    """
    try:
        from goatools.obo_parser import GODag
        from goatools.mapslim import mapslim
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("goslim_category_map needs the [goslim] extra "
                          "(pip install 'pathwaytheme[goslim]')") from exc

    godag = GODag(str(obo_path))
    slimdag = GODag(str(slim_obo_path))
    out: dict[str, str] = {}
    for gid in map(str, go_ids):
        if gid not in godag:
            continue
        _direct, allslim = mapslim(gid, godag, slimdag)
        if not allslim:
            continue
        # pick the slim term with the shortest name as the broad representative
        best = sorted(allslim, key=lambda s: (slimdag[s].depth, s))[0]
        out[gid] = slimdag[best].name
    return pd.Series(out, dtype=object)
