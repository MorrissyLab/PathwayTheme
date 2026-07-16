"""Enrichment-backend interface, registry, caching, and dispatcher.

Every backend consumes a ``FeatureMatrix`` and returns a ``ScoreMatrix`` with
the identical ``terms x samples`` shape, so nothing downstream branches on the
backend.  Results are cached to a TSV keyed on backend + geneset, mirroring the
notebook's "reuse saved results" behaviour.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

from ..config import EnrichmentConfig
from ..contracts import FeatureMatrix, ScoreMatrix


class EnrichmentBackend(ABC):
    """Base class for enrichment backends."""

    name: str = "base"

    @abstractmethod
    def run(self, fm: FeatureMatrix, config: EnrichmentConfig) -> pd.DataFrame:
        """Return a ``terms x samples`` score DataFrame (columns align to fm.samples)."""
        ...


_REGISTRY: dict[str, type[EnrichmentBackend]] = {}


def register(cls: type[EnrichmentBackend]) -> type[EnrichmentBackend]:
    _REGISTRY[cls.name] = cls
    return cls


def get_backend(name: str) -> EnrichmentBackend:
    key = name.lower()
    if key not in _REGISTRY:
        # trigger lazy imports of built-in backends (tolerate optional deps)
        for mod in ("ssgsea", "enrichr", "goslim"):
            try:
                __import__(f"{__package__}.{mod}")
            except Exception:
                pass
    if key not in _REGISTRY:
        raise ValueError(f"Unknown enrichment backend {name!r}; "
                         f"available: {sorted(_REGISTRY)}")
    return _REGISTRY[key]()


def _cache_path(config: EnrichmentConfig) -> Path | None:
    if not config.cache_dir:
        return None
    return Path(config.cache_dir) / f"{config.geneset}_{config.backend}_scores.tsv"


def run_enrichment(fm: FeatureMatrix, config: EnrichmentConfig) -> ScoreMatrix:
    """Run (or load from cache) the configured backend -> ScoreMatrix."""
    cache = _cache_path(config)
    if cache is not None and cache.exists():
        data = pd.read_csv(cache, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
        data.columns = data.columns.astype(str)
        # align metadata to the cached columns
        cols = [c for c in data.columns if c in fm.metadata.table.index]
        data = data.loc[:, cols]
        return ScoreMatrix(data, fm.metadata.align_to(cols),
                           backend=config.backend, geneset=config.geneset)

    backend = get_backend(config.backend)
    data = backend.run(fm, config)
    data.columns = data.columns.astype(str)
    data.index = data.index.astype(str)
    # keep only real samples, drop all-NaN rows/cols
    data = data.dropna(axis=0, how="all").dropna(axis=1, how="all")
    cols = [c for c in data.columns if c in fm.metadata.table.index]
    data = data.loc[:, cols]

    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(cache, sep="\t")

    return ScoreMatrix(data, fm.metadata.align_to(cols),
                       backend=config.backend, geneset=config.geneset)
