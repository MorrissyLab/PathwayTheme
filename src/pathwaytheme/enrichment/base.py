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

    #: True when a sample's scores depend only on that sample -- ssGSEA ranks
    #: genes within a sample, Enrichr takes that sample's top-n genes, GO-slim
    #: scores that sample's own profile.  Only then may the cache be filled in
    #: sample by sample; a cross-sample-normalised backend must set this False
    #: so a partial cache hit recomputes everything.
    per_sample: bool = True

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
    # The missing-value policy changes the scores, so it has to be part of the
    # key -- but only when it is non-default, so existing caches stay valid.
    suffix = ""
    if getattr(config, "missing", "as_is") != "as_is":
        suffix += f"_miss-{config.missing}"
    if getattr(config, "min_observed_fraction", 0.0) > 0:
        suffix += f"_obs{config.min_observed_fraction:g}"
    return (Path(config.cache_dir)
            / f"{config.geneset}_{config.backend}{suffix}_scores.tsv")


def _prepare_features(fm: FeatureMatrix, config: EnrichmentConfig,
                      *, verbose: bool = True) -> FeatureMatrix:
    """Resolve missing values per ``enrichment.missing`` before scoring."""
    from .coverage import apply_missing_policy
    data = apply_missing_policy(
        fm.data, getattr(config, "missing", "as_is"),
        getattr(config, "min_observed_fraction", 0.0), verbose=verbose)
    # Identity, not shape: filling missing values leaves the shape untouched,
    # so a shape comparison here would silently discard the filled matrix and
    # hand the backend the original -- letting gseapy's own zero-fill win.
    if data is fm.data:
        return fm
    return FeatureMatrix(data, fm.metadata)


def _coverage_for(fm: FeatureMatrix, config: EnrichmentConfig,
                  *, verbose: bool = True):
    """Gene-set coverage against the matrix the backend will actually see."""
    if not config.gmt_path:
        return None
    from ..io.genesets import parse_gmt
    from .coverage import geneset_coverage, warn_on_dropped
    try:
        gs = parse_gmt(config.gmt_path)
    except Exception:
        return None
    cov = geneset_coverage(fm.features, gs, config.min_gene_set_size,
                           config.max_gene_set_size)
    warn_on_dropped(cov, geneset=config.geneset, verbose=verbose)
    return cov


def _read_cache(cache: Path) -> pd.DataFrame:
    data = (pd.read_csv(cache, sep="\t", index_col=0)
            .apply(pd.to_numeric, errors="coerce"))
    data.columns = data.columns.astype(str)
    data.index = data.index.astype(str)
    return data


def run_enrichment(fm: FeatureMatrix, config: EnrichmentConfig) -> ScoreMatrix:
    """Run (or load from cache) the configured backend -> ScoreMatrix.

    The cache is keyed on backend + geneset only, so two configs that differ in
    which samples they select share one file.  A requested sample that is not
    in the cache is therefore scored and added, never silently dropped: with a
    per-sample backend a sample's scores do not depend on which other samples
    were run, so filling the cache in is exact.
    """
    cache = _cache_path(config)
    backend = get_backend(config.backend)
    wanted = [str(c) for c in fm.data.columns]
    meta = fm.metadata           # fm is narrowed below on a partial cache hit

    # Resolve missing values and account for gene-set coverage up front, so a
    # cache hit reports the same coverage a fresh run would.
    fm = _prepare_features(fm, config)
    coverage = _coverage_for(fm, config)

    cached: pd.DataFrame | None = None
    if cache is not None and cache.exists():
        cached = _read_cache(cache)
        have = [c for c in wanted if c in cached.columns]
        if len(have) == len(wanted):
            data = cached.loc[:, wanted].dropna(axis=0, how="all")
            return ScoreMatrix(data, meta.align_to(wanted),
                               backend=config.backend, geneset=config.geneset,
                               coverage=coverage)
        missing = [c for c in wanted if c not in cached.columns]
        if not getattr(backend, "per_sample", True):
            cached = None                       # must recompute as a whole
        else:
            print(f"      cache: {len(have)} of {len(wanted)} samples hit, "
                  f"scoring {len(missing)} more")
            fm = FeatureMatrix(fm.data.loc[:, missing],
                               fm.metadata.align_to(missing))

    data = backend.run(fm, config)
    data.columns = data.columns.astype(str)
    data.index = data.index.astype(str)
    # keep only real samples, drop all-NaN rows/cols
    data = data.dropna(axis=0, how="all").dropna(axis=1, how="all")
    data = data.loc[:, [c for c in data.columns if c in fm.metadata.table.index]]

    if cached is not None:
        # outer join on terms: a newly scored sample may carry a term the
        # cached samples did not, and vice versa
        data = cached.join(data, how="outer")

    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(cache, sep="\t")

    cols = [c for c in wanted if c in data.columns]
    missing = [c for c in wanted if c not in data.columns]
    if missing:
        raise RuntimeError(
            f"{len(missing)} requested sample(s) produced no scores, "
            f"e.g. {missing[:3]}")
    data = data.loc[:, cols].dropna(axis=0, how="all")
    return ScoreMatrix(data, meta.align_to(cols),
                       backend=config.backend, geneset=config.geneset,
                       coverage=coverage)
