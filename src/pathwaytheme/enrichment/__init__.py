"""Enrichment backends.  Every backend returns the same ScoreMatrix shape."""

from .base import EnrichmentBackend, get_backend, run_enrichment

__all__ = ["EnrichmentBackend", "get_backend", "run_enrichment"]
