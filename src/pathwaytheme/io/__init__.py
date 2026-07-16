"""Input adapters: raw omic / matrices -> FeatureMatrix + SampleMetadata."""

from .base import InputAdapter, load_input
from .genesets import parse_gmt, gene_universe

__all__ = ["InputAdapter", "load_input", "parse_gmt", "gene_universe"]
