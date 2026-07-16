"""PathwayTheme — modular omic-enrichment + PCA pathway-analysis pipeline.

Stages (each connected by a plain-DataFrame contract, see :mod:`pathwaytheme.contracts`):

    raw omic / matrix  ->  [io]         ->  FeatureMatrix (features x samples)
    FeatureMatrix      ->  [enrichment] ->  ScoreMatrix   (terms x samples)
    ScoreMatrix        ->  [grouping]   ->  Grouping      (sample -> label)
    ScoreMatrix        ->  [pca]        ->  PCAResult      (per scope)
    PCAResult          ->  [viz]        ->  figures + tables

The public entry point is :func:`pathwaytheme.pipeline.run`.
"""

from .contracts import FeatureMatrix, SampleMetadata, ScoreMatrix, Grouping, PCAResult
from .config import PipelineConfig
from .api import (load_matrix, load_h5ad, enrich, group, pca, figures, run_pipeline)

__version__ = "0.1.0"

__all__ = [
    # high-level API (one call per stage, or the whole pipeline)
    "load_matrix", "load_h5ad", "enrich", "group", "pca", "figures", "run_pipeline",
    # contracts + config (for full control)
    "FeatureMatrix", "SampleMetadata", "ScoreMatrix", "Grouping", "PCAResult",
    "PipelineConfig",
    "__version__",
]
