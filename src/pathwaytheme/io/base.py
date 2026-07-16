"""Input-adapter interface + dispatcher.

An adapter turns a source (a matrix file, or a directory of .h5ad) into a
``(FeatureMatrix)`` — features x samples plus per-sample metadata.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from ..config import InputConfig
from ..contracts import FeatureMatrix


class InputAdapter(ABC):
    """Base class for all input adapters."""

    @abstractmethod
    def load(self, config: InputConfig) -> FeatureMatrix:
        ...


def load_input(config: InputConfig) -> FeatureMatrix:
    """Dispatch to the adapter named by ``config.kind``."""
    kind = config.kind.lower()
    if kind == "matrix":
        from .matrix import MatrixAdapter
        return MatrixAdapter().load(config)
    if kind == "h5ad":
        from .h5ad import H5adAdapter
        return H5adAdapter().load(config)
    raise ValueError(f"Unknown input kind {config.kind!r}; expected 'matrix' or 'h5ad'")
