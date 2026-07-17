"""Differential pathway analysis (parallel to PCA)."""

from .core import run_diff
from .stats import benjamini_hochberg, fit_fdist, trigamma_inverse

__all__ = ["run_diff", "benjamini_hochberg", "fit_fdist", "trigamma_inverse"]
