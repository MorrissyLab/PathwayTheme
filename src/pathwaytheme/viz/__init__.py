"""Figure + table generation for PCA pathway analysis."""

from .pca_figures import render_pca_figures
from .tables import write_pca_tables

__all__ = ["render_pca_figures", "write_pca_tables"]
