"""Figure + table generation for PCA pathway analysis."""

from .pca_figures import render_pca_figures
from .tables import write_pca_tables
from .qc import render_sanity_heatmap
from .diff_figures import render_diff_figures, draw_volcano, draw_top_bars
from .diff_tables import write_diff_table, write_category_table
from .diagnostics import (
    draw_workflow, render_workflow, WORKFLOW_STAGES,
    draw_scree, render_scree,
    draw_pc_loadings, render_pc_loadings,
    draw_pc_pairs, render_pc_pairs,
    draw_pca_scatter, render_pca_scatter,
    draw_metadata_attribution, render_metadata_attribution,
    draw_level_attribution, render_level_attribution,
    draw_theme_heatmap, render_theme_heatmap, draw_theme_reduction,
    draw_effect_concordance, render_effect_concordance,
)

__all__ = [
    "render_pca_figures", "write_pca_tables",
    "render_sanity_heatmap",
    "render_diff_figures", "write_diff_table", "write_category_table",
    "draw_volcano", "draw_top_bars",
    "draw_workflow", "render_workflow", "WORKFLOW_STAGES",
    "draw_scree", "render_scree",
    "draw_pc_loadings", "render_pc_loadings",
    "draw_pc_pairs", "render_pc_pairs",
    "draw_pca_scatter", "render_pca_scatter",
    "draw_metadata_attribution", "render_metadata_attribution",
    "draw_level_attribution", "render_level_attribution",
    "draw_theme_heatmap", "render_theme_heatmap", "draw_theme_reduction",
    "draw_effect_concordance", "render_effect_concordance",
]
