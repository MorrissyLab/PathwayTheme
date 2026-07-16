"""Deterministic colour palettes + shared drawing constants.

Ports _class_colour_lut / _celltype_colour_lut / R_BWR / _truncate from
_pca_rollout.py so the same label always maps to the same colour.
"""

from __future__ import annotations

import matplotlib.colors as mc

from ..config import R_BWR_COLORS

R_BWR = mc.LinearSegmentedColormap.from_list("R_bwr", R_BWR_COLORS)

_EXTRAS = ["#1f77b4", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf",
           "#ff7f0e", "#bcbd22", "#7f6d6d", "#3d8b8b", "#c2a44b", "#a4a4a4"]

# Known biological cell-type colours (kept for MOH_SM parity; harmless elsewhere)
_CELL_TYPE_COLORS = {
    "tumor": "#d62728", "tumor_proliferating": "#a02020", "proliferating": "#ff9896",
    "myeloid_proliferating": "#fb7858", "low_quality": "#7f7f7f", "low_numbers": "#bdbdbd",
    "fibroblast": "#8c564b", "CAF": "#a0522d", "myofibroblast": "#c49c94",
    "endothelial": "#1f77b4", "endothelial_lymphatic": "#4ba6d5", "vascular": "#155b8f",
    "pericyte": "#9edae5", "endo-peri": "#17becf", "smooth_muscle_cell": "#5b8db5",
    "macrophage": "#2ca02c", "monocyte": "#5fbb5f", "myeloid": "#1f8a1f",
    "dendritic_cell": "#8de08d", "pDC": "#b9ddb9", "mononuclear_phagocyte": "#3d9b3d",
    "T-cell": "#9467bd", "CD8-T": "#7e4faf", "B-cell": "#c5b0d5", "NK_cell": "#5e3a8a",
    "Treg": "#a986d3", "MAIT": "#dac4ec", "lymphocyte": "#6f3a93", "leukocyte": "#b48fd1",
    "skeletal_muscle": "#e377c2", "muscle_cell": "#f5b6e0", "neuronal_glial": "#ff7f0e",
    "neuron": "#ffae5d", "Schwann_cells": "#d65a05", "adipocyte": "#bcbd22",
    "erythrocyte": "#a3070b", "pneumocyte": "#6a8eb0",
}


def truncate(term: str, n: int = 80) -> str:
    t = str(term)
    return t if len(t) <= n else t[:n - 3] + "..."


def class_colour_lut(labels) -> dict:
    """Palette for the primary comparison labels; tumor red, low_quality grey."""
    base = {"tumor": "#d62728", "low_quality": "#7f7f7f", "": "#dddddd"}
    out, k = {}, 0
    for c in labels:
        if c in base:
            out[c] = base[c]
        elif c in out:
            continue
        else:
            out[c] = _EXTRAS[k % len(_EXTRAS)]; k += 1
    return out


def celltype_colour_lut(labels) -> dict:
    """Palette for a secondary label track, reusing known biological colours."""
    out, k = {}, 0
    for t in labels:
        if t in _CELL_TYPE_COLORS:
            out[t] = _CELL_TYPE_COLORS[t]
        elif not t:
            out[t] = "#dddddd"
        elif t in out:
            continue
        else:
            out[t] = _EXTRAS[k % len(_EXTRAS)]; k += 1
    return out
