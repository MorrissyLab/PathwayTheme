"""Shared matplotlib defaults for the figures drawn outside the package."""
from __future__ import annotations

from functools import wraps

import matplotlib as mpl

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e6e5e1"
MUTED = "#9a9992"


def apply_style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,           # editable text in the PDFs
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelcolor": INK,
        "axes.edgecolor": INK_2,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.fontsize": 6,
        "legend.frameon": False,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.2,
        "text.color": INK,
    })


def styled(draw):
    """Draw a figure under :func:`apply_style` without leaking the settings to
    figures drawn later in the same session."""
    @wraps(draw)
    def wrapper(*args, **kwargs):
        with mpl.rc_context():
            apply_style()
            return draw(*args, **kwargs)
    return wrapper
