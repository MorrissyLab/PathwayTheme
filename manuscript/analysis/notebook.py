"""Small helpers the notebooks share."""
from __future__ import annotations

import os

import pandas as pd
from IPython.display import Image, Markdown, display

from analysis import paths


def setup() -> None:
    """Work from the case-study root, where the configuration paths are relative
    to, and show wide tables in full."""
    os.chdir(paths.ROOT)
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.width", 200)
    print(f"working in {paths.ROOT.name}/")


def show_figure(stem: str, width: int = 900) -> None:
    """Display a figure as the PNG written to ``submission/figures``, so what is
    shown is exactly the file that ships with the paper."""
    hits = sorted(paths.SFIG.glob(f"{stem}*.png"))
    if not hits:
        raise FileNotFoundError(f"no submission/figures/{stem}*.png")
    for png in hits:
        display(Markdown(f"`{png.relative_to(paths.ROOT).as_posix()}`"))
        display(Image(filename=str(png), width=width))


def show_table(df: pd.DataFrame, rows: int = 10, name: str | None = None) -> None:
    """Display the head of a table with its size."""
    if name:
        display(Markdown(f"**{name}** — {len(df):,} rows × {df.shape[1]} columns"))
    display(df.head(rows))
