"""Write differential-analysis and category-summary tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..contracts import DiffResult


def write_diff_table(result: DiffResult, out_dir: str | Path, prefix: str) -> list[Path]:
    """Write the long differential table (all comparisons stacked)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{prefix}_differential_{result.method}.tsv"
    result.table.to_csv(p, sep="\t", index=False)
    return [p]


def write_category_table(summary: pd.DataFrame, out_dir: str | Path,
                         prefix: str) -> list[Path]:
    """Write a category roll-up summary table."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{prefix}_category_summary.tsv"
    summary.to_csv(p, sep="\t", index=False)
    return [p]
