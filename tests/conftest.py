"""Shared synthetic fixtures — no external data or network needed."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathwaytheme.contracts import FeatureMatrix, SampleMetadata, ScoreMatrix


@pytest.fixture
def rng():
    return np.random.default_rng(0)


@pytest.fixture
def synthetic_gmt(tmp_path, rng):
    """A tiny GO-style .gmt: 15 terms, each 8-18 genes drawn from 200 genes."""
    genes = [f"G{i:03d}" for i in range(200)]
    path = tmp_path / "synthetic.gmt"
    lines = []
    for t in range(15):
        k = int(rng.integers(8, 18))
        members = list(rng.choice(genes, size=k, replace=False))
        lines.append(f"TERM_{t:02d}\tdesc\t" + "\t".join(members))
    path.write_text("\n".join(lines), encoding="utf-8")
    return path, genes


@pytest.fixture
def feature_matrix(rng, synthetic_gmt):
    """200 genes x 12 samples, with 2 latent groups baked in for structure."""
    _, genes = synthetic_gmt
    n_samples = 12
    base = rng.normal(0, 1, size=(len(genes), n_samples))
    # inject group structure: first 6 samples get +2 on genes 0-40
    base[:40, :6] += 2.0
    cols = [f"S{i:02d}" for i in range(n_samples)]
    data = pd.DataFrame(base, index=genes, columns=cols)
    meta = pd.DataFrame({
        "sample_id": ["A"] * 6 + ["B"] * 6,
        "cluster_id": [str(i) for i in range(6)] * 2,
        "grp": ["hi"] * 6 + ["lo"] * 6,
        "n_cells": rng.integers(20, 500, size=n_samples),
    }, index=cols)
    return FeatureMatrix(data, SampleMetadata(meta))


@pytest.fixture
def score_matrix(rng):
    """40 terms x 12 samples continuous score matrix + metadata."""
    cols = [f"S{i:02d}" for i in range(12)]
    data = pd.DataFrame(rng.normal(0, 1, size=(40, 12)),
                        index=[f"TERM_{i:02d}" for i in range(40)], columns=cols)
    data.iloc[:10, :6] += 3.0  # structure
    meta = pd.DataFrame({
        "sample_id": ["A"] * 6 + ["B"] * 6,
        "cluster_id": [str(i) for i in range(6)] * 2,
        "grp": ["hi"] * 6 + ["lo"] * 6,
        "n_cells": rng.integers(20, 500, size=12),
    }, index=cols)
    return ScoreMatrix(data, SampleMetadata(meta), backend="synthetic", geneset="SYN")
