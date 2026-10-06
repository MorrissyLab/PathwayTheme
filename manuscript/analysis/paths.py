"""Where the case study reads its inputs and writes its outputs.

Every location is relative to one root, the ``manuscript`` directory.  Set
the environment variable ``MANUSCRIPT_ROOT`` to work on a copy of it instead --
for example to regenerate everything into a fresh directory and compare the
result with the shipped ``submission/``.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("MANUSCRIPT_ROOT")
            or os.environ.get("SARCOMA_ROOT")      # the directory's former name
            or Path(__file__).resolve().parents[1])

# inputs
RAW = ROOT / "data" / "raw"
DELIVERED = RAW / "delivered"
BROHL = RAW / "brohl_supplement"
GDC = RAW / "gdc_clinical"
GENESETS = ROOT / "data" / "genesets"
CONFIGS = ROOT / "configs"

# intermediate
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
PIPE = RESULTS / "pipeline"          # one directory per `pathwaytheme run`
SCORES = RESULTS / "scores"          # the cached ssGSEA score matrices
TABLES = RESULTS / "tables"
MYCN = ROOT / "nbl_mycn"             # inputs and runs of the MYCN analysis

# outputs: what the paper shows
SUB = ROOT / "submission"
SFIG = SUB / "figures"
STAB = SUB / "tables"


def pca_dir(run: str, geneset: str) -> Path:
    """The decomposition a pipeline run wrote for one gene-set collection."""
    return PIPE / run / geneset / "sample_pca" / "__all__"
