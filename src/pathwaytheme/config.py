"""Configuration objects for the pipeline.

Everything the pipeline needs is expressed as nested dataclasses so a run is
fully described by one ``PipelineConfig`` (loadable from / dumpable to YAML).
Defaults reproduce the MOH_SM ssGSEA -> PCA settings.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict, is_dataclass
from pathlib import Path
from typing import Any, Optional

# ── PCA constants (ported verbatim from _pca_rollout.py) ──────────────────
K_MAX = 10                    # upper bound on PCs to keep / display
TOP_PCS_PER_CLUSTER = 3       # PCs combined into a cluster's signature
TOP_PATHWAYS_PER_CLUST = 25   # top pathways per cluster (signature panels)
TOP_PATHWAYS_PER_PC = 50      # pathways tabulated per PC in the TSV

# Blue-white-red, matching R colorRamp2(c(-2,0,2), c("blue","white","red"))
R_BWR_COLORS = ["#0000FF", "#FFFFFF", "#FF0000"]


@dataclass
class InputConfig:
    """How to load the omic data into a FeatureMatrix + SampleMetadata."""

    kind: str = "matrix"                 # "matrix" | "h5ad"
    matrix_path: Optional[str] = None    # features x samples table (csv/tsv/parquet)
    metadata_path: Optional[str] = None  # samples x attributes table
    sep: str = "\t"
    # --- h5ad adapter options ---
    h5ad_dir: Optional[str] = None       # directory of .h5ad files
    h5ad_paths: list[str] = field(default_factory=list)
    cluster_col: str = "seurat_clusters"
    sample_col: str = "Sample_name"
    aggregate: str = "mean"              # per-cluster pseudobulk aggregation


@dataclass
class EnrichmentConfig:
    """Which backend + which gene sets."""

    backend: str = "ssgsea"              # "ssgsea" | "enrichr" | "goslim"
    geneset: str = "GO_BP"               # label used for output naming
    gmt_path: Optional[str] = None       # .gmt file (ssgsea / goslim local mode)
    min_gene_set_size: int = 5
    max_gene_set_size: int = 5000
    threads: int = 8
    cache_dir: Optional[str] = None      # reuse computed score matrix if present
    # --- ssgsea ---
    weight: float = 0.25
    # --- enrichr ---
    enrichr_library: Optional[str] = None      # e.g. "GO_Biological_Process_2023"
    gene_selection: str = "top_n"              # "top_n" | "threshold"
    top_n: int = 200                           # genes per sample for enrichr
    threshold: Optional[float] = None
    enrichr_score: str = "-log10_padj"         # "-log10_padj" | "combined_score"
    organism: str = "human"
    # --- goslim ---
    obo_path: Optional[str] = None             # go-basic.obo
    slim_obo_path: Optional[str] = None        # goslim_generic.obo
    gene2go_path: Optional[str] = None         # NCBI gene2go
    id2sym_path: Optional[str] = None          # optional symbol<->GeneID map
    goslim_score: str = "fraction"             # "fraction" | "count"


@dataclass
class GroupingConfig:
    """How the comparison labels (and optional PCA scope) are derived."""

    mode: str = "target"                 # "target" | "existing" | "auto"
    target_col: Optional[str] = None     # metadata column for "target"/"existing"
    scope_col: Optional[str] = None      # partition into independent PCA runs
    # secondary colouring track (optional, e.g. cell type alongside class)
    secondary_col: Optional[str] = None
    size_col: Optional[str] = None       # metadata column used for dot sizing (n_cells)
    display_col: Optional[str] = None    # metadata column for observation display labels
    # --- auto-clustering ---
    method: str = "kmeans"               # "kmeans" | "hierarchical"
    n_clusters: Optional[int] = None     # fixed k; None -> auto via silhouette
    k_min: int = 2
    k_max: int = 10


@dataclass
class PCAConfig:
    k_max: int = K_MAX
    top_pcs_per_cluster: int = TOP_PCS_PER_CLUSTER
    top_pathways_per_cluster: int = TOP_PATHWAYS_PER_CLUST
    top_pathways_per_pc: int = TOP_PATHWAYS_PER_PC
    min_observations: int = 3            # need >=3 columns for a meaningful PCA
    min_features: int = 10
    random_state: int = 42


@dataclass
class VizConfig:
    make_figures: bool = True
    make_tables: bool = True
    dpi: int = 120
    formats: list[str] = field(default_factory=lambda: ["pdf"])


@dataclass
class PipelineConfig:
    input: InputConfig = field(default_factory=InputConfig)
    enrichment: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    grouping: GroupingConfig = field(default_factory=GroupingConfig)
    pca: PCAConfig = field(default_factory=PCAConfig)
    viz: VizConfig = field(default_factory=VizConfig)
    output_dir: str = "pathwaytheme_output"

    # ── serialization ────────────────────────────────────────────────────
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_yaml(self, path: str | Path) -> None:
        import yaml
        Path(path).write_text(yaml.safe_dump(self.to_dict(), sort_keys=False),
                              encoding="utf-8")

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "PipelineConfig":
        return cls(
            input=InputConfig(**(d.get("input") or {})),
            enrichment=EnrichmentConfig(**(d.get("enrichment") or {})),
            grouping=GroupingConfig(**(d.get("grouping") or {})),
            pca=PCAConfig(**(d.get("pca") or {})),
            viz=VizConfig(**(d.get("viz") or {})),
            output_dir=d.get("output_dir", "pathwaytheme_output"),
        )

    @classmethod
    def from_yaml(cls, path: str | Path) -> "PipelineConfig":
        import yaml
        d = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
        return cls.from_dict(d)


def _ensure_dataclass(obj):  # pragma: no cover - defensive helper
    assert is_dataclass(obj)
    return obj
