# PathwayTheme

Turn **any omic data** into a pathway-level summary:

**enrichment (ssGSEA / EnrichR / GoSlim) → grouping → PCA → figures + tables.**

📄 **How the PCA step works (the math, step by step): [docs/METHODS.md](docs/METHODS.md).**

## What it does

```
   your data            pathway scores        groups            PCA            figures + tables
   ┌──────────┐         ┌────────────┐        ┌────────┐        ┌──────┐       ┌──────────────┐
   │ .h5ad    │  enrich │ per-sample │  group │ compare│   PCA  │ trends│  viz  │ 11 plots     │
   │ or matrix│ ───────▶│ pathway    │ ──────▶│  by ...│ ──────▶│ across│ ─────▶│  3 tables    │
   └──────────┘         │ enrichment │        └────────┘        └──────┘       └──────────────┘
```

1. **Load** a single-cell `.h5ad` (summarised per cluster) or a
   `features × samples` table (genes, proteins, etc.).
2. **Enrich** each sample against gene sets, using one of:
   - **ssGSEA** — a continuous score per pathway per sample.
   - **EnrichR** — over-representation of each sample's top genes.
   - **GoSlim** — coarse GO-slim category scores.
3. **Group** — choose what to compare / colour by:
   - a **metadata column** (e.g. `treatment`, `cell type`, `classification`),
   - **existing** cluster labels, or
   - **automatic** clustering when you have no labels.
   Optionally run **one PCA per sample** (or per any grouping).
4. **PCA + figures** — see which pathways drive the variation, and which
   pathways define each group, as ready-to-use PDFs and tables.
   → full calculation walkthrough in [docs/METHODS.md](docs/METHODS.md).

## Install

```bash
uv venv --python 3.12
uv pip install -e ".[all]"
```

(`[all]` adds single-cell `.h5ad` input, GoSlim, and nicer plot labels. Plain
`uv pip install -e .` is enough for ssGSEA/EnrichR from a matrix.)

## Usage

### Run the whole thing in one call

From the command line with a config file:

```bash
uv run pathwaytheme init -o my_config.yaml     # create a template to edit
uv run pathwaytheme run my_config.yaml
```

Or from Python:

```python
import pathwaytheme as pt

pt.run_pipeline("my_config.yaml")

# ...or without a file:
pt.run_pipeline(**{
    "input.matrix_path": "expr.tsv",
    "input.metadata_path": "meta.tsv",
    "enrichment.backend": "ssgsea",
    "enrichment.gmt_path": "go_bp.gmt",
    "grouping.target_col": "treatment",
    "output_dir": "results",
})
```

See `examples/moh_sm.yaml` (single-cell → ssGSEA) and
`examples/matrix_enrichr.yaml` (matrix → EnrichR) for full configs.

### Or run it step by step

Each step is one function, and its output feeds the next:

```python
import pathwaytheme as pt

data    = pt.load_matrix("expr.tsv", metadata="meta.tsv")   # or pt.load_h5ad("h5ad_folder/")
scores  = pt.enrich(data, backend="ssgsea", gmt="go_bp.gmt", geneset="GO_BP")
groups  = pt.group(scores, mode="target", target_col="treatment")
results = pt.pca(scores, groups)
pt.figures(results, "results", geneset="GO_BP")
```

Change the method by changing one line — the rest stays the same:

```python
scores = pt.enrich(data, backend="enrichr", gmt="go_bp.gmt", top_n=200)     # over-representation
scores = pt.enrich(data, backend="goslim",  gmt="goslim.gmt")               # GO-slim categories

groups = pt.group(scores, mode="auto", n_clusters=4)                        # cluster automatically
groups = pt.group(scores, mode="target", target_col="cell_type",
                  scope_col="sample_id")                                    # one PCA per sample
```

### Look at the numbers yourself

```python
r = pt.pca(scores, groups)[0]
r.scores               # samples positioned on each principal component
r.loadings             # how much each pathway contributes to each component
r.variance_explained   # how much variation each component captures
r.signatures           # the pathways that define each group
```

## What you get

For each group (or sample), written to
`results/<geneset>/sample_pca/<group>/`:

- **11 figures** — scores heatmap, PC1–PC2 biplots, a pairwise-PC grid,
  per-group pathway signatures, per-component top pathways, and pathway
  loading heatmaps.
- **3 tables** — variance explained, top pathways per component, and the
  per-group pathway signatures.

## Common options

| Where | Option | Meaning |
|---|---|---|
| `enrich` | `backend` | `ssgsea` · `enrichr` · `goslim` |
| `enrich` | `gmt` | path to a gene-set `.gmt` file |
| `enrich` | `top_n` / `threshold` | how EnrichR/GoSlim pick each sample's genes |
| `group` | `mode` | `target` · `existing` · `auto` |
| `group` | `target_col` | metadata column to compare/colour by |
| `group` | `scope_col` | run a separate PCA within each value of this column |
| `pca` | `size_col` | metadata column that sets dot sizes (e.g. cell counts) |
