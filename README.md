# PathwayTheme

Modular pipeline: **any omic data → enrichment (ssGSEA / EnrichR / GoSlim) →
grouping (target column / existing labels / auto-clustering) → PCA pathway
analysis → summary figures + tables.**

See [PLAN.md](PLAN.md) for the design rationale.

## Workflow

Four stages, each connected by one plain-DataFrame contract, so any stage can
be swapped without touching the others:

```
 raw omic / matrix          FeatureMatrix        ScoreMatrix          Grouping        PCAResult        figures
 ┌───────────────┐          (features×samples)   (terms×samples)      (sample→label)  (per scope)      + tables
 │ .h5ad         │  ┌─────┐ ───────────────────▶ ┌────────────┐ ────▶ ┌──────────┐ ──▶ ┌─────┐ ───────▶ 11 PDFs
 │ matrix (tsv…) │─▶│ io  │                       │ enrichment │       │ grouping │     │ PCA │           3 TSVs
 └───────────────┘  └─────┘                       └────────────┘       └──────────┘     └─────┘
                     load      ssgsea/enrichr/goslim   target/existing/auto   z-score→PCA→signatures
```

1. **io** — load a single-cell `.h5ad` (→ per-cluster pseudobulk) or a
   pre-quantified `features × samples` matrix.
2. **enrichment** — score every sample against gene sets. All three backends
   return the identical `terms × samples` matrix:
   - `ssgsea` — continuous per-sample NES (gseapy).
   - `enrichr` — over-representation of each sample's gene list (`-log10 padj`).
   - `goslim` — coarse GO-slim category scores.
3. **grouping** — decide what to compare/colour by:
   - `target` — a metadata column (e.g. `classification`, `treatment`).
   - `existing` — cluster labels already in the data.
   - `auto` — cluster the samples (KMeans/hierarchical, silhouette-picked *k*).
   - optional **scope** column → one independent PCA per group (e.g. per sample).
4. **pca** — z-score features → PCA → per-observation pathway signatures.
5. **viz** — 11 figures (scores heatmap, biplots, pairwise grid, signatures,
   PC-loading bars/heatmaps, union views) + 3 TSVs, per scope.

## Install (uv)

```bash
uv venv --python 3.12
uv pip install -e ".[all]"       # h5ad input + GoSlim + nicer figure labels
uv pip install -e ".[all,dev]"   # + pytest for the test suite
```

Extras: `h5ad` (single-cell input), `goslim` (goatools), `viz`
(adjustText + distinctipy), `dev` (pytest). Core install (`uv pip install -e .`)
covers ssGSEA + EnrichR from a matrix. Run the tests with `uv run pytest -q`.

## Usage

### Option A — run the whole pipeline (one call)

From the command line with a YAML config:

```bash
uv run pathwaytheme run examples/moh_sm.yaml
uv run pathwaytheme init -o my_config.yaml     # write a template config to edit
```

Or from Python — pass a YAML path, or set fields inline:

```python
import pathwaytheme as pt

pt.run_pipeline("examples/moh_sm.yaml")

# ...or configure inline (dotted keys), no YAML needed:
pt.run_pipeline(**{
    "input.kind": "matrix", "input.matrix_path": "expr.tsv",
    "input.metadata_path": "meta.tsv",
    "output_dir": "results",
})
```

### Option B — call each stage yourself

Each stage is a single function; results flow straight into the next one:

```python
import pathwaytheme as pt

fm      = pt.load_matrix("expr.tsv", metadata="meta.tsv")   # or pt.load_h5ad("h5ad_dir/")
scores  = pt.enrich(fm, backend="ssgsea", gmt="go_bp.gmt", geneset="GO_BP")
groups  = pt.group(scores, mode="target", target_col="grp", scope_col="sample_id")
results = pt.pca(scores, groups, size_col="n_cells")
pt.figures(results, "results", geneset="GO_BP")
```

Swap one line to change method — nothing downstream changes:

```python
# EnrichR over-representation instead of ssGSEA:
scores = pt.enrich(fm, backend="enrichr", gmt="go_bp.gmt", top_n=200)

# GoSlim categories:
scores = pt.enrich(fm, backend="goslim", gmt="goslim.gmt", goslim_score="fraction")

# auto-cluster the samples instead of using a metadata column:
groups = pt.group(scores, mode="auto", n_clusters=4)   # or omit n_clusters to auto-pick k
```

Each `pt.*` wrapper accepts the extra options of its config
(`pathwaytheme.config`); drop down to the stage modules
(`pathwaytheme.io`, `.enrichment`, `.grouping`, `.pca`, `.viz`) for full control.

### Inspecting results in code

```python
results = pt.pca(scores, groups)
r = results[0]
r.scores          # observations × PCs
r.loadings        # PCs × pathways
r.variance_explained
r.signatures      # long table: cluster, pathway, signature_score, direction
```

## Outputs

Per PCA scope, written to `<output_dir>/<geneset>/sample_pca/<scope>/`:

| Kind | Files |
|---|---|
| Figures (11 PDF) | scores heatmap · PC1-PC2 biplot ×2 · pairwise grid ×2 · per-cluster signatures · PC-loading bars · pathway×PC heatmap + clustermap · union clustermap + signature |
| Tables (3 TSV) | `variance_explained` · `top_pathways_per_pc` · `per_cluster_signatures` |

## Pipeline stages (modules)

| Stage | Module | Contract out |
|---|---|---|
| Input | `pathwaytheme.io` | `FeatureMatrix` (features × samples) |
| Enrichment | `pathwaytheme.enrichment` | `ScoreMatrix` (terms × samples) |
| Grouping | `pathwaytheme.grouping` | `Grouping` (sample → label) |
| PCA | `pathwaytheme.pca` | `PCAResult` (per scope) |
| Viz | `pathwaytheme.viz` | PDFs + TSVs |

## Verified against MOH_SM

The ported stages were checked against the original MOH_SM outputs:

- **ssGSEA**: pooled 26-sample pseudobulk reproduces the golden NES matrix
  (Pearson r = 1.000000).
- **PCA**: variance explained, PC loadings, and per-cluster signatures are
  bit-exact vs the golden rollout TSVs (diffs ~1e-16).
- **Figures + tables**: identical 14-file output set; the 3 TSVs are
  byte-identical to golden.
