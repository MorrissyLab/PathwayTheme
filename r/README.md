# PathwayTheme for R

The R implementation of [PathwayTheme](../README.md). Same pipeline, same
stages, same output files:

**enrichment (ssGSEA / over-representation / GoSlim) → grouping → PCA (+ differential) → figures + tables.**

It is a native R package — no Python, no `reticulate`. It reads and writes the
same YAML configs, the same `.gmt` files and the same TSVs as the Python
package, and it is checked against the Python implementation numerically (see
[Parity](#parity)).

## Install

```r
# install.packages("remotes")
remotes::install_local("r")          # from the repository root
```

Dependencies are CRAN-only: `ggplot2`, `ggrepel`, `gridExtra`, `scales`, `yaml`,
plus base `grid`/`stats`. Optional extras: `hdf5r` for `.h5ad` input, `enrichR`
for the online Enrichr service, `limma` and `testthat` for the tests.

## Usage

### Run the whole thing in one call

```r
library(pathwaytheme)

pt_run_pipeline("my_config.yaml")

# ...or without a file:
pt_run_pipeline(
  "input.matrix_path"    = "expr.tsv",
  "input.metadata_path"  = "meta.tsv",
  "enrichment.backend"   = "ssgsea",
  "enrichment.gmt_path"  = "go_bp.gmt",
  "grouping.target_col"  = "treatment",
  output_dir             = "results")
```

The YAML configs in [`examples/`](../examples) drive either implementation
unchanged — the section names, field names and defaults are identical.

From the command line:

```sh
Rscript -e 'pathwaytheme::pt_main()' init -o my_config.yaml
Rscript -e 'pathwaytheme::pt_main()' run my_config.yaml
```

### Or run it step by step

Each stage is one function, and its output feeds the next:

```r
data    <- pt_load_matrix("expr.tsv", metadata = "meta.tsv")   # or pt_load_h5ad("h5ad_dir/")
data    <- pt_filter_samples(data, "tumor_type", "ASPS")       # optional: keep one group
scores  <- pt_enrich(data, backend = "ssgsea", gmt = "go_bp.gmt", geneset = "GO_BP")
groups  <- pt_group(scores, mode = "target", target_col = "treatment")
results <- pt_pca(scores, groups)
pt_figures(results, "results", geneset = "GO_BP")
```

Test **which pathways differ between groups**, then roll the hits up into broad
categories:

```r
d <- pt_diff(scores, groups, method = "moderated_t")   # welch | mannwhitney | moderated_t
pt_significant(d, 0.05)                                # per-pathway effect, p, FDR

summary <- pt_summarize_categories(d, "gobp_categories.tsv",
                                   significance_col = "fdr")

pt_sanity_heatmap(scores, "qc.pdf", label_col = "treatment")   # full-matrix QC heatmap
```

Change the method by changing one line — the rest stays the same:

```r
scores <- pt_enrich(data, backend = "enrichr", gmt = "go_bp.gmt", top_n = 200)
scores <- pt_enrich(data, backend = "goslim",  gmt = "goslim.gmt")

groups <- pt_group(scores, mode = "auto", n_clusters = 4)
groups <- pt_group(scores, mode = "target", target_col = "cell_type",
                   scope_col = "sample_id")            # one PCA per sample
```

### Look at the numbers yourself

```r
r <- pt_pca(scores, groups)[[1]]
r$scores               # samples positioned on each principal component
r$loadings             # how much each pathway contributes to each component
r$variance_explained   # how much variation each component captures
r$signatures           # the pathways that define each group
```

## What you get

Identical to the Python package: **11 figures + 3 tables** per PCA scope under
`results/<geneset>/sample_pca/<scope>/`, a gene-set coverage table, and — when
enabled — the differential table, a volcano and top-pathway figure per
comparison, a category summary, and the QC sanity heatmap.

## Name mapping

R function names are prefixed `pt_` (several of the Python names — `diff`,
`pca`, `figures` — collide with base R).

| Python | R |
|---|---|
| `pt.load_matrix` / `pt.load_h5ad` | `pt_load_matrix()` / `pt_load_h5ad()` |
| `pt.filter_samples` | `pt_filter_samples()` |
| `pt.enrich` | `pt_enrich()` |
| `pt.group` | `pt_group()` |
| `pt.pca` | `pt_pca()` |
| `pt.diff` | `pt_diff()` |
| `pt.summarize_categories` | `pt_summarize_categories()` |
| `pt.sanity_heatmap` | `pt_sanity_heatmap()` |
| `pt.figures` | `pt_figures()` |
| `pt.run_pipeline` | `pt_run_pipeline()` |
| `DiffResult.significant(a)` | `pt_significant(d, a)` |
| `PCAResult.scores` etc. | `r$scores` etc. |

## Parity

The R implementation does not wrap the Python one, so every stage is checked
against it numerically on shared inputs. To reproduce:

```sh
.venv/Scripts/python.exe r/parity/run_python.py   # write reference outputs
Rscript r/parity/compare.R                        # run R, compare, report
```

All 38 checks agree to machine precision (largest observed absolute difference
`2e-13`, in the PCA scores): ssGSEA scores, over-representation, GO-slim,
gene-set coverage, PCA variance/scores/loadings/signatures, all three
differential methods with all three contrast forms, the limma variance model,
BH correction, and the category roll-up. The report is written to
`r/parity/r/parity_report.tsv`.

Three details had to be matched deliberately, and are the things to re-check if
the two ever drift:

- **ssGSEA** is the ssGSEAProjection area-under-the-running-sum score, with
  `NES = ES / (max(ES) - min(ES))` taken over the whole result table. The area
  is evaluated in closed form (`sum_i RES_i = sum_j (N - j + 1) * t_j`) so cost
  scales with gene-set size rather than with the number of features.
- **Rank ties.** NumPy and pandas sort descending by reversing a *stable
  ascending* sort, so tied values come out in reverse of their original order.
  `order(decreasing = TRUE)` does the opposite. On a matrix that is ~45% exact
  zeros this moves ssGSEA scores by ~0.4 — the package uses the NumPy
  convention throughout (`.pt_order_desc()`).
- **PCA signs.** Components are flipped so each one's largest-magnitude loading
  is positive, matching scikit-learn's `svd_flip(..., u_based_decision = FALSE)`.
  Without this, scores and loadings would differ by an arbitrary sign per PC.

### Known differences

- **Auto-clustering** (`pt_group(mode = "auto")`) uses R's `kmeans` and
  `hclust(method = "ward.D2")` rather than scikit-learn's. Cluster *identities*
  are arbitrary in both, and the chosen partition can differ; the labelled
  modes (`target`, `existing`) are exact.
- **Figures** carry the same content, panels and colours, but are drawn with
  grid/ggplot2 rather than matplotlib, so they are not pixel-identical.
- **`.h5ad` input** is read with `hdf5r`, covering dense, CSR and CSC `X` and
  both the modern and legacy `obs` layouts. Exotic AnnData encodings are not
  supported; export to a matrix in that case.
- **GO-slim in OBO mode** uses this package's own OBO reader and `pt_mapslim()`
  rather than `goatools`.

## Development

```r
roxygen2::roxygenise("r")                                  # regenerate NAMESPACE + man
testthat::test_dir("r/tests/testthat", package = "pathwaytheme")
```
