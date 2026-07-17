# Differential pathway analysis + category roll-up

**Input:** an enrichment score matrix + a grouping (which samples are which
group).
**Output:** a per-pathway comparison table (effect, p-value, FDR) for each
contrast, optionally summarized into broad categories.

This is the statistical branch that runs *alongside* PCA on the same
`ScoreMatrix` + `Grouping`. Where PCA is unsupervised (finds axes of variation),
the differential step is supervised: it asks *which pathways differ between two
groups you name*. It is the Python analogue of the `limma` /
`wilcox_pathway_stats` step in the reference ASPS R workflow.

---

## 1. What gets compared

Given group labels (from a metadata column or a `Grouping`), the step builds one
or more two-group contrasts:

| you specify | contrasts run |
|---|---|
| nothing | **one-vs-rest** for every group (`hi_vs_rest`, `lo_vs_rest`, …) |
| `reference="lo"` | every other group vs the reference (`hi_vs_lo`, …) |
| `contrasts=[["hi","lo"]]` | exactly those pairs |

For each contrast the score matrix is split into a *case* column set and a
*reference* column set. Groups smaller than `min_group_size` (default 2) are
skipped.

---

## 2. The three tests

All three take the pathway × sample scores for the two groups and return, per
pathway, an **effect** (= mean_case − mean_reference), a **statistic**, and a
**p-value**. FDR is Benjamini–Hochberg across pathways *within each contrast*.

### `welch` — Welch's t-test (default)
Unequal-variance two-sample t (`scipy.stats.ttest_ind(..., equal_var=False)`).
Parametric, robust to unequal group sizes/variances, no cross-pathway
assumptions. A good default for continuous ssGSEA scores.

### `mannwhitney` — Mann–Whitney U
Non-parametric rank test (`scipy.stats.mannwhitneyu`, two-sided). Matches the
`wilcox_pathway_stats` helper in the reference workflow. Use when scores are
skewed or you prefer a distribution-free test.

### `moderated_t` — limma empirical-Bayes moderated t
A faithful port of `limma`'s two-group moderated t-test (Smyth, 2004), for
closest parity with the ASPS notebook. Python has no `limma`, so the variance
model is reimplemented directly:

1. **Per-pathway pooled model.** For pathway *i* with group means
   $\bar{x}_{i}$, $\bar{y}_{i}$ and sizes $n_1, n_2$:

   $$
   \beta_i=\bar{x}_i-\bar{y}_i,\qquad
   s_i^2=\frac{\sum(x-\bar{x}_i)^2+\sum(y-\bar{y}_i)^2}{n_1+n_2-2},\qquad
   d_i=n_1+n_2-2,
   $$

   with the (pathway-independent) unscaled SD
   $v=\sqrt{1/n_1+1/n_2}$.

2. **Fit the variance prior** (`fit_fdist`, a port of `limma::fitFDist`): the
   residual variances $s_i^2$ are fit to a scaled-$F$ distribution by the method
   of moments on $\log s_i^2$, yielding prior degrees of freedom $d_0$ and prior
   variance $s_0^2$. The prior-df solve uses a port of
   `limma::trigammaInverse` (Newton iteration on the trigamma function).

3. **Shrink each variance** toward the prior:

   $$
   \tilde{s}_i^2=\frac{d_0\,s_0^2+d_i\,s_i^2}{d_0+d_i}.
   $$

4. **Moderated t** on $d_0+d_i$ degrees of freedom:

   $$
   t_i=\frac{\beta_i}{v\,\sqrt{\tilde{s}_i^2}},\qquad
   p_i=2\,\Pr\!\big(T_{d_0+d_i}>|t_i|\big).
   $$

   When the prior is degenerate ($d_0\to\infty$, i.e. the residual variances are
   near-constant), $\tilde{s}_i^2\to s_0^2$ and the null uses the normal
   distribution.

Shrinking each pathway's variance toward the pooled prior stabilises small
groups — a pathway with a spuriously tiny within-group variance no longer
produces a huge t-statistic. This is why moderated t is preferred over a plain
t-test when group sizes are small (the usual ASPS case).

> **Verification.** `fit_fdist` recovers a known $(d_0,s_0^2)$ from simulated
> scaled-$F$ variances, `trigamma_inverse` inverts trigamma to 1e-6, and the
> moderated statistic correlates >0.95 with the ordinary pooled t while
> shrinking extreme values — see `tests/test_diff.py`.

---

## 3. The output table

One long table, all contrasts stacked. Columns:

| column | meaning |
|---|---|
| `comparison` | e.g. `hi_vs_lo`, `ASPS_vs_rest` |
| `pathway` | the term |
| `method` | `welch` / `mannwhitney` / `moderated_t` |
| `n_case`, `n_reference` | group sizes |
| `mean_case`, `mean_reference` | group mean scores |
| `effect` | `mean_case − mean_reference` |
| `statistic` | t / U / moderated-t |
| `p_value`, `fdr` | raw and BH-adjusted |
| `direction` | `up` (effect ≥ 0) / `down` |

Figures per contrast: a **volcano** (effect vs −log10 FDR) and a **top-pathway
bar** panel (strongest hits by FDR).

---

## 4. Category roll-up (downstream)

Optionally, the pathway-level table is summarized into broad **categories** — a
term → category mapping (a TSV, or a GO-ID → GO-slim map via goatools). For each
`(comparison, category)` it reports the mean/median/sum of the effect, the
number of pathways, and how many go up vs down.

> **Note on GoSlim.** This is the *downstream* use of GO-slim: collapsing GO:BP
> **results** into broad biological themes — distinct from PathwayTheme's GoSlim
> *enrichment backend*, which uses GO-slim *upstream* to build the score matrix
> that enters PCA/diff. Both are supported; they answer different questions.

---

## 5. In one line

$$
\underbrace{S}_{\text{pathways}\times\text{samples}}\ +\ \text{groups}
\ \xrightarrow[\text{two-group test}]{}\ (\text{effect},\,p,\,\text{FDR})\text{ per pathway per contrast}
\ \xrightarrow[\text{map}]{}\ \text{category summary.}
$$
