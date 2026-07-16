# The PCA step: from enrichment + grouping to the summary

**Input:** an enrichment score matrix + a target/cluster grouping.
**Output:** the PCA summary — variance explained, PC scores, PC loadings
(top pathways per PC), and a per-cluster pathway signature.

This document covers only the PCA step (`pathwaytheme.pca` + `pathwaytheme.viz`).
How the score matrix itself is produced (ssGSEA / EnrichR / GoSlim) is not
covered here — the PCA step treats it as a given.

Notation:

- $S$ — enrichment scores, terms (pathways) $\times$ columns. $S_{ij}$ = score
  of pathway $i$ in column $j$.
- The **columns are the observations** the PCA positions — in the single-cell
  case, one column per *cluster*.
- $p$ = number of pathways (variables), $n$ = number of observations in one
  analysis.
- Defaults: $K_{\max}=10$, top PCs per cluster $=3$, top pathways per
  cluster $=25$, top pathways per PC $=50$.

---

## 0. What "enrichment + target/cluster" means as input

Two things enter the step:

1. **The score matrix $S$** — the enrichment result, pathways $\times$ columns.
2. **The grouping**, which supplies up to two roles derived from a
   target/cluster column:

   | role | what it is | how the PCA uses it |
   |---|---|---|
   | **observation** | the column identity (e.g. each **cluster**) | each column becomes one point the PCA positions |
   | **scope** *(optional)* | a column that partitions the data (e.g. `sample_id`) | run **one independent PCA per group** — e.g. one PCA per sample, over its clusters |
   | **label / target** *(optional)* | the comparison variable (e.g. cell class, treatment) | overlaid as colour on the outputs; used to interpret the axes, **not** to compute them |

Important: **PCA is unsupervised.** The target label never enters the math — it
is attached to the points afterwards so you can read the result ("PC1 separates
tumor from immune"). The cluster/column identity *is* what gets positioned; the
scope decides which columns are compared together.

For one scope we take the sub-matrix $X \in \mathbb{R}^{p\times n}$ (all
pathways $\times$ the $n$ columns in that scope), with missing values set to
$0$. A scope is analysed only if it has $n\ge 3$ observations and $\ge 10$
variable pathways.

---

## 1. Standardise each pathway (z-score across observations)

PCA finds directions of *variation*, so each pathway is put on a common scale.
For pathway $i$ across the $n$ observations:

$$
\mu_i = \frac1n\sum_{j} X_{ij},\qquad
\sigma_i = \sqrt{\frac1n\sum_{j}(X_{ij}-\mu_i)^2}\ \ (\texttt{ddof=0}),\qquad
Z_{ij} = \frac{X_{ij}-\mu_i}{\sigma_i}.
$$

Pathways that are constant across this scope ($\sigma_i\le10^{-12}$) are
dropped — they contribute no variance. After this, every retained pathway has
mean $0$ and SD $1$ across the observations.

---

## 2. Orient: observations are the points

Transpose so each observation is a row (a point in pathway space):

$$M = Z^{\top}\in\mathbb{R}^{n\times p}.$$

---

## 3. PCA

Keep

$$K=\min(n-1,\;K_{\max}),\qquad K_{\max}=10$$

components ($n$ points span at most $n-1$ dimensions; capped at 10). The SVD
$M = U\Sigma V^{\top}$ gives the three summary arrays:

**PC scores** — where each observation sits on each axis:

$$T = U\Sigma\in\mathbb{R}^{n\times K},\qquad T_{ck}=\text{position of observation }c\text{ on PC }k.$$

**PC loadings** — how much each pathway defines each axis (unit-length rows):

$$W\in\mathbb{R}^{K\times p},\qquad W_{ki}=\text{weight of pathway }i\text{ in PC }k,\qquad \sum_i W_{ki}^2=1.$$

**Variance explained** — how much of the total spread each axis captures:

$$
\lambda_k=\frac{s_k^2}{n-1},\qquad
\text{var\_explained}_k=\frac{\lambda_k}{\sum_m\lambda_m}.
$$

(The sign of a PC is arbitrary but reproducible for a fixed input.)

---

## 4. Per-cluster pathway signature

This is the main interpretive output: *for each cluster (observation), which
pathways make it distinct?* It reconstructs the part of the cluster's profile
that PCA captures, using only that cluster's dominant axes.

For observation $c$:

1. **Its top 3 PCs** $P_c$ = the PCs with the largest $|T_{ck}|$ (where $c$
   moves the most).
2. **Signature score** for pathway $i$:

$$
\boxed{\ \mathrm{sig}(c,i)=\sum_{k\in P_c} T_{ck}\,W_{ki}\ }
$$

   Each term is signed: positive when the cluster's direction on PC $k$ pulls
   pathway $i$ up, negative when it pulls it down.
3. **Rank** pathways by $|\mathrm{sig}(c,i)|$; keep the top 25, sign preserved:
   **+ = enriched (up), − = depleted (down)** in that cluster.

---

## 5. Top pathways per PC

To label the axes themselves, each PC $k$'s pathways are ranked by $|W_{ki}|$
and the top 50 are recorded with signed loading and direction. Reading the
positive vs negative ends tells you what biology the axis contrasts.

---

## 6. The summary that comes out

**Tables**
- `variance_explained` — var explained + cumulative per PC (§3).
- `top_pathways_per_pc` — top-50 pathways per PC by $|W_{ki}|$ (§5).
- `per_cluster_signatures` — top-25 signed pathways per cluster (§4).

**Figures** (the target/cluster label is overlaid as colour throughout)
- **Scores heatmap** — the $T$ matrix (observations $\times$ PCs), rows
  clustered; a colour strip shows each observation's target label.
- **PC1–PC2 biplot** — observations at $(T_{c1},T_{c2})$, coloured by target,
  dot size $\propto\log(1+\text{size})$ (e.g. cell count); the 10 pathways with
  the largest $\sqrt{W_{1i}^2+W_{2i}^2}$ drawn as loading arrows.
- **Pairwise grid** — scatter of every pair among the top-4 PCs.
- **Per-cluster signatures** — one bar panel per cluster (§4), red up / blue down.
- **PC loadings** — top ±25 loadings per PC (§5).
- **Pathway × PC heatmaps** — the 80 highest-$\max_k|W_{ki}|$ pathways, ordered
  and hierarchically clustered on their loading vectors.
- **Union views** — union of each PC's top ±10 pathways, arranged by loading
  similarity and by home-PC + sign.

---

## 7. In one line

$$
\underbrace{S}_{\text{pathways}\times\text{clusters}}\ +\ \text{grouping}
\ \xrightarrow[\text{z-score}]{}\ Z
\ \xrightarrow[\text{transpose}]{}\ M
\ \xrightarrow[\text{PCA}]{}\ (T,\,W,\,\lambda)
\ \xrightarrow[\S4\text{–}\S5]{}\ \text{signatures} + \text{PC summary} + \text{figures.}
$$

The cluster/column identities are the points; the scope decides what is
compared together; the target label colours and explains the result but does
not change it.
