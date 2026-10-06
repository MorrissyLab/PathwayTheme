"""Gene-set coverage accounting and missing-value policy.

Two things the pipeline must not do silently, both of which first showed up
when the same gene sets were run over a mass-spectrometry proteome rather than
a transcriptome:

**Sets can vanish.**  A ``.gmt`` set whose members are largely absent from the
feature matrix falls below ``min_gene_set_size`` and is dropped by the backend.
The score table then has fewer rows than the ``.gmt`` has sets, with nothing to
say which are missing or why.  That is tolerable when it never happens; it is
not tolerable when comparing two modalities, because coverage differs between
them and the two tables would silently describe different set universes.
``geneset_coverage`` records the fate of every requested set.

**Values can be missing.**  A TMT proteome is distributed as signed log-ratios
with ~20% of entries absent.  ssGSEA is rank-based, and an unstated fill puts
those entries at some arbitrary point in the ranking -- for a log-ratio matrix,
filling with 0 lands them at the median.  ``apply_missing_policy`` makes the
choice explicit, reports it, and leaves the default a no-op so existing
results are unchanged.
"""

from __future__ import annotations

import warnings

import pandas as pd

#: fates a requested gene set can meet
SCORED = "scored"
NO_OVERLAP = "dropped_no_overlap"
TOO_SMALL = "dropped_below_min_size"
TOO_LARGE = "dropped_above_max_size"

MISSING_POLICIES = ("as_is", "zero", "low")


def geneset_coverage(features, gene_sets: dict[str, list[str]],
                     min_size: int, max_size: int) -> pd.DataFrame:
    """Per requested set: how many members are measurable, and its fate.

    ``features`` is the feature index of the matrix actually handed to the
    backend.  Matching is done on the upper-cased symbol, mirroring
    ``harmonize_case``, so the count here is the count the backend will see.
    """
    present = {str(f).upper() for f in features}
    rows = []
    for name, genes in gene_sets.items():
        uniq = {str(g).upper() for g in genes}
        n_present = len(uniq & present)
        if n_present == 0:
            status = NO_OVERLAP
        elif n_present < min_size:
            status = TOO_SMALL
        elif n_present > max_size:
            status = TOO_LARGE
        else:
            status = SCORED
        rows.append({
            "pathway": name,
            "n_genes": len(uniq),
            "n_present": n_present,
            "fraction_present": (n_present / len(uniq)) if uniq else 0.0,
            "status": status,
        })
    out = pd.DataFrame(rows).sort_values(
        ["status", "pathway"], kind="stable").reset_index(drop=True)
    return out


def warn_on_dropped(coverage: pd.DataFrame, *, geneset: str = "",
                    verbose: bool = True) -> pd.DataFrame:
    """Emit a warning naming the sets that will not be scored.  Returns them."""
    dropped = coverage[coverage["status"] != SCORED]
    if dropped.empty or not verbose:
        return dropped
    label = f" ({geneset})" if geneset else ""
    lines = [f"        {r.pathway}: {r.n_present}/{r.n_genes} members measured "
             f"[{r.status}]" for r in dropped.head(10).itertuples()]
    more = "" if len(dropped) <= 10 else f"\n        ... and {len(dropped) - 10} more"
    warnings.warn(
        f"{len(dropped)} of {len(coverage)} gene sets{label} will not be "
        f"scored:\n" + "\n".join(lines) + more,
        stacklevel=2,
    )
    return dropped


def apply_missing_policy(data: pd.DataFrame, policy: str,
                         min_observed_fraction: float = 0.0,
                         *, verbose: bool = True) -> pd.DataFrame:
    """Drop under-observed features, then resolve the remaining missing values.

    ``min_observed_fraction`` drops any feature measured in fewer than that
    fraction of samples (0.0, the default, drops nothing beyond all-missing
    rows).  ``policy`` then decides the rest:

    ``as_is``  leave them for the backend, warning that something else will
               decide (gseapy fills with 0).  The default -- it changes no
               existing result.
    ``zero``   fill with 0 explicitly, making the current de-facto behaviour a
               stated choice rather than a side effect.
    ``low``    fill with one unit below each sample's own observed minimum, so
               an unmeasured feature ranks below every measured one.
    """
    if policy not in MISSING_POLICIES:
        raise ValueError(f"enrichment.missing must be one of {MISSING_POLICIES}, "
                         f"got {policy!r}")
    out = data
    n_feat_before = out.shape[0]

    frac_obs = out.notna().mean(axis=1)
    keep = frac_obs > 0 if min_observed_fraction <= 0 else frac_obs >= min_observed_fraction
    if not keep.all():
        out = out.loc[keep]
        if verbose:
            print(f"      missing data: dropped {n_feat_before - out.shape[0]} of "
                  f"{n_feat_before} features observed in <"
                  f"{max(min_observed_fraction, 0.0):.0%} of samples")

    n_missing = int(out.isna().sum().sum())
    if n_missing == 0:
        return out
    frac = n_missing / out.size
    if policy == "as_is":
        if verbose:
            warnings.warn(
                f"{n_missing} missing value(s) ({frac:.1%} of the matrix) reach "
                f"the enrichment backend with enrichment.missing='as_is'; "
                f"gseapy fills these with 0, which for a log-ratio matrix is "
                f"the middle of the ranking. Set enrichment.missing='low' to "
                f"rank unmeasured features below measured ones, or 'zero' to "
                f"state the current behaviour explicitly.",
                stacklevel=2)
        return out
    if policy == "zero":
        if verbose:
            print(f"      missing data: filled {n_missing} value(s) "
                  f"({frac:.1%}) with 0")
        return out.fillna(0.0)
    # policy == "low"
    floor = out.min(axis=0, skipna=True) - 1.0
    out = out.fillna(floor)
    if verbose:
        print(f"      missing data: filled {n_missing} value(s) ({frac:.1%}) "
              f"one unit below each sample's observed minimum")
    return out
