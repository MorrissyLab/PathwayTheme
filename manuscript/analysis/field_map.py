"""Table S21: which analysis each annotated field enters, and why not where not.

For each of the 29 analysed fields: whether it is tested against the
components (the attribution grid of Supplementary Methods 6), against the other
variables (the association matrix of Methods 8.2), and level by level against
the technical variables (the confounding grid of Methods 8.3).  The three
recorded columns that are not analysed are listed too, so the table accounts
for every column a reader meets.

Membership is read from the shipped tables so it cannot drift from them.  Only
the reason a field is absent from the attribution grid is written here by hand,
because that is the part no table records.  Run after
:func:`analysis.annotation.run_all`.
"""
from __future__ import annotations

import pandas as pd

from .paths import PROC, STAB

BIOLOGICAL = ["sample_type", "diagnosis", "sex", "tissue_of_origin", "race",
              "ethnicity", "disease_at_diagnosis", "age_at_diagnosis_years",
              "vital_status"]

# Why a field is not in the attribution grid (Methods 6).  The grid runs inside
# the pipeline, on data/processed/metadata.tsv, so a field joined afterwards was
# never a candidate; where there is also a reason not to add it, that reason is
# given as well.
WHY = {
    "flowcell":
        "run sheet; excluded by configuration -- most of its levels fall "
        "below the three-sample floor",
    "sex":
        "joined after the pipeline; also derived rather than recorded, so "
        "adding it would change the grid's scope",
    "tissue_of_origin":
        "joined after the pipeline; recorded on the annotated normals only",
    "race": "joined after the pipeline; observed at a single chemistry",
    "ethnicity": "joined after the pipeline; observed at a single chemistry",
    "vital_status": "joined after the pipeline; observed at a single chemistry",
    "age_at_diagnosis_years":
        "joined after the pipeline; observed at a single chemistry",
    "disease_at_diagnosis":
        "joined after the pipeline; observed at a single chemistry",
}
# The alignment metrics' agreement with each other is the technical x technical
# block of the pairwise association matrix (Methods 8.2, Fig. S10).
QC_WHY = ("joined after the pipeline; one of fifteen alignment metrics, which "
          "Methods 8.2 shows are not fifteen independent variables")

# recorded columns that are not among the analysed fields
NOT_ANALYSED = (
    ("batch", "a re-spelling of sample type crossed with sequencing "
              "centre, which with one centre is sample type again"),
    ("sequencing_centre", "constant on the analysed cohort"),
    ("tumour/normal", "a re-spelling of sample type; not carried past "
                      "the build (Table S2)"),
)


def table_s21() -> pd.DataFrame:
    """Table S21: every annotated field against every analysis -- block,
    source, coverage, and whether it enters the attribution grid, the
    association matrix and the confounding grid, with the reason wherever it
    is left out of the attribution grid."""
    s11 = pd.read_csv(STAB / "tableS11_annotation_matrix.tsv", sep="\t")
    grid = set(pd.read_csv(STAB / "tableS3_metadata_pc_HALLMARK.tsv",
                           sep="\t")["variable"])
    s13a = pd.read_csv(STAB / "tableS13a_annotation_association.tsv", sep="\t")
    s13b = pd.read_csv(STAB / "tableS13b_confounding_levels.tsv", sep="\t")
    md = set(pd.read_csv(PROC / "metadata.tsv", sep="\t", nrows=1).columns)

    assoc = ((set(s13a["variable_a"]) | set(s13a["variable_b"]))
             - {"(excluded)"})
    lvl = set(s13b["biological"]) - {"(all)"}
    tech = set(s13b["technical"])

    fields = [c for c in s11.columns if c != "sample_id"]
    rows = []
    for f in fields:
        n = int(s11[f].notna().sum())
        in_grid = f in grid
        if in_grid:
            why = ""
        elif f.startswith("qc_"):
            why = QC_WHY
        else:
            why = WHY.get(f, "")
            if not why:
                raise ValueError(f"no reason recorded for {f!r} being outside "
                                 "the attribution grid")
        rows.append({
            "field": f,
            "block": "biological" if f in BIOLOGICAL else "technical",
            "source": ("recorded (Methods 9.1) and inferred from expression "
                       "(Methods 2)") if f == "sex"
                      else "run sheet or derived here" if f in md
                      else "joined (Methods 9)",
            "n_covered": n,
            "pct_of_901": round(100 * n / len(s11), 1),
            "attribution_grid": "yes" if in_grid else "no",
            "association_matrix": "yes" if f in assoc else "no",
            "confounding_grid": ("as a biological level" if f in lvl else
                                 "as a technical variable" if f in tech
                                 else "no"),
            "not_in_attribution_grid_because": why,
        })

    for f, why in NOT_ANALYSED:
        rows.append({"field": f, "block": "recorded, not analysed",
                     "source": "run sheet", "n_covered": pd.NA,
                     "pct_of_901": pd.NA, "attribution_grid": "no",
                     "association_matrix": "no", "confounding_grid": "no",
                     "not_in_attribution_grid_because": why})

    out = pd.DataFrame(rows)
    STAB.mkdir(parents=True, exist_ok=True)
    out.to_csv(STAB / "tableS21_field_analysis_map.tsv", sep="\t", index=False)

    analysed = out[out.block != "recorded, not analysed"]
    print(f"   tableS21: {len(analysed)} analysed fields + "
          f"{len(out) - len(analysed)} recorded but not analysed")
    return out


def run_all() -> None:
    """Write Table S21.  Needs Tables S3, S11, S13a and S13b."""
    table_s21()
