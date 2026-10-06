"""Link Brohl et al. 2021 supplementary Table S1 and TARGET clinical data.

Brohl AS, Sindiri S, Wei JS, *et al.* "Immuno-transcriptomic profiling of
extracranial pediatric solid malignancies." *Cell Rep* 2021;37(8):110047.
DOI 10.1016/j.celrep.2021.110047 · PMID 34818552 · PMCID PMC8642810.

That paper describes the cohort the two delivered files come from, so its
Table S1 carries per-sample annotation the run sheet does not.  Nothing here
re-analyses their results; only the sample annotation is taken.

``NIHMS1759318-supplement-2.xlsx`` (Table S1) has four sheets, each with a title
row above the real header:

    S1A  Patient Samples     patient / sample id, diagnosis, library type, and
                             **tissue sub type** -- which is where the normals'
                             tissue of origin lives (Normal-cerebellum, ...)
    S1B  RNA-Seq_QC          RIN plus alignment/library QC metrics
    S1C  Column description  documentation for S1B; not linked
    S1D  OS_Clinical         Gender / Race / Ethnicity / Disease at diagnosis /
                             age / survival -- **osteosarcoma only**, keyed by
                             TARGET USI

Because S1D is osteosarcoma-only, open TARGET clinical tables from the GDC
(``data/raw/gdc_clinical/``) are joined as a second source, keyed on the TARGET
USI that the neuroblastoma and osteosarcoma samples carry as their run-sheet
``patient_id``.  That extends race, ethnicity, vital status and age to the
TARGET neuroblastomas; the two sources are compared where they overlap.  Sex is
not taken from the GDC: its ``gender`` column is present and empty in all five
projects.  Sex is instead recorded from S1D where it exists and inferred from
Y-linked expression elsewhere, with the inferred call checked against S1D.

Values are canonicalised against declared vocabularies before merging, because
the two sources spell the same categories differently and would otherwise
produce two levels for one category.
"""
from __future__ import annotations

import re
import sys

import numpy as np
import pandas as pd

from . import paths

COORDS = (paths.DELIVERED /
          "sarcoma_plus_NciLandscape.pca_coords.log2_tpmPlusOne.N_33967.260705T1319.tsv")
TABLE_S1 = paths.BROHL / "NIHMS1759318-supplement-2.xlsx"

# sheet -> (id columns in preference order, {source column: our name})
SHEETS: dict[str, tuple[list[str], dict[str, str]]] = {
    "Table S1A_Patient Samples": (
        ["SampleID", "Patient ID "],
        {"Tissue Sub Type": "brohl_tissue_subtype",
         "Diagnosis": "brohl_diagnosis",
         "Diagnosis Abbreviation": "brohl_diagnosis_code",
         "RNAseq Library Type": "brohl_library_type"},
    ),
    "Table S1B.RNA-Seq_QC": (
        ["SampleID", "Sample.Data.ID", "Patient.ID"],
        {"Tissue_type": "brohl_tissue_type",
         "RNA_RIN": "brohl_rin",
         "Mapping Rate": "qc_mapping_rate",
         "Unique Rate of Mapped": "qc_unique_rate",
         "Duplication Rate of Mapped": "qc_duplication_rate",
         "Exonic Rate": "qc_exonic_rate",
         "Intronic Rate": "qc_intronic_rate",
         "Intragenic Rate": "qc_intragenic_rate",
         "Intergenic Rate": "qc_intergenic_rate",
         "Fragment Length Mean": "qc_fragment_length_mean",
         "Read Length": "qc_read_length",
         "Estimated Library Size": "qc_library_size",
         "Mapped": "qc_mapped",
         "Mapped Unique": "qc_mapped_unique",
         "Base Mismatch Rate": "qc_base_mismatch_rate",
         "Mean Per Base Cov.": "qc_mean_per_base_cov",
         "rRNA": "qc_rrna"},
    ),
    "Table S1D. OS_Clinical": (
        ["TARGET USI/Patient ID"],
        {"Gender": "sex",
         "Race": "race",
         "Ethnicity": "ethnicity",
         "Disease at diagnosis": "disease_at_diagnosis",
         "Age at Diagnosis in Days": "age_at_diagnosis_days",
         # S1D's outcome columns -- first event, overall and event-free
         # survival, primary tumour site, histologic response -- are not taken.
         # Nothing in the study tests an outcome, so joining them would put
         # columns into the annotation that no result and no provenance row
         # accounts for.
         "Vital Status": "vital_status"},
    ),
}

NUMERIC_OUT = {"brohl_rin", "age_at_diagnosis_days"}

# Open TARGET clinical from the GDC, one file per project, keyed by
# ``case_submitter_id`` = ``TARGET-<project>-<USI>``.  The neuroblastoma and
# osteosarcoma samples carry that USI as their run-sheet ``patient_id``, so
# these files reach samples Brohl's Table S1D does not: S1D is
# osteosarcoma-only, and without this join race, ethnicity, vital status and
# age are blank for every other diagnosis.
#
# `gender` is deliberately not taken from here: the column exists in all five
# files and is empty in all five.
GDC_FIELDS = {"race": "race", "ethnicity": "ethnicity",
              "vital_status": "vital_status",
              "age_at_diagnosis": "age_at_diagnosis_days",
              "primary_diagnosis": "gdc_primary_diagnosis"}

# The two sources spell the same categories differently -- Brohl writes
# "Black or African American" and "Not Hispanic or Latino", the GDC writes them
# lower-case, and S1D is internally inconsistent about
# "Non-metastatic (confirmed)" vs "(Confirmed)".  Left alone that produces two
# levels for one category, which would silently split every count and every
# contrast.  Canonical spellings are therefore declared, matched
# case-insensitively; anything unrecognised is passed through and reported so a
# new level cannot be absorbed unnoticed.
CANONICAL = {
    "race": ["White", "Black or African American", "Asian",
             "Native Hawaiian or Other Pacific Islander",
             "American Indian or Alaska Native", "Unknown", "Not Reported"],
    "ethnicity": ["Hispanic or Latino", "Not Hispanic or Latino", "Unknown",
                  "Not Reported"],
    "vital_status": ["Alive", "Dead", "Unknown", "Not Reported"],
    "disease_at_diagnosis": ["Metastatic", "Metastatic (confirmed)",
                             "Non-metastatic (confirmed)", "Unknown"],
    "sex": ["Male", "Female", "Unknown"],
}
_CANON = {f: {v.lower(): v for v in vals} for f, vals in CANONICAL.items()}

# Brohl S1A writes skeletal muscle two ways -- "Normal-Skeletal muscle" and
# "Normal-muscle" -- which would otherwise be two levels of the tissue field,
# splitting the largest normal tissue in half.
TISSUE_ALIASES = {"Muscle": "Skeletal muscle"}

# Run-sheet chemistry vocabulary <- Brohl S1A's.  Used only to fill samples the
# run sheet does not cover; where both record a value the run sheet is kept.
CHEMISTRY = {"PolyA": "polya", "PolyA Stranded": "polya_stranded",
             "Ribozero Stranded": "ribozero"}


def _canon(field: str, value, unrecognised: set) -> str:
    """Canonical spelling for a controlled-vocabulary field."""
    s = str(value).strip()
    table = _CANON.get(field)
    if table is None or not s:
        return s
    hit = table.get(s.lower())
    if hit is None:
        unrecognised.add((field, s))
        return s
    return hit


def _same_value(a, b) -> bool:
    """Do two raw strings mean the same thing? Numerics compared as numbers."""
    sa, sb = str(a).strip(), str(b).strip()
    try:
        return abs(float(sa) - float(sb)) < 1e-9
    except (TypeError, ValueError):
        return sa.lower() == sb.lower()


def _norm_id(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def candidate_keys(sample_id: str, patient_id) -> list[str]:
    """Progressively looser join keys for one of our samples, best first.

    Our ids are ``<name>_<T|N>_<flowcell>`` (``ASPS001tumor_T_D1PR6ACXX``);
    Brohl's are the bare name (``ASPS001``).  The run sheet's own ``patient_id``
    is tried first where present, since that is the same identifier space.
    """
    keys: list[str] = []
    if pd.notna(patient_id):
        keys.append(_norm_id(patient_id))
    sid = str(sample_id)
    keys.append(_norm_id(sid))
    m = re.match(r"^(.*?)_[TN]_[A-Za-z0-9]+$", sid)
    if m:
        keys.append(_norm_id(m.group(1)))
    keys.append(_norm_id(sid.split("_")[0]))
    # Trailing words the run sheet appends to the biosample name.  "muscle" is
    # here because four normals are named NS128muscle / NS129muscle /
    # NS132muscle / NS134muscle while Brohl S1A calls them NS128 ... NS134 --
    # without it those four are the only normals in the cohort with no external
    # annotation at all.
    for head in list(keys):
        for suffix in ("tumor", "tumour", "normal", "seq", "dtp", "cellline",
                       "muscle"):
            if head.endswith(suffix) and len(head) > len(suffix):
                keys.append(head[: -len(suffix)])
    # A trailing passage number on a cell line: FUUR1p28 is Brohl's FUUR1.
    # Tried last, so an exact id always wins over this.
    for head in list(keys):
        m2 = re.fullmatch(r"(.*[a-z0-9]{3,})p\d{1,3}", head)
        if m2:
            keys.append(m2.group(1))
    out, seen = [], set()
    for k in keys:
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


def _load_sheet(sheet: str) -> pd.DataFrame:
    """Read a Table S1 sheet, skipping its title row."""
    df = pd.read_excel(TABLE_S1, sheet_name=sheet, dtype=str, header=1)
    return df.dropna(axis=1, how="all")


def _fill_from_second_source(meta: pd.DataFrame,
                             annot: pd.DataFrame) -> list[dict]:
    """Complete the run sheet's technical fields from Brohl S1 and the flowcell.

    The run sheet has a row for the tumours and cell lines and no row at all
    for the normals, so every field it carries is blank on the normals *by
    construction* -- not at random, and not because a measurement was lost.
    Two of those fields are in the analysis (`library_type`, `rin`), and a
    third (`seq_kit`) is also blank on many tumours and cell lines.

    Three fills, each with a stated basis, each recorded in a `*_source` column
    so a reader can drop back to the delivered values:

      library_type  Brohl S1A "RNAseq Library Type", mapped onto the run
                    sheet's own vocabulary.  Blanks only: where both record a
                    value the run sheet wins, and the disagreements are counted
                    rather than reconciled.
      rin           Brohl S1B "RNA_RIN", blanks only.  The two sources give an
                    identical number on almost every sample both cover, so
                    they are the same measurement rather than two scales.
      seq_kit       propagated along the flowcell.  A flowcell is one physical
                    run, so every sample on it was sequenced with one kit; every
                    flowcell that carries a kit at all carries exactly one, with
                    no conflicts, and that check is what makes this safe.

    Nothing is imputed from a model: each filled value is either recorded
    elsewhere for that same sample, or fixed by a physical fact about the run it
    was on.
    """
    raw = pd.read_csv(COORDS, sep="\t", dtype=str).replace("NA", pd.NA)
    raw = raw.set_index("sample_id").reindex(meta.index)
    rows: list[dict] = []

    # -- library chemistry
    base = raw["Library Type"]
    other = annot["brohl_library_type"].map(CHEMISTRY).reindex(meta.index)
    both = base.notna() & other.notna()
    disagree = int((base[both] != other[both]).sum())
    meta["library_type"] = base.fillna(other)
    meta["library_type_source"] = np.where(
        base.notna(), "run sheet",
        np.where(meta["library_type"].notna(), "Brohl S1A", "not recorded"))
    gained = int((base.isna() & other.notna()).sum())
    print(f"   library_type {int(base.notna().sum())} -> "
          f"{int(meta['library_type'].notna().sum())} (+{gained} from Brohl "
          f"S1A); {disagree}/{int(both.sum())} disagree where both record one, "
          "run sheet kept")
    rows.append({"kind": "fill", "field": "library_type",
                 "from_primary": int(base.notna().sum()), "filled": gained,
                 "total": int(meta["library_type"].notna().sum()),
                 "fill_source": "Brohl S1A RNAseq Library Type",
                 "disagreements_where_both": disagree,
                 "n_where_both": int(both.sum())})

    # -- RIN
    base = pd.to_numeric(raw["RIN"], errors="coerce")
    other = pd.to_numeric(annot["brohl_rin"], errors="coerce").reindex(meta.index)
    both = base.notna() & other.notna()
    same = int((base[both] == other[both]).sum())
    meta["rin"] = base.fillna(other)
    meta["rin_source"] = np.where(
        base.notna(), "run sheet",
        np.where(meta["rin"].notna(), "Brohl S1B", "not recorded"))
    gained = int((base.isna() & other.notna()).sum())
    print(f"   rin {int(base.notna().sum())} -> "
          f"{int(meta['rin'].notna().sum())} (+{gained} from Brohl S1B); "
          f"identical on {same}/{int(both.sum())} of the overlap")
    rows.append({"kind": "fill", "field": "rin",
                 "from_primary": int(base.notna().sum()), "filled": gained,
                 "total": int(meta["rin"].notna().sum()),
                 "fill_source": "Brohl S1B RNA_RIN",
                 "disagreements_where_both": int(both.sum()) - same,
                 "n_where_both": int(both.sum())})

    # -- sequencing kit, along the flowcell
    base = raw["Seq kit"]
    fc = meta["flowcell"]
    tab = pd.crosstab(fc, base)
    one_kit = {f: r.idxmax() for f, r in tab.iterrows() if int((r > 0).sum()) == 1}
    clash = [f for f, r in tab.iterrows() if int((r > 0).sum()) > 1]
    if clash:
        print(f"   WARNING {len(clash)} flowcell(s) carry more than one kit; "
              f"not propagated: {clash}", file=sys.stderr)
    meta["seq_kit"] = base.fillna(fc.map(one_kit))
    meta["seq_kit_source"] = np.where(
        base.notna(), "run sheet",
        np.where(meta["seq_kit"].notna(), "flowcell", "not recorded"))
    gained = int((base.isna() & meta["seq_kit"].notna()).sum())
    print(f"   seq_kit {int(base.notna().sum())} -> "
          f"{int(meta['seq_kit'].notna().sum())} (+{gained} propagated along "
          f"{len(one_kit)} single-kit flowcells, {len(clash)} conflicting)")
    rows.append({"kind": "fill", "field": "seq_kit",
                 "from_primary": int(base.notna().sum()), "filled": gained,
                 "total": int(meta["seq_kit"].notna().sum()),
                 "fill_source": f"flowcell ({len(one_kit)} single-kit runs)",
                 "disagreements_where_both": len(clash),
                 "n_where_both": int(fc.notna().sum())})
    return rows


def _check_inferred_sex(meta: pd.DataFrame,
                        annot: pd.DataFrame) -> tuple[pd.Series, pd.Series, dict]:
    """Score the Y-chromosome sex call against the sexes that are recorded.

    Only osteosarcomas carry a recorded sex (Brohl S1D), and they are the only
    test available, so they are used as one: the call has to reproduce every
    recorded value before the field it produces is allowed into a figure.  Sex
    is then reported as a single column, recorded where a source records it and
    inferred elsewhere, with the origin of each value kept beside it.
    """
    inferred = meta["sex_inferred"].reindex(annot.index)
    recorded = (annot["sex"] if "sex" in annot.columns
                else pd.Series(pd.NA, index=annot.index, dtype=object))
    both = recorded.notna() & inferred.notna()
    agree = int((recorded[both] == inferred[both]).sum())
    stat = {"kind": "check", "field": "sex",
            "from_primary": int(recorded.notna().sum()),
            "filled": int((recorded.isna() & inferred.notna()).sum()),
            "total": int(recorded.fillna(inferred).notna().sum()),
            "fill_source": "Y-linked expression, checked against Brohl S1D",
            "disagreements_where_both": int(both.sum()) - agree,
            "n_where_both": int(both.sum())}
    print(f"   sex: recorded for {stat['from_primary']}, the inferred call "
          f"agrees on {agree}/{int(both.sum())} of them; +{stat['filled']} "
          f"inferred, {stat['total']} of {len(annot)} now assigned")
    if both.sum() and agree < both.sum():
        print("   WARNING inferred sex disagrees with the recorded value on "
              f"{int(both.sum()) - agree} sample(s)", file=sys.stderr)
    merged = recorded.fillna(inferred)
    source = pd.Series(np.where(recorded.notna(), "Brohl S1D",
                                np.where(merged.notna(), "Y-linked expression",
                                         "unassigned")), index=annot.index)
    return merged, source, stat


def _refresh_provenance(annot: pd.DataFrame, meta: pd.DataFrame) -> None:
    """Rewrite the Table S2 rows that only the external join can fill.

    :func:`analysis.build.build_matrix` writes provenance.tsv first, so its
    rows for externally joined columns are placeholders until here.
    """
    path = paths.PROC / "provenance.tsv"
    if not path.exists():
        return
    prov = pd.read_csv(path, sep="\t")

    toi = (annot["tissue_of_origin"].dropna()
           if "tissue_of_origin" in annot.columns else pd.Series(dtype=object))
    n_norm = int((meta["sample_type"] == "NORMAL").sum())
    row = prov["item"] == "tissue of origin for the normals"
    if row.any():
        prov.loc[row, "source"] = (
            "Brohl Supplementary Table S1A, joined here - "
            f"{len(toi)} of the {n_norm} normals carry a tissue of origin.  "
            "The tumour-versus-normal contrast as run pools the normals and "
            "does not use it")
        prov.loc[row, "unit"] = f"{len(toi)} samples over {toi.nunique()} tissues"
        prov.loc[row, "access"] = "external"

    n_run_sheet = int(meta["sample_type"].isin(["TUMOR", "CELL_LINE"]).sum())
    filled = [c for c in ("library_type", "rin", "seq_kit") if c in meta.columns]
    prov = pd.concat([prov, pd.DataFrame([
        {"item": "run-sheet fields completed from a second source",
         "source": f"the run sheet stops at the {n_run_sheet} tumours and cell "
                   "lines; each field below is completed from Brohl Table S1 or "
                   "from the flowcell, where the two agree on the samples they "
                   "share (Methods 9.2).  The pipeline runs on the completed "
                   "columns, not the delivered ones",
         "unit": "; ".join(f"{c} -> {int(meta[c].notna().sum())} of "
                           f"{len(meta)}" for c in filled),
         "access": "delivered + external"}])], ignore_index=True)
    prov.to_csv(path, sep="\t", index=False)
    print(f"   refreshed the external rows of {path}")


def link_supplement() -> pd.DataFrame:
    """Join Brohl Table S1 and the GDC TARGET clinical tables to the cohort.

    Needs ``metadata.tsv`` and ``provenance.tsv`` from
    :func:`analysis.build.build_matrix`.  Rewrites ``metadata.tsv`` with
    library type, RIN and sequencing kit completed from a second source,
    checks the inferred sex against the recorded one, and writes into
    ``data/processed``:

        clinical_annotation.tsv     sample_id -> every field that matched
        supplement_link_report.tsv  which sheet and id column produced each
                                    field, the fills, and the coverage of each
                                    field over the cohort

    It also completes the external rows of ``provenance.tsv`` (Table S2).
    Returns the clinical annotation.
    """
    if not TABLE_S1.exists():
        raise FileNotFoundError(
            f"Brohl Table S1 missing at {TABLE_S1}; run "
            "analysis.fetch.fetch_inputs() first")

    meta = pd.read_csv(paths.PROC / "metadata.tsv", sep="\t",
                       dtype={"sample_id": str, "patient_id": str})
    ours = meta["sample_id"].tolist()
    keys_for = {r.sample_id: candidate_keys(r.sample_id, r.patient_id)
                for r in meta.itertuples()}
    unrecognised: set[tuple[str, str]] = set()

    wanted = [n for _, (_, fmap) in SHEETS.items() for n in fmap.values()]
    annot = pd.DataFrame(index=pd.Index(ours, name="sample_id"),
                         columns=wanted, dtype=object)
    report: list[dict] = []

    # ---- Brohl Table S1, sheet by sheet ----------------------------------
    for sheet, (idcols, fmap) in SHEETS.items():
        df = _load_sheet(sheet)
        print(f"-> {sheet}: {df.shape[0]} rows x {df.shape[1]} cols")
        missing = [c for c in list(idcols) + list(fmap) if c not in df.columns]
        if missing:
            print(f"   WARNING columns absent, skipped: {missing}",
                  file=sys.stderr)
        idcols = [c for c in idcols if c in df.columns]
        fields = {src: tgt for src, tgt in fmap.items() if src in df.columns}
        if not idcols or not fields:
            continue
        for idc in idcols:
            lut: dict[str, pd.Series] = {}
            for _, row in df.iterrows():
                k = _norm_id(row[idc])
                if k and k not in lut:
                    lut[k] = row
            hits = 0
            for sid in ours:
                for key in keys_for[sid]:
                    row = lut.get(key)
                    if row is None:
                        continue
                    took = False
                    for src, tgt in fields.items():
                        v = row[src]
                        if (pd.notna(v) and str(v).strip()
                                and pd.isna(annot.at[sid, tgt])):
                            annot.at[sid, tgt] = _canon(tgt, v, unrecognised)
                            took = True
                    if took:
                        hits += 1
                    break
            print(f"   via {idc:24s} matched {hits:4d}/{len(ours)}")
            if hits:
                report.append({"kind": "match", "sheet": sheet,
                               "id_column": idc, "samples_matched": hits,
                               "fields": ",".join(sorted(fields.values()))})

    # ---- open GDC TARGET clinical, joined on the USI ----------------------
    gdc: dict[str, pd.Series] = {}
    projects: dict[str, str] = {}
    for f in sorted(paths.GDC.glob("clinical_TARGET-*.tsv")):
        d = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
        proj = f.stem.replace("clinical_", "")
        for _, row in d.iterrows():
            mm = re.fullmatch(r"TARGET-\d+-([A-Z0-9]{6})",
                              str(row.get("case_submitter_id", "")))
            if mm and mm.group(1) not in gdc:
                gdc[mm.group(1)] = row
                projects[mm.group(1)] = proj
    if gdc:
        for tgt in GDC_FIELDS.values():
            if tgt not in annot.columns:
                annot[tgt] = pd.NA
        annot["gdc_project"] = pd.NA
        filled = agree = conflict = hits = 0
        conflicts: list[tuple] = []
        for sid, pid in zip(meta["sample_id"], meta["patient_id"]):
            row = gdc.get(str(pid))
            if row is None:
                continue
            hits += 1
            annot.at[sid, "gdc_project"] = projects[str(pid)]
            for src, tgt in GDC_FIELDS.items():
                v = row.get(src)
                if pd.isna(v) or not str(v).strip():
                    continue
                v = _canon(tgt, v, unrecognised)
                have = annot.at[sid, tgt]
                if pd.isna(have):
                    annot.at[sid, tgt] = v
                    filled += 1
                elif _same_value(have, v):
                    agree += 1
                else:
                    conflict += 1
                    conflicts.append((sid, tgt, have, v))
        print(f"-> GDC TARGET clinical ({len(gdc)} USIs across "
              f"{len(set(projects.values()))} projects)")
        print(f"   matched {hits}/{len(ours)} on run-sheet patient_id = USI; "
              f"filled {filled} blanks, {agree} agreed with Brohl, "
              f"{conflict} conflicted")
        if conflict:
            print(f"   WARNING {conflict} value(s) genuinely disagree between "
                  f"Brohl Table S1D and GDC; Brohl kept", file=sys.stderr)
            for sid, tgt, have, v in conflicts[:10]:
                print(f"     {sid} {tgt}: Brohl={have!r} GDC={v!r}",
                      file=sys.stderr)
        report.append({"kind": "match", "sheet": "GDC TARGET clinical",
                       "id_column": "case_submitter_id (USI)",
                       "samples_matched": hits,
                       "fields": ",".join(sorted(GDC_FIELDS.values()))})

    if unrecognised:
        print(f"   NOTE {len(unrecognised)} value(s) outside the declared "
              f"vocabularies, passed through as-is:")
        for f_, v_ in sorted(unrecognised):
            print(f"     {f_}: {v_!r}")

    annot = annot.dropna(axis=1, how="all")
    for c in annot.columns:
        if c in NUMERIC_OUT or c.startswith("qc_"):
            annot[c] = pd.to_numeric(annot[c], errors="coerce")
    if annot.notna().sum().sum() == 0:
        raise RuntimeError("no external annotation matched any sample")

    # a derived age in years, since days is what the source records
    if "age_at_diagnosis_days" in annot.columns:
        annot["age_at_diagnosis_years"] = (
            annot["age_at_diagnosis_days"] / 365.25).round(2)

    # Tissue of origin, for the normals only.  S1A writes it as
    # "Normal-cerebellum"; a tumour's or cell line's tissue is its diagnosis, so
    # those stay blank rather than being filled with "Tumor", which would add
    # levels to the variable and read as if they were recorded.
    if "brohl_tissue_subtype" in annot.columns:
        sub = annot["brohl_tissue_subtype"].astype("object")
        is_norm = sub.notna() & sub.str.startswith("Normal", na=False)
        tissue = pd.Series(
            sub.where(is_norm).str.replace(r"^Normal-\s*", "", regex=True),
            index=annot.index).str.capitalize()
        annot["tissue_of_origin"] = tissue.replace(TISSUE_ALIASES)

    # ---- one sex column: recorded where recorded, inferred elsewhere -----
    meta = meta.set_index("sample_id")
    annot["sex_recorded"] = annot.get("sex")
    annot["sex"], annot["sex_source"], sex_stat = _check_inferred_sex(meta, annot)

    # ---- complete the run sheet's technical fields ------------------------
    print("-> completing the run sheet")
    fills = _fill_from_second_source(meta, annot) + [sex_stat]
    meta.to_csv(paths.PROC / "metadata.tsv", sep="\t")
    print(f"   rewrote {paths.PROC / 'metadata.tsv'} with the completed fields")

    annot.to_csv(paths.PROC / "clinical_annotation.tsv", sep="\t")

    cov = pd.DataFrame({
        "kind": "coverage", "field": annot.columns,
        "n_annotated": [int(annot[c].notna().sum()) for c in annot.columns],
        "pct_of_cohort": [round(float(annot[c].notna().mean() * 100), 1)
                          for c in annot.columns],
        "n_levels": [int(annot[c].nunique(dropna=True)) for c in annot.columns],
    })
    print(f"\ncoverage over our {len(annot)} samples:")
    print(cov[["field", "n_annotated", "pct_of_cohort", "n_levels"]]
          .to_string(index=False))
    pd.concat([pd.DataFrame(report), pd.DataFrame(fills), cov],
              ignore_index=True).to_csv(
        paths.PROC / "supplement_link_report.tsv", sep="\t", index=False)

    _refresh_provenance(annot, meta)

    if "tissue_of_origin" in annot.columns:
        norm = annot["brohl_tissue_subtype"].dropna()
        norm = norm[norm.str.startswith("Normal")]
        print(f"\nnormals with a tissue of origin: {len(norm)} across "
              f"{norm.nunique()} tissues")
        print("  " + ", ".join(f"{k} {v}" for k, v in
                               norm.value_counts().items()))

    print("\nwrote", paths.PROC / "clinical_annotation.tsv")
    return annot


def run_all() -> None:
    """Join the external annotation (needs :func:`analysis.build.run_all`)."""
    link_supplement()
