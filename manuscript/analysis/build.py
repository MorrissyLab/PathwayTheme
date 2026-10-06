"""Turn the delivered TPM matrix and run sheet into a PathwayTheme input pair.

Inputs, both in ``data/raw/delivered``:

    sarcoma_plus_NciLandscape.tpm.N_33967.260709T0925.tsv
        Ensembl id x sample, TPM.
    sarcoma_plus_NciLandscape.pca_coords.log2_tpmPlusOne.N_33967.260705T1319.tsv
        One row per sample: the data provider's own gene-level PC1-PC3,
        followed by the sequencing run sheet.  Only the run sheet is read; no
        analysis in this study uses a gene-level decomposition.

Two exclusion rules are applied here, and both are recorded in provenance.tsv:

  * Samples from any sequencing centre other than the public RNAseq Landscape
    cohort.  Thirteen delivered samples come from two other centres and are
    not part of the published resource, so they are dropped and the analysed
    cohort is the Landscape samples alone.  With one centre remaining,
    ``sequencing_centre`` is constant and ``batch`` -- a run-sheet column whose
    values are literally "<Sample Type>@<Sequencing Centre>" -- carries nothing
    that ``sample_type`` does not.  Neither is tested anywhere downstream.
  * Samples whose expression is more than 50% zeros.  The rule is applied, but
    it triggers on no analysed sample: the one delivered sample above it is
    already removed by the centre exclusion.

Two annotation fields are derived rather than read, because the delivered
column is unusable as delivered:

  * ``flowcell``.  The run sheet's "Run id" names a flowcell for the in-house
    samples and a concatenation of SRR run accessions -- unique per sample --
    for the TARGET osteosarcomas taken from SRA.  The accessions are discarded
    and the flowcell is read off the sample id instead, which also supplies it
    for the normals, who have no run-sheet row at all.
  * ``sex_inferred``.  The run sheet records no sex, and the external clinical
    tables record it only for osteosarcomas.  Y-linked expression assigns
    almost every sample, and :func:`analysis.link.link_supplement` checks the
    call against the recorded values before the field is used.

Nothing is imputed and nothing is batch-corrected: the analysis downstream
*measures* how much of the pathway-space variance the technical variables
carry, which requires leaving them in.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import paths

TPM = paths.DELIVERED / "sarcoma_plus_NciLandscape.tpm.N_33967.260709T0925.tsv"
COORDS = (paths.DELIVERED /
          "sarcoma_plus_NciLandscape.pca_coords.log2_tpmPlusOne.N_33967.260705T1319.tsv")

# GENCODE v36, taken from a GDC STAR-Counts file: it carries gene_id /
# gene_name / gene_type, so the Ensembl -> symbol map is local and versioned
# rather than fetched from a web service at build time.
GENCODE_MAP = paths.RAW / "gencode_v36_gene_map.tsv"

# A sample whose expression is mostly zeros cannot be scored meaningfully by any
# enrichment backend.  The worst analysed sample is near 31% zeros, so the
# threshold is not near anything real.
MAX_ZERO_FRACTION = 0.50

# Only the public RNAseq Landscape cohort is analysed.  Samples from other
# sequencing centres are dropped here rather than filtered downstream, so every
# artefact in the study is built from public samples alone.
PUBLIC_CENTRE = "RNAseq_Landscape_Manuscript"

# Analysis name <- run-sheet column.  Everything not listed is constant, an
# internal identifier, or empty for every sample.
KEEP = {
    "diagnosis": "OncoTree_Code",          # harmonised code; blank on the normals
    "diagnosis_label": "Diagnosis",        # spelled out, for figures
    "sample_type": "Sample Type",          # TUMOR / CELL_LINE / (blank = normal)
    "library_type": "Library Type",        # polya / ribozero / ribominus
    "rin": "RIN",
    "seq_kit": "Seq kit",
    "sequencing_centre": "Sequencing Centre",
    "batch": "Batch",
    "run_id": "Run id",
    "patient_id": "Patient ID",
}

# Cell-line panels large enough to contrast against their own tumours.
MATCHED_SUBTYPES = ("ES", "NBL", "RMS")

# Male-specific genes on the Y chromosome that survive the protein-coding
# filter.  XIST is the usual partner marker and is a lncRNA, so it is not in
# this matrix; the Y-linked panel alone separates the two groups completely on
# the samples where sex is recorded, so the missing partner costs nothing here.
Y_GENES = ("RPS4Y1", "DDX3Y", "UTY", "KDM5D", "USP9Y", "EIF1AY", "NLGN4Y")

# Thresholds on the mean log2(TPM+1) over Y_GENES.  On the samples with a
# recorded sex the two groups are separated by an empty interval -- female max
# 0.154, male min 1.009 -- so these are placed inside that gap rather than
# fitted to it, and the samples that land between them are left unassigned
# instead of being pushed to the nearer side.
Y_FEMALE_MAX = 0.20
Y_MALE_MIN = 1.00

# The run sheet's "Run id" is the flowcell for the samples sequenced in house,
# but for the TARGET osteosarcomas re-processed from SRA it is a concatenation
# of SRR run accessions, unique to each sample.  A per-sample unique value is
# not a batch: left in, it turns `flowcell` into one-sample strata that inflate
# every association it takes part in.  Those are set to NA.  The flowcell is
# instead recovered from the sample id, whose last field is the flowcell for
# every sample that has one -- including the normals, which have no run-sheet
# row at all.
FLOWCELL_RE = r"[A-Z0-9]{9,10}"


def gencode_map() -> pd.DataFrame:
    """Unversioned Ensembl gene id -> GENCODE v36 symbol and gene type."""
    if not GENCODE_MAP.exists():
        raise FileNotFoundError(f"no GENCODE reference at {GENCODE_MAP}")
    g = pd.read_csv(GENCODE_MAP, sep="\t",
                    usecols=["gene_id", "gene_name", "gene_type"])
    g = g[g["gene_id"].astype(str).str.startswith("ENSG")].copy()
    # strip the version suffix; the delivered matrix is unversioned
    g["ens"] = g["gene_id"].str.split(".").str[0]
    # a handful of PAR_Y genes share an unversioned id
    return g.drop_duplicates("ens").set_index("ens")


def _sparsity_rule_text(failed: pd.Series, n: int, worst: float,
                        delivered_over: pd.Series, delivered_type: pd.Series) -> str:
    """Provenance wording for the >50%-zeros rule."""
    if len(failed):
        return (f"EXCLUDED - {len(failed)} sample(s): " + ", ".join(failed.index))
    text = (f"NOT TRIGGERED - no sample of the {n} exceeds "
            f"{MAX_ZERO_FRACTION:.0%} zeros (worst {worst:.1%}).")
    if len(delivered_over) == 1:
        sid = delivered_over.index[0]
        text += (f"  The one sample in the delivered file that does, at "
                 f"{delivered_over.iloc[0]:.1%}, is a {delivered_type[sid]} from a "
                 "non-public centre and is already removed by the exclusion above")
    elif len(delivered_over):
        text += (f"  The {len(delivered_over)} samples in the delivered file that "
                 "do are from non-public centres and are already removed by the "
                 "exclusion above")
    return text


def build_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the expression matrix and sample annotation the pipeline runs on.

    Writes into ``data/processed``:

        expression_log2.tsv   HGNC symbol x sample, log2(TPM+1), protein-coding
        metadata.tsv          one row per sample, the analysis-ready annotation
        cohort_summary.tsv    diagnosis x sample type
        qc_samples.tsv        per-sample zero fraction and library type
        provenance.tsv        where every column came from, and what is missing

    ``provenance.tsv`` is Table S2 before its external rows are completed by
    :func:`analysis.link.link_supplement`.  Returns ``(expression, metadata)``.
    """
    paths.PROC.mkdir(parents=True, exist_ok=True)
    for p in (TPM, COORDS):
        if not p.exists():
            raise FileNotFoundError(f"missing input: {p}")

    print("1. gene annotation (GENCODE v36, local)")
    ann = gencode_map()
    print(f"   {len(ann)} genes, {int(ann.gene_type.eq('protein_coding').sum())} "
          "protein-coding")

    print("2. TPM matrix")
    X = pd.read_csv(TPM, sep="\t", index_col=0).astype("float32")
    print(f"   {X.shape[0]} genes x {X.shape[1]} samples")
    a = ann.reindex(X.index)
    unmapped = int(a["gene_name"].isna().sum())
    print(f"   unmapped ensembl ids: {unmapped}")

    pc = a["gene_type"].eq("protein_coding").to_numpy()
    X = X[pc]
    X.index = pd.Index(a.loc[pc, "gene_name"].to_numpy(), name="gene")
    # symbols are unique within protein-coding GENCODE v36, so no collapsing
    # step is needed; checked rather than assumed
    if X.index.has_duplicates:
        dup = sorted(set(X.index[X.index.duplicated()]))
        raise ValueError(f"duplicate protein-coding symbols: {dup[:10]}")
    print(f"   protein-coding, symbol-indexed: {X.shape[0]} genes")

    # zero fraction over every delivered sample, so the provenance can say
    # which delivered samples the sparsity rule would have caught
    zero_all = pd.Series((X.to_numpy() == 0).mean(axis=0), index=X.columns)

    print("3. run sheet, restricted to the public cohort")
    raw = pd.read_csv(COORDS, sep="\t", dtype=str).replace("NA", pd.NA)
    raw = raw.set_index("sample_id")
    missing = [c for c in KEEP.values() if c not in raw.columns]
    if missing:
        raise KeyError(f"run-sheet columns absent: {missing}")
    centre = raw["Sequencing Centre"].reindex(X.columns)
    other = sorted(centre[centre.ne(PUBLIC_CENTRE) | centre.isna()].index)
    by_centre = centre.loc[other].value_counts().to_dict() if other else {}
    if other:
        print(f"   dropping {len(other)} sample(s) from non-public centres: "
              f"{by_centre}")
        X = X.drop(columns=other)
    print(f"   {X.shape[1]} samples from {PUBLIC_CENTRE}")

    print("4. per-sample sparsity")
    zero_fraction = pd.Series((X.to_numpy() == 0).mean(axis=0), index=X.columns)
    print(f"   zero fraction: median {zero_fraction.median():.3f}, "
          f"max {zero_fraction.max():.3f}")
    failed = zero_fraction[zero_fraction > MAX_ZERO_FRACTION]
    if len(failed):
        print(f"   dropping {len(failed)} sample(s) above "
              f"{MAX_ZERO_FRACTION:.0%} zeros: {list(failed.index)}")
        X = X.drop(columns=list(failed.index))
    delivered_over = zero_all.loc[other][zero_all.loc[other] > MAX_ZERO_FRACTION]

    print("5. log2(TPM+1)")
    X = np.log2(X + 1.0)

    print("6. metadata")
    meta = pd.DataFrame(index=X.columns)
    meta.index.name = "sample_id"
    for new, old in KEEP.items():
        meta[new] = raw[old].reindex(meta.index)
    meta["rin"] = pd.to_numeric(meta["rin"], errors="coerce")
    meta["zero_fraction"] = zero_fraction.reindex(meta.index).round(5)

    # flowcell, reconstructed -- see FLOWCELL_RE above
    tail = pd.Series(meta.index, index=meta.index).str.rsplit("_", n=1).str[-1]
    fc = raw["Run id"].reindex(meta.index)
    n_srr = int(fc.astype(str).str.contains("SRR", na=False).sum())
    fc = fc.mask(fc.astype(str).str.contains("SRR", na=False))
    from_id = tail.where(tail.str.fullmatch(FLOWCELL_RE, na=False)
                         & ~tail.str.startswith("SRR"))
    n_from_id = int((fc.isna() & from_id.notna()).sum())
    meta["flowcell"] = fc.fillna(from_id)
    meta["flowcell_source"] = np.where(
        fc.notna(), "run sheet",
        np.where(meta["flowcell"].notna(), "sample id", "not recorded"))
    print(f"   flowcell: {int(meta['flowcell'].notna().sum())} samples, "
          f"{int(meta['flowcell'].nunique())} runs "
          f"({n_srr} SRA-accession values dropped, "
          f"{n_from_id} recovered from the sample id)")

    # sex, inferred from Y-chromosome expression -- see Y_GENES above.  The
    # cohort records a sex only for osteosarcomas, which is too little to test
    # anything against; the expression matrix carries the answer for almost
    # every sample, and the linking step checks the call against the recorded
    # values before the field is used anywhere.
    present = [g for g in Y_GENES if g in X.index]
    if len(present) < 4:
        raise ValueError(f"only {len(present)} Y-linked markers present; "
                         "sex cannot be inferred")
    yscore = X.loc[present].mean(axis=0)
    meta["y_chr_expression"] = yscore.reindex(meta.index).round(4)
    sex = pd.Series(pd.NA, index=meta.index, dtype=object)
    sex[yscore > Y_MALE_MIN] = "Male"
    sex[yscore < Y_FEMALE_MAX] = "Female"
    meta["sex_inferred"] = sex
    n_amb = int(sex.isna().sum())
    print(f"   sex inferred from {len(present)} Y-linked genes: "
          f"{int(sex.eq('Male').sum())} male, {int(sex.eq('Female').sum())} "
          f"female, {n_amb} between the thresholds and left unassigned")

    # A blank Sample Type means normal tissue, not "unknown": those rows carry
    # no run-sheet entry at all beyond the sequencing centre, their sample ids
    # are NS0xx / *normal / *muscle, and their sparsity is tissue-like rather
    # than cell-line-like.  Making it explicit is what gives the cohort a
    # tumour-versus-normal contrast.
    n_blank = int(meta["sample_type"].isna().sum())
    meta["sample_type"] = meta["sample_type"].fillna("NORMAL")
    print(f"   blank Sample Type -> NORMAL: {n_blank} sample(s)")

    # matched contrast groups, so tumour-versus-cell-line can be asked inside
    # one diagnosis instead of across the whole cohort
    grp = pd.Series("other", index=meta.index, dtype=object)
    for dx in MATCHED_SUBTYPES:
        for st, tag in (("TUMOR", "Tumor"), ("CELL_LINE", "CellLine")):
            m = meta["diagnosis"].eq(dx) & meta["sample_type"].eq(st)
            grp[m] = f"{dx}_{tag}"
    meta["contrast_group"] = grp

    X.to_csv(paths.PROC / "expression_log2.tsv", sep="\t")
    meta.to_csv(paths.PROC / "metadata.tsv", sep="\t")
    print(f"   {meta.shape[0]} samples x {meta.shape[1]} annotation columns")

    print("7. summaries")
    cohort = pd.crosstab(meta["diagnosis"].fillna("unlabelled"),
                         meta["sample_type"].fillna("unknown"), dropna=False)
    cohort["total"] = cohort.sum(axis=1)
    cohort = cohort.sort_values("total", ascending=False)
    cohort.to_csv(paths.PROC / "cohort_summary.tsv", sep="\t")
    print(cohort.to_string())

    qc = pd.DataFrame({
        "zero_fraction": meta["zero_fraction"],
        "library_type": meta["library_type"],
        "sample_type": meta["sample_type"],
        "rin": meta["rin"],
    })
    qc.to_csv(paths.PROC / "qc_samples.tsv", sep="\t")

    # a run-sheet inconsistency worth surfacing rather than silently resolving
    mislabelled = raw.index[raw["Diagnosis"].eq("Normal")
                            & raw["Sample Type"].isin(["TUMOR", "CELL_LINE"])]
    mislabelled = [s for s in mislabelled if s in meta.index]

    tumor_normal = raw["Tumor|Normal"].reindex(meta.index)

    prov = pd.DataFrame([
        {"item": "expression",
         "source": f"delivered {TPM.name}",
         "unit": "TPM -> log2(TPM+1), protein-coding only",
         "access": "delivered"},
        {"item": "gene identifiers",
         "source": "Ensembl -> HGNC symbol via GENCODE v36 (GDC STAR-Counts "
                   "header already in this repo)",
         "unit": f"{X.shape[0]} unique symbols, {unmapped} unmapped ids",
         "access": "local"},
        {"item": "diagnosis / sample type / library type / RIN / batch / "
                 "sequencing centre",
         "source": f"delivered {COORDS.name} (sequencing run sheet)",
         "unit": "categorical, RIN numeric", "access": "delivered"},
        {"item": "gene-level PC1-PC3",
         "source": "NOT USED - the delivered file carries the data provider's "
                   "own gene-level PC1-PC3; they are read past and no analysis "
                   "in this study uses a gene-level decomposition",
         "unit": "-", "access": "not used"},
        {"item": "sequencing_centre",
         "source": f"delivered {COORDS.name}; CONSTANT after the exclusion "
                   f"below - every analysed sample is {PUBLIC_CENTRE}",
         "unit": "single level, therefore not testable and not tested",
         "access": "delivered"},
        {"item": "batch",
         "source": "delivered run-sheet column whose values are literally "
                   "'<Sample Type>@<Sequencing Centre>'; blank on the normals, "
                   "whose Sample Type field is blank.  With one centre it is "
                   "sample_type re-spelled, so it carries no information of its "
                   "own and is excluded from every figure and every test",
         "unit": "recorded, never tested", "access": "delivered"},
        {"item": "samples from non-public sequencing centres",
         "source": "EXCLUDED - " + (
             "; ".join(f"{k} n={v}" for k, v in sorted(by_centre.items()))
             if other else "none")
                   + ". Not part of the published RNAseq Landscape resource",
         "unit": f"{len(other)} sample(s)", "access": "excluded"},
        {"item": "flowcell",
         "source": "delivered run sheet 'Run id' where it names a flowcell, "
                   f"otherwise the last field of the sample id.  The {n_srr} "
                   "TARGET osteosarcomas carry a concatenation of SRR run "
                   "accessions there instead, unique per sample; those are set "
                   f"to NA rather than treated as {n_srr} one-sample batches",
         "unit": f"{int(meta['flowcell'].nunique())} runs over "
                 f"{int(meta['flowcell'].notna().sum())} samples "
                 f"({n_from_id} recovered from the sample id, {n_srr} "
                 "SRA-accession values discarded)",
         "access": "delivered + derived"},
        {"item": "sex_inferred",
         "source": "derived: mean log2(TPM+1) over "
                   + ", ".join(present)
                   + f"; > {Y_MALE_MIN} male, < {Y_FEMALE_MAX} female, "
                     "between the two left unassigned.  The run sheet records "
                     "no sex, and a recorded sex exists for only 85 samples "
                     "(Brohl Table S1D), too few to test against; the call is "
                     "checked against those 85 before it is used",
         "unit": f"{int(meta['sex_inferred'].notna().sum())} of {len(meta)} "
                 f"assigned, {n_amb} unassigned",
         "access": "derived"},
        {"item": "zero_fraction",
         "source": "derived: fraction of protein-coding genes at exactly 0 TPM",
         "unit": "continuous, 0-1; the sparsity variable the association test "
                 "can be run against",
         "access": "derived"},
        {"item": f"sample(s) above {MAX_ZERO_FRACTION:.0%} zeros",
         "source": _sparsity_rule_text(failed, len(meta), zero_fraction.max(),
                                       delivered_over, raw["Sample Type"]),
         "unit": f"worst retained {zero_fraction.max():.3f}",
         "access": "excluded" if len(failed) else "not triggered"},
        {"item": "normal tissue",
         "source": f"derived: {n_blank} samples carry no Sample Type, which "
                   "denotes normal tissue (ids NS0xx / *normal / *muscle, "
                   "tissue-like sparsity, no run-sheet entry)",
         "unit": "sample_type == NORMAL", "access": "derived"},
        {"item": "tissue of origin for the normals",
         "source": "EXTERNAL - not in the delivered material; joined from "
                   "Brohl Supplementary Table S1A at the linking step, "
                   "which rewrites this row with the coverage it achieves.  "
                   "The tumour-versus-normal contrast as run pools the normals "
                   "and does not use it",
         "unit": "filled at the linking step", "access": "external"},
        {"item": "Tumor|Normal",
         "source": "RECORDED, NEVER TESTED - a run-sheet column whose two "
                   "levels are tumour and cell line, blank on the normals: it "
                   "is sample_type with one level removed, so against "
                   "sample_type it would contribute a guaranteed maximum "
                   "rather than a finding.  Not carried past the build",
         "unit": f"{int(tumor_normal.notna().sum())} of {len(meta)} non-blank, "
                 f"{int(tumor_normal.nunique())} levels",
         "access": "delivered"},
        {"item": "diagnosis for the normals",
         "source": "NOT APPLICABLE - left as NA, so the diagnosis association "
                   "test skips them rather than being given a label",
         "unit": "-", "access": "missing"},
        {"item": "Diagnosis == 'Normal' on non-normal rows",
         "source": "RUN-SHEET INCONSISTENCY - "
                   + (", ".join(mislabelled) if mislabelled else "none")
                   + "; Sample Type is used instead, and these stay tumour / "
                     "cell line",
         "unit": "-", "access": "flagged"},
    ])
    prov.to_csv(paths.PROC / "provenance.tsv", sep="\t", index=False)

    print("\n   library type x sparsity")
    print(meta.groupby("library_type")["zero_fraction"]
          .agg(["count", "median"]).round(3).to_string())
    print("\n   matched contrast groups")
    print(meta["contrast_group"].value_counts().to_string())
    print(f"\n  wrote {paths.PROC / 'expression_log2.tsv'} and metadata.tsv")
    return X, meta


def run_all() -> None:
    """Build the analysis inputs."""
    build_matrix()
