# Case study: paediatric solid tumours, normal tissues and cell lines

This directory reproduces every figure, table and quoted number of the
PathwayTheme Application Note from the public input data: 901 transcriptomes
(654 tumours over fifteen diagnoses, 118 normal tissues, 129 cell lines), ten
`pathwaytheme run` calls, one main figure, twelve supplementary figures and
twenty-one supplementary tables.

```
manuscript/
  notebooks/   the reproduction, run in order 00 -> 05
  analysis/    the code the notebooks call: one function per figure or table
  configs/     the ten pipeline configurations -- the record of what was run
  data/raw/        inputs (downloaded by notebook 00, checked against MANIFEST.sha256)
  data/genesets/   the 50 hallmark processes, GO:BP, and the GO-slim theme map
  data/processed/  the built matrix and sample annotation (written by notebook 01)
  results/     cached pathway scores and full pipeline output (written by notebook 02)
  nbl_mycn/    sample lists and runs of the MYCN analysis (written by notebook 02)
  submission/  the figures, tables and configurations of the paper
```

## Requirements

- Python ≥ 3.10 and the package with its notebook extra, from the repository root:
  `pip install -e ".[paper]"` (or `uv sync --extra paper`).
- R with the Bioconductor package `limma`, for Fig. S12 / Table S16 (the
  moderated *t* against `limma::eBayes`) and the `limma::camera` comparison of
  Supplementary Methods §12. Set `RSCRIPT` to the `Rscript` executable if it is
  not on the `PATH`.

## Inputs

| File (under `data/raw/`) | Source |
|---|---|
| `delivered/sarcoma_plus_NciLandscape.tpm.N_33967.260709T0925.tsv` | TPM, 33,967 Ensembl genes × 915 delivered samples |
| `delivered/sarcoma_plus_NciLandscape.pca_coords.log2_tpmPlusOne.N_33967.260705T1319.tsv` | the sequencing run sheet for the same samples |
| `gencode_v36_gene_map.tsv` | GENCODE v36 Ensembl → symbol and gene type |
| `brohl_supplement/NIHMS1759318-supplement-2.xlsx` | Brohl *et al.*, *Cell Rep* 2021;37:110047, Table S1 |
| `brohl_supplement/NIHMS1759318-supplement-3.xlsx` | Brohl *et al.* 2021, Table S2 |
| `gdc_clinical/clinical_TARGET-*.tsv` | open TARGET clinical tables, NCI Genomic Data Commons |

The inputs are deposited as one archive (DOI: *to be added*). Notebook
`00_inputs` downloads it into `data/raw/` and checks every file against
`data/raw/MANIFEST.sha256`; if the files are already in place it only checks
them.

## Run

Open the notebooks in Jupyter and run them in order, or execute all six
unattended from this directory:

```bash
for nb in notebooks/0*.ipynb; do
  jupyter nbconvert --to notebook --execute --inplace "$nb"
done
```

| Notebook | Produces |
|---|---|
| `00_inputs` | the verified inputs |
| `01_cohort` | the expression matrix (18,573 protein-coding genes × 901 samples, log2(TPM + 1)) and the sample annotation, joined to the external clinical tables |
| `02_pathway_runs` | the ten pipeline runs and the tables gathered from them |
| `03_unsupervised` | Figure 1, Figs. S1–S4; Tables S3, S4, S7–S9, S12, S15, S20 |
| `04_supervised` | Figs. S5–S8, S12; Tables S5, S6, S16–S19; the term-dependence and `camera` figures of Supplementary Methods §12 |
| `05_annotation` | Figs. S9–S11; Tables S1, S2, S10, S11, S13a, S13b, S14, S21 |

The GO:BP scoring pass in notebook 02 is the slow step (about 40 minutes on a
laptop); every later run reads the cached scores. Every figure is written as
PDF and PNG under `submission/figures/`, and every table as TSV under
`submission/tables/`.

To regenerate into a separate copy and compare it with the shipped outputs, set
`MANUSCRIPT_ROOT` to a directory holding `configs/`, `data/raw/` and
`data/genesets/` before starting Jupyter; every path the code reads or writes
then resolves under it.

## The cohort

Of the 915 delivered samples, the 14 from two non-public sequencing centres are
excluded, which also removes the only xenografts. The rest is:

| Sample type | n | Note |
|---|---|---|
| Tumour | 654 | 15 OncoTree diagnoses, plus 4 teratomas without a code |
| Cell line | 129 | ES 43, NBL 39, RMS 32, OS 9, ASPS 3, SYNS 3 |
| Normal | 118 | no run-sheet row; tissue of origin for 117 from Brohl Table S1A |

The technical variables are recorded rather than inferred, which is what makes
the cohort a test of separating biological from technical structure:

| Variable | Levels / range | Role |
|---|---|---|
| `sample_type` | TUMOR / CELL_LINE / NORMAL (a blank run-sheet Sample Type is normal tissue) | biological |
| `diagnosis` | 15 OncoTree codes | biological |
| `sex` | 85 recorded, 807 inferred from seven Y-linked genes, 9 unassigned | biological |
| `library_type` | polyA 575, ribo-zero 265, ribo-minus 36, polyA-stranded 25 | technical |
| `rin` | 4.8–10.0, 631 samples | technical |
| `seq_kit` | four kits, 596 samples | technical |
| `flowcell` | 92 runs, 780 samples | technical |
| `zero_fraction` | 0.025–0.307: the fraction of genes at zero TPM | technical |

Library chemistry, RIN and kit are completed from Brohl Table S1 and from the
flowcell where the run sheet is blank, and every completed value carries a
source column (`data/processed/provenance.tsv`, Table S2).

### Limits of the annotation

1. **The normals have no run-sheet row**, so every run-sheet field is blank on
   exactly those 118 samples. Brohl Table S1 closes most of it, but the tumours
   have no tissue of origin, so a tissue theme in the tumour-versus-normal
   contrast stays ambiguous between malignancy and tissue composition.
2. **Library chemistry is partly confounded with diagnosis.** Wilms tumour and
   CCSK are entirely ribo-depleted, hepatoblastoma entirely polyA; Ewing sarcoma
   and neuroblastoma span both chemistries, which is what lets chemistry resolve
   as a component of its own (Fig. S11).
3. **The twenty technical variables are not twenty independent probes**: 29 of
   their 190 pairs are associated at 0.70 or above (Fig. S10).
4. **Six fields are recorded on a small, non-random slice**: tissue of origin on
   the normals, disease at diagnosis on 85 osteosarcomas, and the four GDC
   demographic fields on 246 TARGET samples, all polyA. Each biological level
   is therefore compared with the samples in which its own field is recorded.
5. **Three rows read "Normal" in Diagnosis but TUMOR or CELL_LINE in Sample
   Type**; Sample Type is used, and the rows are named in Table S2.

`dev/` and `exploratory/` hold the authors' working tools and analyses that the
paper does not report; nothing in the reproduction depends on them.
