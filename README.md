# GLA variant analysis pipeline

Reproducible pipeline for the review "Генетическая гетерогенность GLA и границы применимости терапий при болезни Фабри" (рабочее название).
Principle: raw files -> scripts -> tables/figures. No manual edits of any file in data/interim, data/processed or results.

## Run
    conda env create -f environment.yml && conda activate gla
    make all          # whole pipeline
    make stage3       # single stage

## Stages (one script = one output)
| Stage | Script | Input (data/raw) | Output |
|---|---|---|---|
| 1 | 01_ingest_clinvar.py | ClinVar variant_summary (GLA) | interim/clinvar.csv |
| 1 | 02_ingest_gnomad.py | gnomAD GLA export | interim/gnomad.csv |
| 1 | 03_ingest_amenability.py | Amicus amenability table | interim/amenability.csv |
| 1 | 04_ingest_literature.py | curated literature CSV | interim/literature.csv |
| 2 | 05_harmonize.py | interim/*.csv, NM_000169 | processed/variants_master.csv + audit_report.md |
| 3 | 06_classify_cpg.py | variants_master | processed/variants_cpg.csv + level-1 tables |
| 4 | 07_structure.py | PDB, residue_map | processed/features_structural.csv |
| 4 | 08_ddg.py | structure + variants | processed/ddg.csv |
| 5 | 09_conservation.py | ortholog FASTA | processed/features_conservation.csv |
| 6 | 10_stats.py | all processed | results/tables/*.csv |
| 6b | 13_multivariable.py | processed features | results/tables/stats_G_* (GEE model, feature correlation, model comparison) |
| 6c | 14_active_site_sensitivity.py | processed features | results/tables/stats_H_active_site_sensitivity.csv |
| 7 | 11_figures.py | processed + results | results/figures/* |

Three analysis levels are kept separate: level 1 (CpG, stage 3), level 2 (structure, stages 4), level 3 (conservation, stage 5).
See docs/schema.md for the master table contract and docs/samples.md for sample definitions.

## Status
| Stage | State |
|---|---|
| 1-2 ingest + harmonise | done, tested. NM_000169.3 GenBank is used for reference-base checks and CDS-derived protein effects |
| 3 CpG (level 1) | done: `results/tables/cpg_enrichment.csv` |
| 4 structure (level 2) | done (3GXT): RSA, interface, ligand distance, glycosylation, active-site category with a cut-off sensitivity at 3.5/4.5/6.0 A. **ΔΔG is deliberately out of scope**: it is not computed and the manuscript says so; `scripts/08_ddg.py` only writes the submission list `results/ddg_input_dynamut2.csv` for a possible later extension |
| 5 conservation (level 3) | done: `scripts/09_conservation.py`, codon-aware MSA of 198 sequences -> 177 orthologs kept (QC in `results/tables/conservation_qc.csv`), Henikoff-weighted conservation; source: Ensembl release 116, downloaded 2026-09-11, already aligned at source |
| 6 statistics (levels 1-3) | done: `docs/stats_report.md`, plan in `docs/stats_plan.md` (tables A3, B3, F = level 3) |
| 6b multivariable model | done: logistic GEE clustered on residue, `results/tables/stats_G_*` |
| 6c active-site sensitivity | done: `results/tables/stats_H_active_site_sensitivity.csv` |
| 7 figures | done: `results/figures/fig1-5` (pdf/png/svg; fig5 = conservation). `fig6_structure.pml` is generated but not rendered here (no PyMOL): check it by eye |

`data/legacy/` holds the earlier merged tables; the pipeline never reads them except for the cross-check in the audit report.
Domain boundaries are intentionally not hard-coded: add them to config.yaml with a citation once confirmed.

## Figures
The drawing code lives in `src/gla/figures/`: `style.py` (palette, rcParams, saving, panel
letters, confidence-interval whiskers) and one module per figure (`fig1_map.py`,
`fig2_spectrum.py`, `fig3_structure.py`, `fig4_cpg.py`, `fig5_conservation.py`).
`scripts/11_figures.py` only selects which to build:

    python scripts/11_figures.py            # all figures
    python scripts/11_figures.py 3 4        # only figures 3 and 4
    python scripts/11_figures.py 1 --png    # one figure, PNG only (fast preview)

Change the look of every figure in `style.py`; change one figure in its own module.
Confidence intervals are drawn by `style.ci_x` / `style.ci_y`, which use capped whiskers.

## Obtaining the amenability table

One input is not redistributed: `data/raw/galafold_amenability.csv`, the Amicus Galafold
amenability table, which is proprietary. Everything else in `data/raw/` is public and is committed.
To run the pipeline from scratch:

1. Download the amenability table from the manufacturer and save it as
   `data/raw/galafold_amenability.csv`.
2. Verify it against the recorded checksum:

       python -c "import hashlib,sys;print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest())" data/raw/galafold_amenability.csv

   and compare with the `sha256` field of that file in `data/raw/MANIFEST.jsonl`. A different
   checksum means a different version of the table, and the coverage figures will differ.
3. `make all`.

Without that file, stages 1-2 fail; `data/processed/` and `results/` are committed, so every
downstream stage and every figure can still be rebuilt and inspected.

## Licensing

Source code: MIT (`LICENSE`). Derived data and figures: CC BY 4.0 (`LICENSE-DATA`).
Third-party inputs keep their own terms; see `LICENSE-DATA` for the per-source statement.

## Citing

See `CITATION.cff`. Cite the paper for the findings and the archived release (Zenodo DOI) for the
code and data.

## Manuscript
See `MANUSCRIPT.md` for the figure/table mapping between this repository and the manuscript,
and for the list of open items.
