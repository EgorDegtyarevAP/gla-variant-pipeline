# LOG (versions, download dates, decisions)
- 2026-10-21: pipeline rebuilt from raw data; legacy GLA_final_merged_v2.csv is NOT used as input (only as a cross-check of the new master table).
- 2026-10-26: stages 1-2 implemented. Inputs: ClinVar export (1990 rows; ClinVar release/date NOT recorded), gnomAD v4.1.1 export (2026-09-11), Amicus Galafold table (version/date NOT recorded; assay + thresholds not in file), 17 ClinVar genomic-only names.
- Decision: gnomAD canonical transcript ENST00000218516.4 assumed = NM_000169.3 (MANE Select) - verify.
- Decision: Amicus rows with "or" alternatives expanded to one row per nucleotide change; the 69 rows on the user's exclusion list are not merged (kept in results/tables/excluded_or_unresolved.csv).
- Decision: Amicus-absent variants get amenability NA (not non_amenable).
- Open: gnomAD lists 9 ClinVar IDs (incl. c.937G>T, c.335G>A) that are missing from the ClinVar export -> re-export from full ClinVar.
- Open: NM_000169 sequence + PDB 3GXT/3GXP needed for stages 3-4 (sandbox has no access to NCBI/Ensembl/RCSB).
- 2026-10-26: NM_000169.3 (GenBank, 10-SEP-2026, MANE Select, matches ENST00000218516.4) and PDB 3GXT (2.7 A, hDGJ complex, chains A/B) supplied. Verified: CDS translation == GenBank translation; exon ends == IVS table; 2440/2440 SNV reference bases match; protein effects computed from the CDS agree with all three sources for SNVs (2 indel notation conflicts remain).
- Decision: CpG sample for level 1 = gnomAD PASS and ClinVar P/LP CDS SNVs only (Amicus-only variants were chosen for in-vitro testing, not observed in nature).
- Decision: RSA from SASA of protein atoms only (Shrake-Rupley, probe 1.4 A, 200 points) on the A+B dimer; maxASA = Tien et al. 2013 (verify). Active site = within 4.5 A of bound ligand (D170 and D231 come out at 3.0 / 2.8 A).
- Open: DDG (DynaMut2, 1724 unique protein changes in results/ddg_input_dynamut2.csv; check batch size limit and sign convention); ortholog FASTA for stage 5; statistics must handle many substitutions per position (pseudo-replication).
- 2026-10-26: stage 6 (levels 1-2) implemented. Units: S_amen collapsed to unique protein change (1808; no status conflicts) with position-clustered permutation; positions for B; windows (20 aa) for E. Seed/permutations in config.yaml. BH within each table.
- 2026-10-26: stage 7 figures 1-4 (matplotlib; palette slots 1-3 validated with --pairs all, contrast WARN for aqua handled by direct labels + legends). fig5 PyMOL script generated, untested.
- 2026-10-26: stage 5 (level 3). Input Human_GLA_orthologues.fa (source/release NOT recorded; Ensembl-style IDs). Filters: coverage >= 0.80, identity to human >= 0.55 (rationale from memory, verify), one sequence per species; 198 -> 177 kept. Positions valid if weighted non-gap >= 0.70 (410/429). Conservation = 1 - H/log2(20) with Henikoff weights (deviation from the earlier 'unweighted' plan: dense mammal sampling; unweighted version is saved and agrees, Spearman 0.99).
- 2026-10-26: level-3 statistics (tables A3, B3, F) and fig5_conservation; conservation track added to fig1; PyMOL script renamed fig6_structure.pml. Stage 5 added to `make all`.
- Open: sensitivity at identity threshold 0.4 not yet run; domain boundaries; DDG.

## 2026-11-09 — reviewer round: reporting fix and two added analyses

- `src/gla/cpg.py`: added `n_cpg_ts_expected = n_obs * p_exp`, derived from the same class-specific
  null proportion used by the test, plus an assertion that the printed expectation and the tested
  null cannot diverge. The manuscript had carried expected counts of 2.2 / 1.0 / 5.8 / 1.1 for the
  four class-specific rows, computed by hand from the overall 0.98%; the correct values under the
  class-specific proportions (missense 0.8772%, synonymous 1.1016%, nonsense 2.1622%) are
  2.0 / 1.2 / 5.2 / 2.4. Folds and p values in the pipeline were already correct and are unchanged.
- `src/gla/cpg.py`: added `power_fold()` and `mde_fold()` (exact binomial), and the columns
  `mde_fold_80pct_power` and `power_at_observed_fold` to `cpg_enrichment.csv`. The clinical
  comparisons reach 80% power only from x2.28 (all CDS SNVs), x2.40 (missense) and x3.23 (nonsense);
  power against the population-sized x7.6 exceeds 0.99.
- `scripts/13_multivariable.py` (new, stage 6b): logistic GEE for amenability clustered on residue,
  exchangeable working correlation, 1697 changes at 385 positions. Distance to ligand loses its
  association once burial and active-site membership are in the model (OR 0.95, p = 0.66);
  conservation, BLOSUM62 and the alternative-residue term add over the structural model
  (Wald chi2 119.6, df 3; QIC 1641 -> 1434). Outputs `stats_G_*`.
- `scripts/14_active_site_sensitivity.py` (new, stage 6c): active-site category recomputed at
  3.5 / 4.5 / 6.0 A (10 / 15 / 20 residues). Non-amenable share among active-site P/LP missense with
  a status is 100% (31/31), 97.6% (41/42) and 98.1% (52/53); direction and significance unchanged.
  Outputs `stats_H_active_site_sensitivity.csv`.
- `dist_ligand` verified independently against `data/raw/3GXT.cif`: agreement is exact (0.000 A) for
  every residue within 8 A of the ligand; D170 3.00 A, D231 2.81 A.
- `environment.yml`: pinned `statsmodels=0.15.0`. `Makefile`: new target `stage6b`.
- Manuscript: thermodynamic stability (ddG) dropped from the claimed scope; Section 2.5 renamed and
  states explicitly that no stability calculation was performed. `scripts/08_ddg.py` is retained,
  but it only prepares the submission list (Table S5) and computes nothing.
- Full pipeline rerun end to end on this date; all manuscript figures reproduce unchanged.

