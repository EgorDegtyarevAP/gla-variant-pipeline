# Statistical plan (levels 1-2). Written before the final run; the first structure/amenability look was exploratory.
Seed and permutation counts: config.yaml (`stats`). Multiple testing: Benjamini-Hochberg within each table (family = all primary tests of that table). Effect sizes and 95% CIs are reported with every p-value.

| id | question | unit | test | statistic |
|---|---|---|---|---|
| A | do non-amenable and amenable protein changes (S_amen, collapsed to unique protein change) differ in structural context? | protein change, clustered by position | permutation of the position->feature assignment (10 000), cluster bootstrap CI (2 000) | AUC (continuous), difference in proportions (binary) |
| B | where do P/LP missense (S_path) and gnomAD PASS missense sit relative to all resolved positions? | position | Mann-Whitney / Fisher exact | AUC, odds ratio |
| C | structural category of S_path variants; share non-amenable per category | variant (descriptive) | Wilson 95% CI | proportion |
| D | ceiling of migalastat coverage among ClinVar P/LP | variant (descriptive) | Wilson 95% CI; lower/upper bound (NA = not in table) | proportion |
| E | regional clustering of P/LP missense relative to mutational opportunity | 20-aa window | max-statistic permutation (10 000), null = same number of variants drawn from all possible missense SNVs | observed vs expected, family-wise p |

Definitions: buried = RSA < 0.2 (convention; 0.1 and 0.3 as sensitivity); active site = heavy atom within 4.5 A of the bound ligand in 3GXT; interface = dSASA(monomer->dimer) >= 1 A^2; RSA/distance are chain means of the A/B dimer.
Known limits: Amicus variants were selected for in-vitro testing, ClinVar for clinical reporting, gnomAD for population presence: each comparison inherits its ascertainment. Position-level features cannot separate substitutions at one position.

## Level 3 (conservation) - added with stage 5
- Features: position-level conservation (1 - H/log2 20, Henikoff-weighted) and frac_human_aa (share of orthologs carrying the human residue, human row excluded); substitution-level BLOSUM62, ortholog frequency of the alternative residue, alt_seen.
- A3: amenable vs non-amenable (same unique protein changes as table A). Position-level features: position-clustered permutation (10 000) + cluster bootstrap CI. Substitution-level features: labels permuted within position (2 000). BH within table.
- B3: positions with P/LP (or gnomAD) missense vs other valid positions (Mann-Whitney at position level).
- F: the two senses of 'non-conserved': (a) spectrum heterogeneity (share of positions with P/LP) and (b) low positional conservation (distribution of conservation, share of positions with <50% orthologs matching human), by conservation quartile.
- Deviation: Henikoff weighting is used although the original plan said 'no weighting'; the unweighted version is reported as sensitivity (features_conservation_unweighted.csv).
- Not done: sensitivity at identity threshold 0.4.
