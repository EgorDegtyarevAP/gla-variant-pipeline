"""Stage 7. Multivariable model of amenability, clustered on residue.
IN: processed/variants_master, features_structural, features_conservation, features_substitution.
OUT: results/tables/stats_G_gee_amenability.csv, stats_G_feature_correlation.csv, stats_G_model_comparison.csv.

Why this stage exists: the level-2 and level-3 features are correlated (RSA vs conservation
Spearman rho = -0.50), so the univariate AUCs in stats_A cannot establish an independent
contribution of any one of them. Substitutions at the same residue share every position-level
feature, so the model clusters on the residue."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, pandas as pd, statsmodels.api as sm
from gla.common import PROCESSED, RESULTS

T = RESULTS / "tables"
CONT = ["rsa", "dist_ligand", "conservation", "blosum62"]

def build():
    M = pd.read_parquet(PROCESSED / "variants_master.parquet")
    S = pd.read_parquet(PROCESSED / "features_structural.parquet")
    C = pd.read_csv(PROCESSED / "features_conservation.csv")
    B = pd.read_parquet(PROCESSED / "features_substitution.parquet")
    R = S[S.structure_resolved].set_index("protein_pos")
    a = M[M.S_amen].drop_duplicates("protein")[["protein", "amenability_status"]].copy()
    a["pos"] = a.protein.str.extract(r"(\d+)")[0].astype(int)
    a = a[a.pos.isin(R.index)].reset_index(drop=True)
    a["y"] = (a.amenability_status == "non_amenable").astype(int)
    a = a.merge(R[["rsa", "dist_ligand", "active_site"]], left_on="pos", right_index=True, how="left")
    a = a.merge(C[["protein_pos", "conservation", "valid"]], left_on="pos", right_on="protein_pos", how="left")
    a = a.merge(B[["protein", "blosum62", "alt_seen"]], on="protein", how="left")
    d = a[a.valid.fillna(False)].dropna(subset=CONT + ["alt_seen"]).copy()
    for c in CONT:
        d["z_" + c] = (d[c] - d[c].mean()) / d[c].std()
    d["alt_seen_i"] = d.alt_seen.astype(int)
    d["active_site_i"] = d.active_site.astype(int)
    return a, d

def fit(d, formula):
    return sm.GEE.from_formula(formula, groups="pos", data=d, family=sm.families.Binomial(),
                               cov_struct=sm.cov_struct.Exchangeable()).fit()

if __name__ == "__main__":
    a, d = build()
    print(f"stats_A sample {len(a)} changes / {a.pos.nunique()} positions; "
          f"model sample {len(d)} changes / {d.pos.nunique()} positions "
          f"(positions whose alignment column fails the validity threshold are dropped)")

    d[CONT + ["active_site_i", "alt_seen_i"]].astype(float).corr(method="spearman") \
        .round(4).to_csv(T / "stats_G_feature_correlation.csv")

    full = "y ~ z_rsa + z_dist_ligand + active_site_i + z_conservation + z_blosum62 + alt_seen_i"
    stru = "y ~ z_rsa + z_dist_ligand + active_site_i"
    mf, ms = fit(d, full), fit(d, stru)

    G = pd.DataFrame({"term": mf.params.index, "coef": mf.params.values,
                      "OR": np.exp(mf.params.values),
                      "OR_ci95_low": np.exp(mf.conf_int()[0].values),
                      "OR_ci95_high": np.exp(mf.conf_int()[1].values),
                      "p": mf.pvalues.values})
    G["scale"] = ["per 1 SD" if t.startswith("z_") else "" for t in G.term]
    G.round(4).to_csv(T / "stats_G_gee_amenability.csv", index=False)
    print(G.round(3).to_string(index=False))

    names = list(mf.params.index)
    added = ["z_conservation", "z_blosum62", "alt_seen_i"]
    L = np.array([[1.0 if n == t else 0.0 for n in names] for t in added])
    w = mf.wald_test(L, scalar=True)
    cmp = pd.DataFrame([{"comparison": "full vs structure-only",
                         "added_terms": "+".join(added),
                         "wald_chi2": float(w.statistic), "df": len(added), "p": float(w.pvalue),
                         "QIC_full": float(mf.qic()[0]), "QIC_structure_only": float(ms.qic()[0])}])
    cmp.round(4).to_csv(T / "stats_G_model_comparison.csv", index=False)
    print(cmp.round(3).to_string(index=False))
