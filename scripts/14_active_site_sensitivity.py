"""Stage 6c. Sensitivity of the active-site category to the ligand-distance cut-off.
IN: processed/variants_master, features_structural, features_conservation, features_substitution.
OUT: results/tables/stats_H_active_site_sensitivity.csv

The 4.5 A cut-off is conventional, not derived, and the ligand it is measured from is migalastat
itself, so the category is not independent of the outcome it is compared against. This stage
repeats the active-site comparison and the multivariable model at 3.5, 4.5 and 6.0 A."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, pandas as pd, statsmodels.api as sm
from gla.common import PROCESSED, RESULTS, load_config
from gla import stats as st

T = RESULTS / "tables"
CUTOFFS = (3.5, 4.5, 6.0)
SEED = load_config()["stats"]["random_seed"]
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from importlib import import_module
build = import_module("13_multivariable").build

if __name__ == "__main__":
    a, d = build()
    S = pd.read_parquet(PROCESSED / "features_structural.parquet")
    M = pd.read_parquet(PROCESSED / "variants_master.parquet")
    R = S[S.structure_resolved].set_index("protein_pos")
    plp = M[M.S_path].copy(); plp["pos"] = plp.protein.str.extract(r"(\d+)")[0].astype(int)
    plp = plp.merge(R[["dist_ligand"]], left_on="pos", right_index=True, how="left")

    upos = np.sort(a.pos.unique()); idx = np.searchsorted(upos, a.pos.to_numpy())
    y = a.y.to_numpy().astype(bool)
    out = []
    for c in CUTOFFS:
        site = (R.dist_ligand <= c)
        fp = site.reindex(upos).to_numpy(float)
        obs, p = st.cluster_permutation(idx, y, fp, st._prop_diff, 10000, np.random.default_rng(SEED))
        lo, hi = st.cluster_bootstrap_ci(idx, y, fp, st._prop_diff, 2000, np.random.default_rng(SEED + 1))
        inn = plp.dist_ligand <= c
        sub = plp[inn].drop_duplicates("protein")
        with_status = sub[sub.amenability_status.notna()]
        d2 = d.copy(); d2["site_i"] = (d2.dist_ligand <= c).astype(int)
        m = sm.GEE.from_formula("y ~ z_rsa + z_dist_ligand + site_i + z_conservation + z_blosum62 + alt_seen_i",
                                groups="pos", data=d2, family=sm.families.Binomial(),
                                cov_struct=sm.cov_struct.Exchangeable()).fit()
        out.append({"cutoff_A": c,
                    "n_active_site_residues": int(site.sum()),
                    "share_of_resolved_positions_pct": round(100 * site.sum() / len(R), 1),
                    "n_PLP_missense_in_site": len(sub),
                    "n_PLP_missense_in_site_with_status": len(with_status),
                    "non_amenable_share_pct": round(100 * (with_status.amenability_status == "non_amenable").mean(), 1)
                                              if len(with_status) else float("nan"),
                    "diff_prop_non_minus_amen": obs, "ci95_low": lo, "ci95_high": hi, "p_cluster_perm": p,
                    "GEE_OR_active_site": float(np.exp(m.params["site_i"])),
                    "GEE_OR_ci95_low": float(np.exp(m.conf_int().loc["site_i", 0])),
                    "GEE_OR_ci95_high": float(np.exp(m.conf_int().loc["site_i", 1])),
                    "GEE_p_active_site": float(m.pvalues["site_i"]),
                    "GEE_OR_dist_ligand": float(np.exp(m.params["z_dist_ligand"])),
                    "GEE_p_dist_ligand": float(m.pvalues["z_dist_ligand"])})
    H = pd.DataFrame(out)
    H.to_csv(T / "stats_H_active_site_sensitivity.csv", index=False)
    print(H.round(4).to_string(index=False))
