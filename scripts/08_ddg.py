"""Stage 4 (level 2), step 2: stability change (DDG).
DynaMut2 is a web server and cannot be reached from the build sandbox, so this step is split:
 (a) writes the submission list  results/ddg_input_dynamut2.csv  (chain A numbering = UniProt/HGVS numbering,
     one row per unique protein change in S_amen U S_path whose residue is resolved in 3GXT);
 (b) if data/raw/dynamut2_results.csv exists (columns: pos, wt, mut, ddg  - kcal/mol, as returned by the server;
     CHECK THE SIGN CONVENTION on the server page and record it in LOG.md), merges it into
     processed/features_ddg.parquet keyed by protein change.
Nothing is estimated locally; missing values stay NA."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from gla.common import PROCESSED, RAW, RESULTS

if __name__ == "__main__":
    RESULTS.mkdir(exist_ok=True)
    M = pd.read_parquet(PROCESSED / "variants_master.parquet")
    S = pd.read_parquet(PROCESSED / "features_structural.parquet").set_index("protein_pos")
    s = M[(M.S_amen | M.S_path)].copy()
    s = s[s.protein.str.match(r"^[A-Z]\d+[A-Z]$")]
    s["pos"] = s.protein.str.extract(r"(\d+)")[0].astype(int)
    s = s[s.pos.map(S.structure_resolved)]
    u = s.drop_duplicates("protein")[["protein", "pos"]].copy()
    u["wt"], u["mut"] = u.protein.str[0], u.protein.str[-1]
    u["chain"] = "A"
    u["dynamut2_line"] = u.chain + " " + u.wt + u.pos.astype(str) + u.mut
    u.sort_values("pos")[["chain", "wt", "pos", "mut", "protein", "dynamut2_line"]].to_csv(RESULTS / "ddg_input_dynamut2.csv", index=False)
    print(f"{len(u)} unique protein changes to submit (positions resolved in 3GXT)")
    f = RAW / "dynamut2_results.csv"
    if f.exists():
        d = pd.read_csv(f); d["protein"] = d.wt + d.pos.astype(str) + d.mut
        d[["protein", "ddg"]].drop_duplicates("protein").to_parquet(PROCESSED / "features_ddg.parquet")
        print("merged", len(d), "DynaMut2 results")
    else:
        print("data/raw/dynamut2_results.csv not found - DDG left NA")
