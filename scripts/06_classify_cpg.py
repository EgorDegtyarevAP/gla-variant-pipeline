"""Stage 3 (level 1: mutational context). IN: processed/variants_master.parquet, raw/NM_000169.3.gb.
OUT: processed/variants_cpg.parquet (hgvs_c -> CpG flags for CDS SNVs),
results/tables/cpg_possible_snvs.csv, results/tables/cpg_enrichment.csv.
Samples (population vs patient ascertainment; Amicus-only variants are NOT used: they were selected for in-vitro testing):
  gnomAD CDS SNVs (all consequences), ClinVar P/LP CDS SNVs (all consequences), S_path (missense P/LP),
  gnomAD missense, and ClinVar P/LP split by consequence."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from gla import reference
from gla.common import PROCESSED, RESULTS
from gla.cpg import possible_snvs, enrichment

if __name__ == "__main__":
    (RESULTS / 'tables').mkdir(parents=True, exist_ok=True)
    ref = reference.load()
    M = pd.read_parquet(PROCESSED / "variants_master.parquet")
    P = possible_snvs(ref)
    P.to_csv(RESULTS / "tables" / "cpg_possible_snvs.csv", index=False)
    M[["hgvs_c"]].merge(P[["hgvs_c", "cpg_site", "cpg_transition", "transition", "junction_ambiguous"]], on="hgvs_c", how="left") \
        .to_parquet(PROCESSED / "variants_cpg.parquet")

    cds_snv = M[M.hgvs_c.isin(set(P.hgvs_c))]
    gn = cds_snv[cds_snv.in_gnomad & cds_snv.gnomad_pass.fillna(False)]
    plp = cds_snv[cds_snv.clinvar_class.isin(["P", "LP", "P/LP"])]
    allc = {"missense", "synonymous", "nonsense", "start_lost", "stop_lost"}
    rows = [
        enrichment(gn.hgvs_c, P, allc, "gnomAD (PASS) CDS SNVs"),
        enrichment(gn.hgvs_c, P, {"missense"}, "gnomAD (PASS) missense"),
        enrichment(gn.hgvs_c, P, {"synonymous"}, "gnomAD (PASS) synonymous"),
        enrichment(plp.hgvs_c, P, allc, "ClinVar P/LP CDS SNVs"),
        enrichment(plp.hgvs_c, P, {"missense"}, "ClinVar P/LP missense (S_path)"),
        enrichment(plp.hgvs_c, P, {"nonsense"}, "ClinVar P/LP nonsense"),
    ]
    T = pd.DataFrame(rows); T.to_csv(RESULTS / "tables" / "cpg_enrichment.csv", index=False)
    pd.set_option("display.width", 220); print(T.drop(columns=["consequences"]).round(4).to_string(index=False))
    print("possible CDS SNVs:", len(P), " CpG sites (positions):", int(P.drop_duplicates('c_pos').cpg_site.sum()),
          " junction-ambiguous positions:", int(P.drop_duplicates('c_pos').junction_ambiguous.sum()))
