"""Stage 5 (level 3). IN: raw/Human_GLA_orthologues.fa (codon-aware nt MSA), raw/NM_000169.3.gb.
OUT: processed/features_conservation.parquet (+csv; one row per protein position 1..429, weighted = main),
processed/features_conservation_unweighted.csv (sensitivity), processed/features_substitution.parquet (BLOSUM62, alt-allele frequency per missense protein change),
results/tables/conservation_qc.csv (every row: coverage, identity, kept/dropped + reason)."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, pandas as pd
from gla import reference, conservation as cons
from gla.common import PROCESSED, RAW, RESULTS, load_config, log_source

if __name__ == "__main__":
    cfg = load_config()["conservation"]; (RESULTS / "tables").mkdir(parents=True, exist_ok=True)
    f = RAW / cfg["alignment"]; log_source(f, "Ortholog codon alignment (user-supplied)", "source/release NOT RECORDED - to be filled (looks like Ensembl Compara IDs)")
    ref = reference.load()
    recs, human = cons.load_alignment(f); cc = cons.codon_columns(human)
    q, prot = cons.qc_table(recs, human, cc, ref.protein_gb, cfg)
    q.drop(columns="score").to_csv(RESULTS / "tables" / "conservation_qc.csv", index=False)
    kept = q[q.status == "kept"].id.tolist(); hp = cons.translate_row(human.seq, cc)
    W = cons.position_metrics(prot, kept, hp, cfg, "henikoff"); U = cons.position_metrics(prot, kept, hp, cfg, "unweighted")
    W.drop(columns="freq").to_parquet(PROCESSED / "features_conservation.parquet"); W.drop(columns="freq").to_csv(PROCESSED / "features_conservation.csv", index=False)
    U.drop(columns="freq").to_csv(PROCESSED / "features_conservation_unweighted.csv", index=False)
    M = pd.read_parquet(PROCESSED / "variants_master.parquet")
    S = cons.substitution_features(W, sorted(set(M[M.variant_class == "missense"].protein)))
    S.to_parquet(PROCESSED / "features_substitution.parquet")
    print(q.status.value_counts().to_string()); print("kept species:", len(kept), "| valid positions:", int(W.valid.sum()), "of", len(W))
    print("median conservation (weighted/unweighted):", round(W.conservation.median(), 3), round(U.conservation.median(), 3),
          "| frac positions with >=90% of orthologs carrying the human residue:", round((W.frac_human_aa >= 0.9).mean(), 3))
    print("invalid positions:", W[~W.valid].protein_pos.tolist()[:40])
