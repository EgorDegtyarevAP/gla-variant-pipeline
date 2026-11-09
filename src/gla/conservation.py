"""Level 3 - evolutionary conservation from the codon-aware ortholog alignment.

Pipeline: (1) rows are mapped to the 429 human codons (alignment columns where the human row has a base; columns where the
human row has a gap are insertions and are dropped); (2) codons are translated, any codon with a gap or non-ACGT is 'missing';
(3) rows are filtered by coverage/identity, one row per species; (4) per-position metrics with Henikoff sequence weights
(dense mammal sampling would otherwise dominate) and, as sensitivity, unweighted; the human row itself is excluded from
the frequencies, so `frac_human_aa` is the share of ORTHOLOGS carrying the human residue.
"""
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import substitution_matrices

AA = "ACDEFGHIKLMNPQRSTVWY"
IDX = {a: i for i, a in enumerate(AA)}
BLOSUM62 = substitution_matrices.load("BLOSUM62")


def load_alignment(path):
    recs = list(SeqIO.parse(path, "fasta"))
    human = [r for r in recs if "_Hsap" in r.id]
    assert len(human) == 1, "exactly one human row expected"
    return recs, human[0]


def codon_columns(human_row):
    cols = [i for i, c in enumerate(str(human_row.seq).upper()) if c != "-"]
    assert len(cols) % 3 == 0
    cc = [cols[i:i + 3] for i in range(0, len(cols), 3)]
    assert all(c[2] - c[0] == 2 for c in cc), "alignment is not codon-contiguous on the human row"
    return cc


def translate_row(seq, codon_cols):
    seq = str(seq).upper(); out = []
    for cc in codon_cols:
        cod = "".join(seq[i] for i in cc)
        out.append("-" if any(c not in "ACGT" for c in cod) else str(Seq(cod).translate()))
    return out


def species_of(rec_id):
    m = re.search(r"_([A-Za-z]+)/", rec_id)
    return m.group(1) if m else rec_id


def qc_table(recs, human, codon_cols, ref_protein, cfg):
    hp = translate_row(human.seq, codon_cols)
    hp_str = "".join(hp)
    assert hp_str == ref_protein[:len(hp_str)], "human alignment row does not translate to NP_000160 / NM_000169.3 CDS"
    rows, prot = [], {}
    for r in recs:
        p = translate_row(r.seq, codon_cols); prot[r.id] = p
        ok = [a != "-" for a in p]
        cov = float(np.mean(ok))
        idn = float(np.mean([a == b for a, b, k in zip(p, hp, ok) if k])) if any(ok) else np.nan
        rows.append({"id": r.id, "species": species_of(r.id), "coverage": cov, "identity": idn, "internal_stops": sum(a == "*" for a in p),
                     "is_human": r is human})
    q = pd.DataFrame(rows)
    q["score"] = q.coverage * q.identity
    q["status"] = "kept"
    q.loc[q.is_human, "status"] = "human (reference, excluded from frequencies)"
    q.loc[(q.status == "kept") & (q.internal_stops > 0), "status"] = "dropped: internal stop"
    q.loc[(q.status == "kept") & (q.coverage < cfg["min_coverage"]), "status"] = "dropped: coverage"
    q.loc[(q.status == "kept") & (q.identity < cfg["min_identity"]), "status"] = "dropped: identity (possible paralog / too divergent)"
    if cfg["one_per_species"]:
        for sp, g in q[q.status == "kept"].groupby("species"):
            if len(g) > 1:
                best = g.score.idxmax()
                q.loc[g.index.difference([best]), "status"] = "dropped: duplicate species (lower identity x coverage)"
    return q, prot


def henikoff_weights(mat):
    """mat: (n_seq, L) array of residue indices, -1 = missing. Henikoff & Henikoff 1994 position-based weights."""
    n, L = mat.shape; w = np.zeros(n)
    for j in range(L):
        col = mat[:, j]; ok = col >= 0
        if ok.sum() == 0:
            continue
        vals, counts = np.unique(col[ok], return_counts=True)
        k = len(vals); cnt = dict(zip(vals, counts))
        w[ok] += np.array([1.0 / (k * cnt[v]) for v in col[ok]])
    return w / w.sum()


def position_metrics(prot, kept_ids, human_prot, cfg, weighting="henikoff"):
    L = len(human_prot)
    mat = np.array([[IDX.get(a, -1) for a in prot[i]] for i in kept_ids])
    w = henikoff_weights(mat) if weighting == "henikoff" else np.full(len(kept_ids), 1 / len(kept_ids))
    rows = []
    for j in range(L):
        col = mat[:, j]; ok = col >= 0
        wj = w * ok
        cover = wj.sum()                                  # weighted share of orthologs covering the position
        if cover == 0:
            rows.append({"protein_pos": j + 1, "aa": human_prot[j], "n_orthologs": 0, "nongap_weight": 0.0})
            continue
        f = np.bincount(col[ok], weights=w[ok], minlength=20); f = f / f.sum()
        nz = f[f > 0]; H = float(-(nz * np.log2(nz)).sum())
        hum = IDX.get(human_prot[j], -1)
        rows.append({"protein_pos": j + 1, "aa": human_prot[j], "n_orthologs": int(ok.sum()), "nongap_weight": float(cover),
                     "entropy_bits": H, "conservation": 1 - H / np.log2(20), "frac_human_aa": float(f[hum]) if hum >= 0 else np.nan,
                     "top_aa": AA[int(np.argmax(f))], "freq": f})
    d = pd.DataFrame(rows)
    d["valid"] = d.nongap_weight >= cfg["min_nongap_weight"]
    return d


def substitution_features(pos_df, protein_changes):
    """Variant-level features for one-letter missense strings like 'G35R'."""
    P = pos_df.set_index("protein_pos")
    out = []
    for pc in protein_changes:
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", pc)
        ref, pos, alt = m[1], int(m[2]), m[3]
        r = P.loc[pos]
        f = r.get("freq")
        out.append({"protein": pc, "blosum62": float(BLOSUM62[ref][alt]),
                    "freq_alt_orthologs": float(f[IDX[alt]]) if isinstance(f, np.ndarray) and alt in IDX else np.nan,
                    "alt_seen": bool(f[IDX[alt]] > 0) if isinstance(f, np.ndarray) and alt in IDX else np.nan})
    return pd.DataFrame(out)
