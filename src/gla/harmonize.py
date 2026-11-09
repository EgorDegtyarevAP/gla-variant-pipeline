"""Stage 2: map to one key (hgvs_c on NM_000169), merge the sources, run integrity checks."""
import re
import numpy as np
import pandas as pd

from . import reference
from .genomic import GlaMapper, COMP
from .hgvs import classify, region, proteins_compatible

EXONS_CDNA = [(1, 194), (195, 369), (370, 547), (548, 639), (640, 801), (802, 999), (1000, 1290)]
AMEN_CRITERION = ("Amicus Galafold amenability table (galafold_amenability.csv); assay and thresholds "
                  "not stated in the file - TO FILL from the table's source before submission")


def build_mapper(gn):
    snv = gn[(gn.g_ref.str.len() == 1) & (gn.g_alt.str.len() == 1)]
    mp = GlaMapper.from_gnomad(list(zip(snv.g_pos, snv.hgvs_c_raw)), EXONS_CDNA)
    bad = [(r.hgvs_c_raw, mp.snv_to_c(r.g_pos, r.g_ref, r.g_alt)) for r in snv.itertuples()
           if mp.snv_to_c(r.g_pos, r.g_ref, r.g_alt) != r.hgvs_c_raw]
    return mp, len(snv), bad


def observed_utr_extent(frames):
    """Largest |c.-N| and c.*N seen in any source (used to separate transcript from flank)."""
    d5 = d3 = 0
    for f in frames:
        for k in f["hgvs_c"].dropna():
            m = re.match(r"^c\.(-|\*)(\d+)", k)
            if m:
                if m[1] == "-":
                    d5 = max(d5, int(m[2]))
                else:
                    d3 = max(d3, int(m[2]))
    return d5, d3


def resolve_unparsed(unp, cv, mp, d5, d3):
    """Map ClinVar genomic-only names to c.; classify as in-transcript / flank / unresolved."""
    out = []
    for r in unp.itertuples():
        base = {"clinvar_id": r.clinvar_id, "name_g": r.name_g}
        if r.g_kind != "snv":
            out.append({**base, "hgvs_c": None, "status": "unresolved_non_snv"}); continue
        key = mp.snv_to_c(r.g_pos, r.g_ref, r.g_alt)
        m = re.match(r"^c\.(-|\*)(\d+)", key)
        flank = bool(m) and ((m[1] == "-" and int(m[2]) > d5) or (m[1] == "*" and int(m[2]) > d3))
        out.append({**base, "hgvs_c": key, "status": "outside_transcript_extent" if flank else "mapped"})
    return pd.DataFrame(out)


_SNV_RE = re.compile(r"^c\.(?P<b>[-*]?\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")


def _cds_annotation(ref, c):
    """Reference-base check (CDS and mRNA UTR) and codon-level protein effect (CDS only)."""
    out = {"ref_ok": None, "protein_cds": None, "codon_no": None, "codon_pos": None}
    m = _SNV_RE.match(c)
    if not m:
        return out
    b = m["b"]
    base = ref.base_at(b)
    if base is not None:
        out["ref_ok"] = base == m["ref"]
    if not (b.startswith("-") or b.startswith("*")) and 1 <= int(b) <= len(ref.cds):
        e = ref.snv_effect(int(b), m["ref"], m["alt"])
        out.update(protein_cds=e["protein"] if e["ref_ok"] else None, codon_no=e["codon_no"], codon_pos=e["codon_pos"])
    return out


def _first(s):
    s = s.dropna()
    return s.iloc[0] if len(s) else np.nan


def build_master(cv, gn, am, unp_res):
    audit = {}
    # --- ClinVar: add mapped genomic-only records
    cv = cv.copy()
    cv["clinvar_hgvs_source"] = np.where(cv.hgvs_c.notna(), "c.change", "unparsed")
    for r in unp_res[unp_res.status == "mapped"].itertuples():
        if r.clinvar_id in set(cv.clinvar_id):
            cv.loc[cv.clinvar_id == r.clinvar_id, ["hgvs_c", "clinvar_hgvs_source"]] = [r.hgvs_c, "mapped_from_g"]
        else:
            cv = pd.concat([cv, pd.DataFrame([{"clinvar_id": r.clinvar_id, "hgvs_c": r.hgvs_c,
                                              "clinvar_class": "annotation_missing",
                                              "clinvar_hgvs_source": "mapped_from_g"}])], ignore_index=True)
    unresolved_cv = cv[cv.hgvs_c.isna()]
    audit["clinvar_unresolved"] = unresolved_cv[["clinvar_id", "hgvs_c_raw"]]
    cvk = cv.dropna(subset=["hgvs_c"])
    audit["clinvar_dup_keys"] = int(cvk.hgvs_c.duplicated().sum())
    cvk = cvk.sort_values(["clinvar_stars"], ascending=False).drop_duplicates("hgvs_c")
    gnk = gn.dropna(subset=["hgvs_c"])
    audit["gnomad_dup_keys"] = int(gnk.hgvs_c.duplicated().sum())

    # --- Amicus: status per key, conflicts
    amg = am.groupby("hgvs_c").agg(
        amenability_status=("amenability_status", lambda s: s.iloc[0] if s.nunique() == 1 else "CONFLICT"),
        amicus_ids=("amicus_id", lambda s: ";".join(s)),
        amicus_variant_type=("amicus_variant_type", "first"),
        amicus_n_alt=("n_alternatives_in_row", "max"),
        amicus_protein=("protein", _first)).reset_index()
    audit["amicus_conflicts"] = amg[amg.amenability_status == "CONFLICT"]

    m = cvk[["hgvs_c", "clinvar_id", "protein", "clinvar_significance", "clinvar_class", "clinvar_review",
             "clinvar_stars", "clinvar_phenotype", "clinvar_last_evaluated", "rsid", "clinvar_origin",
             "clinvar_hgvs_source"]].rename(columns={"protein": "protein_clinvar"})
    g = gnk[["hgvs_c", "gnomad_id", "g_pos", "g_ref", "g_alt", "protein", "gnomad_vep", "gnomad_ac", "gnomad_an",
             "gnomad_af", "gnomad_hom", "gnomad_hemi", "gnomad_pass", "gnomad_flags", "gnomad_faf_group",
             "gnomad_faf", "gnomad_clinvar_id", "gnomad_clinvar_class", "cadd", "revel", "sift", "polyphen", "spliceai", "pangolin",
             "phylop"]].rename(columns={"protein": "protein_gnomad"})
    a = amg.rename(columns={"amicus_protein": "protein_amicus"})
    M = m.merge(g, on="hgvs_c", how="outer").merge(a, on="hgvs_c", how="outer")
    for c in ("clinvar_id", "gnomad_clinvar_id"):
        M[c] = M[c].astype("Int64")
    M["in_clinvar"] = M.clinvar_id.notna()
    M["in_gnomad"] = M.gnomad_id.notna()
    M["in_amicus"] = M.amicus_ids.notna()
    M["sources"] = (M.in_clinvar.map({True: "clinvar;", False: ""}) + M.in_gnomad.map({True: "gnomad;", False: ""})
                    + M.in_amicus.map({True: "amicus;", False: ""})).str.rstrip(";")
    M["n_sources"] = M[["in_clinvar", "in_gnomad", "in_amicus"]].sum(axis=1)

    # --- reference check + CDS-derived protein effect for every SNV that lies in the mRNA
    ref = reference.load()
    ann = M.hgvs_c.map(lambda c: _cds_annotation(ref, c)).apply(pd.Series)
    for col in ann.columns:
        M[col] = ann[col]
    M["ref_ok"] = M["ref_ok"].astype("boolean")

    # --- protein: CDS-derived for SNVs in the CDS (authoritative); else priority ClinVar > gnomAD > Amicus
    def pick(r):
        if isinstance(r["protein_cds"], str):
            return r["protein_cds"]
        for c in ("protein_clinvar", "protein_gnomad", "protein_amicus"):
            if isinstance(r[c], str) and r[c]:
                return r[c]
        return ""
    M["protein"] = M.apply(pick, axis=1)
    rep = M[["protein_clinvar", "protein_gnomad", "protein_amicus"]].fillna("")
    M["protein_conflict"] = [
        not all(proteins_compatible(x, ref_p) for x in row)
        for row, ref_p in zip(rep.itertuples(index=False), M.protein.fillna(""))]
    audit["protein_conflicts"] = M[M.protein_conflict][["hgvs_c", "protein", "protein_clinvar", "protein_gnomad", "protein_amicus"]]
    audit["ref_mismatch"] = M[M.ref_ok.eq(False)][["hgvs_c", "sources"]]

    M["variant_class"] = [("ref_mismatch" if ok is False else classify(c, p))
                          for c, p, ok in zip(M.hgvs_c, M.protein, M.ref_ok.astype(object).where(M.ref_ok.notna(), None))]
    M["region"] = M.hgvs_c.map(region)
    M["amenability_status"] = M.amenability_status.replace("CONFLICT", np.nan)
    M["amenability_criterion"] = np.where(M.in_amicus, AMEN_CRITERION, None)
    M["variant_key"] = "GLA:" + M.hgvs_c

    # --- samples (docs/samples.md)
    M["S_all"] = True
    M["S_mis"] = M.variant_class.eq("missense")
    M["S_amen"] = M.S_mis & M.amenability_status.isin(["amenable", "non_amenable"])
    M["S_path"] = M.S_mis & M.clinvar_class.isin(["P", "LP", "P/LP"])
    return M.sort_values("hgvs_c").reset_index(drop=True), audit
