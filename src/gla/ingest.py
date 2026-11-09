"""Stage 1: read each raw source, rename columns, type them, add only *source-local* derived fields.
No cross-source logic here (that is stage 2)."""
import re
import numpy as np
import pandas as pd

from .common import RAW, RAW_FILES
from .hgvs import norm_c, parse_protein, split_alternatives, ivs_to_c

CANONICAL_ENST = "ENST00000218516.4"   # GLA canonical transcript in gnomAD v4.1.1 (MANE Select NM_000169.3 assumed; verify)

CLINVAR_CLASS = {
    "Pathogenic": "P", "Likely pathogenic": "LP", "Pathogenic/Likely pathogenic": "P/LP",
    "Uncertain significance": "VUS", "Likely benign": "LB", "Benign": "B", "Benign/Likely benign": "B/LB",
    "Conflicting classifications of pathogenicity": "conflicting",
}
STARS = {
    "practice guideline": 4, "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1, "criteria provided, conflicting classifications": 1,
    "no assertion criteria provided": 0,
}


def clinvar_class(sig):
    if sig in CLINVAR_CLASS:
        return CLINVAR_CLASS[sig]
    if isinstance(sig, str) and ";" in sig:
        return "conflicting"
    return "not_classified"


def read_clinvar():
    d = pd.read_csv(RAW / RAW_FILES["clinvar"], dtype=str)
    out = pd.DataFrame({
        "clinvar_id": d["variant_id"].astype(int),
        "hgvs_c_raw": d["c.change"],
        "hgvs_p_raw": d["p.change"],
        "clinvar_significance": d["clinical_significance"].replace("-", np.nan),
        "clinvar_review": d["review_status"].replace("-", np.nan),
        "clinvar_phenotype": d["phenotype"],
        "clinvar_last_evaluated": pd.to_datetime(d["last_evaluated"], format="%b %d, %Y", errors="coerce"),
        "rsid": d["rsid"].replace("-1", np.nan),
        "clinvar_origin": d["origin"],
    })
    out["clinvar_class"] = out["clinvar_significance"].map(clinvar_class)
    out["clinvar_stars"] = out["clinvar_review"].map(STARS).fillna(0).astype(int)
    out["hgvs_c"] = out["hgvs_c_raw"].map(norm_c)
    out["protein"] = out["hgvs_p_raw"].map(parse_protein)
    return out


def read_clinvar_unparsed():
    d = pd.read_csv(RAW / RAW_FILES["clinvar_unparsed"])
    d = d.rename(columns={"VariationID": "clinvar_id", "Name": "name_g"})
    m = d["name_g"].str.extract(r"^NC_000023\.11:g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
    d["g_pos"] = pd.to_numeric(m["pos"])
    d["g_ref"], d["g_alt"] = m["ref"], m["alt"]
    d["g_kind"] = np.where(m["pos"].notna(), "snv", "complex_or_unresolved")
    return d


def read_gnomad():
    d = pd.read_csv(RAW / RAW_FILES["gnomad"])
    n_all = len(d)
    other = d[d["Transcript"] != CANONICAL_ENST]
    c = d[d["Transcript"] == CANONICAL_ENST].copy()
    out = pd.DataFrame({
        "gnomad_id": c["gnomAD ID"],
        "g_pos": c["Position"], "g_ref": c["Reference"], "g_alt": c["Alternate"],
        "hgvs_c_raw": c["Transcript Consequence"],
        "hgvs_p_raw": c["Protein Consequence"],
        "gnomad_vep": c["VEP Annotation"],
        "gnomad_ac": c["Allele Count"], "gnomad_an": c["Allele Number"], "gnomad_af": c["Allele Frequency"],
        "gnomad_hom": c["Homozygote Count"], "gnomad_hemi": c["Hemizygote Count"],
        "gnomad_filter_joint": c["Filters - joint"], "gnomad_flags": c["Flags"],
        "gnomad_faf_group": c["GroupMax FAF group"], "gnomad_faf": c["GroupMax FAF frequency"],
        "gnomad_clinvar_id": c["ClinVar Variation ID"],
        "gnomad_clinvar_class": c["ClinVar Germline Classification"],
        "cadd": c["cadd"], "revel": c["revel_max"], "sift": c["sift_max"], "polyphen": c["polyphen_max"],
        "spliceai": c["spliceai_ds_max"], "pangolin": c["pangolin_largest_ds"], "phylop": c["phylop"],
    })
    out["gnomad_pass"] = out["gnomad_filter_joint"].eq("PASS")
    out["hgvs_c"] = out["hgvs_c_raw"].map(norm_c)
    out["protein"] = out["hgvs_p_raw"].map(parse_protein)
    meta = {"rows_total": n_all, "rows_canonical": len(out), "rows_other_transcripts": len(other),
            "other_transcripts": other["Transcript"].value_counts().to_dict()}
    return out, meta


def read_amicus():
    """One output row per *nucleotide change*. Rows listing alternatives ('c.70T>C or c.70T>A') are expanded."""
    main = pd.read_csv(RAW / RAW_FILES["amicus"])
    excl = pd.read_csv(RAW / RAW_FILES["amicus_excluded"])
    rows, unresolved = [], []
    for r in main.itertuples(index=False):
        base = {"amicus_id": r.variant_id, "amenability_status": r.amenability_status,
                "amicus_variant_type": r.variant_type, "amicus_protein_raw": r.protein_change,
                "amicus_legacy_nt": r.legacy_nt_change}
        if r.variant_id in set(excl["variant_id"]):
            unresolved.append({**base, "reason": f"excluded_{r.variant_type}_user_list"})
            continue
        parts = split_alternatives(r.legacy_nt_change)
        keys = [(p, norm_c(p)) for p in parts]
        for raw, key in keys:
            if key is None:
                unresolved.append({**base, "reason": "unparseable_or_unlocalised_nt_change", "amicus_nt_part": raw})
            else:
                rows.append({**base, "hgvs_c_raw": raw, "hgvs_c": key, "n_alternatives_in_row": len(parts),
                             "protein": parse_protein(r.protein_change)})
    # a row that listed several nt changes may yield duplicate normalised keys (e.g. delGAT or del)
    ok = pd.DataFrame(rows).drop_duplicates(["amicus_id", "hgvs_c"])
    return ok, pd.DataFrame(unresolved), len(main), len(excl)
