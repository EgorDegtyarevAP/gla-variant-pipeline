"""Invariants of the built master table. Run after `make stage1 stage2`."""
import re, pandas as pd, pytest
from gla.common import PROCESSED, INTERIM
from gla.harmonize import build_mapper

pytestmark = pytest.mark.skipif(not (PROCESSED / "variants_master.parquet").exists(), reason="run stage 2 first")
M = pd.read_parquet(PROCESSED / "variants_master.parquet") if (PROCESSED / "variants_master.parquet").exists() else None

def test_keys_unique_and_wellformed():
    assert M.hgvs_c.is_unique and M.variant_key.is_unique
    assert M.hgvs_c.str.match(r"^c\.[-*]?\d+([+-]\d+)?").all()
    assert not M.hgvs_c.str.contains(r"\.\d+\.0|\.0[ACGT]").any(), "float formatting leaked into a key"

def test_snv_ref_alt_differ():
    snv = M.hgvs_c.str.extract(r"([ACGT])>([ACGT])$").dropna()
    assert (snv[0] != snv[1]).all()

def test_no_status_conflict_and_na_not_nonamenable():
    assert set(M.amenability_status.dropna()) <= {"amenable", "non_amenable"}
    assert M.loc[~M.in_amicus, "amenability_status"].isna().all()      # absent from table => NA, never non_amenable

def test_samples_nested():
    assert M.S_amen.sum() <= M.S_mis.sum() <= M.S_all.sum()
    assert (M.loc[M.S_amen, "variant_class"] == "missense").all()
    assert M.loc[M.S_path, "clinvar_class"].isin(["P", "LP", "P/LP"]).all()

def test_start_codon_ref_bases():
    """c.1-3 is ATG by definition: any SNV there must have the right reference base."""
    for pos, ref in ((1, "A"), (2, "T"), (3, "G")):
        s = M.hgvs_c[M.hgvs_c.str.match(rf"^c\.{pos}[ACGT]>")]
        assert s.str.extract(rf"^c\.{pos}([ACGT])>")[0].eq(ref).all()

def test_mapper_reproduces_gnomad():
    gn = pd.read_parquet(INTERIM / "gnomad.parquet")
    _, n, bad = build_mapper(gn)
    assert n > 700 and bad == []

def test_gnomad_clinvar_ids_agree():
    gn = pd.read_parquet(INTERIM / "gnomad.parquet").dropna(subset=["gnomad_clinvar_id"])
    cv = pd.read_parquet(INTERIM / "clinvar.parquet").dropna(subset=["hgvs_c"]).set_index("clinvar_id").hgvs_c
    k = gn.gnomad_clinvar_id.astype(int).map(cv)
    assert (k.dropna() == gn.hgvs_c[k.notna()]).all()
