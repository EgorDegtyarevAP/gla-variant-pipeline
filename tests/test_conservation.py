import pandas as pd
import pytest
from gla import reference
from gla.common import PROCESSED, RESULTS

@pytest.fixture(scope="module")
def C():
    return pd.read_parquet(PROCESSED / "features_conservation.parquet").set_index("protein_pos")

def test_covers_all_residues(C):
    assert list(C.index) == list(range(1, 430))

def test_ranges(C):
    v = C[C.valid]
    assert v.conservation.between(0, 1).all() and v.frac_human_aa.between(0, 1).all()
    assert len(v) > 380

def test_catalytic_residues_conserved(C):
    for pos in (170, 231):          # catalytic D170 / D231 (3GXT, <=3 A from ligand)
        assert C.loc[pos, "valid"] and C.loc[pos, "conservation"] >= 0.95

def test_human_row_is_reference_protein():
    ref = reference.load()
    from gla import conservation as cs
    from gla.common import RAW
    recs, hs = cs.load_alignment(RAW / "Human_GLA_orthologues.fa")
    cols = cs.codon_columns(hs)
    assert "".join(cs.translate_row(hs.seq, cols)) == ref.protein_gb[:len(cols)]
    assert len(cols) == 429

def test_qc_statuses():
    q = pd.read_csv(RESULTS / "tables" / "conservation_qc.csv")
    assert set(q.status) >= {"kept"} and (q.status == "kept").sum() > 100
