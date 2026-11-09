import pandas as pd, pytest
from gla import reference
from gla.hgvs import UTR5_LEN, UTR3_LEN, EXON_ENDS, EXON_STARTS
from gla.common import PROCESSED
from gla.cpg import possible_snvs

REF = reference.load()

def test_reference_constants_match_genbank():
    c = REF.self_check()
    assert (c["utr5"], c["utr3"]) == (UTR5_LEN, UTR3_LEN)
    assert c["exon_ends_c"] == EXON_ENDS
    assert [e + 1 for e in EXON_ENDS] == EXON_STARTS
    assert c["protein_len"] == 429

def test_glycosylation_sites_are_asn():
    assert all(REF.protein_gb[i - 1] == "N" for i in REF.glyco_codons)

def test_possible_snv_space():
    P = possible_snvs(REF)
    assert len(P) == 3 * 1290
    assert P.drop_duplicates("c_pos").cpg_site.sum() == 2 * REF.cds.count("CG")      # no CpG crosses a junction in this CDS
    assert not P.junction_ambiguous.any()
    assert P.cpg_transition.sum() == P.drop_duplicates("c_pos").cpg_site.sum()       # exactly one transition per CpG base

needs = pytest.mark.skipif(not (PROCESSED / "features_structural.parquet").exists(), reason="run stage 4 first")

@needs
def test_structure_table():
    T = pd.read_parquet(PROCESSED / "features_structural.parquet").set_index("protein_pos")
    assert len(T) == 429 and T.index.is_unique
    assert T.loc[1:31, "structure_resolved"].sum() == 0 and T.loc[32:421, "structure_resolved"].all()
    assert T.loc[[170, 231], "dist_ligand"].max() < 4.0           # catalytic Asp residues contact the bound ligand
    assert T.loc[[170, 231], "active_site"].all()
    assert abs(T.rsa_chainA.corr(T.rsa_chainB)) > 0.95            # two chains of the same dimer agree
    assert (T.rsa.dropna() >= 0).all() and T.rsa.dropna().max() < 1.6
