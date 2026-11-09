import pandas as pd, pytest
from gla.hgvs import norm_c, ivs_to_c, parse_protein, classify, split_alternatives, proteins_compatible
from gla.common import RAW

def test_ivs_matches_amicus_own_conversions():
    """Amicus excluded file carries both IVS and c. spellings: independent check of the exon-boundary table."""
    assert ivs_to_c("IVS1-2A>G") == "c.195-2A>G"
    assert ivs_to_c("IVS6-10G>A") == "c.1000-10G>A"
    assert ivs_to_c("IVS4-1G>T") == "c.640-1G>T"

@pytest.mark.parametrize("raw,key", [
    ("c.18delA", "c.18del"), ("c.123_126dupCATG", "c.123_126dup"), ("c.194G>C (putative splicing site)", "c.194G>C"),
    ("c.134_138delTGCACinsGCTCG", "c.134_138delinsGCTCG"), ("IVS6-2A>G", "c.1000-2A>G"), ("c.7C>G", "c.7C>G")])
def test_norm_c(raw, key):
    assert norm_c(raw) == key

@pytest.mark.parametrize("raw", ["c.1209AAG[1]", "c.70T>C or c.70T>A", "g.?_?del", "c.?_?dup", None, float("nan")])
def test_norm_c_rejects(raw):
    assert norm_c(raw) is None

def test_split_alternatives():
    assert split_alternatives("c.126G>A or c.126G>C or c.126G>T") == ["c.126G>A", "c.126G>C", "c.126G>T"]

@pytest.mark.parametrize("raw,p1", [("p.Gly35Arg", "G35R"), ("L3V", "L3V"), ("W340X", "W340*"), ("p.Asn336=", "N336N"),
                                    ("p.Asp427Asp", "D427D"), ("p.Gly334Ter", "G334*"), ("X430Q", "*430Q"), ("UNKNOWN", "")])
def test_parse_protein(raw, p1):
    assert parse_protein(raw) == p1

@pytest.mark.parametrize("c,p1,cls", [
    ("c.100A>C", "N34H", "missense"), ("c.1019G>A", "W340*", "nonsense"), ("c.1008C>T", "N336N", "synonymous"),
    ("c.1A>G", "M1V", "start_lost"), ("c.194+64T>G", "", "intronic"), ("c.1000-1G>A", "", "splice_site"),
    ("c.-12G>A", "", "utr5"), ("c.*5T>C", "", "utr3"), ("c.1021del", "E341fs", "frameshift"),
    ("c.1055_1057dup", "A352_M353insT", "inframe_indel"), ("c.1116_1117delinsAT", "G373C", "mnv_missense"),
    ("c.539_547+9del", "", "splice_spanning_indel"), ("c.1288T>C", "*430Q", "stop_lost"), ("c.500G>A", "", "cds_snv_unannotated")])
def test_classify(c, p1, cls):
    assert classify(c, p1) == cls

def test_protein_compat():
    assert proteins_compatible("E341fs", "E341Nfs*6")
    assert proteins_compatible("A348fs", "L347fs*27")          # repeat shifted by tools
    assert not proteins_compatible("G35R", "G35S")
    assert not proteins_compatible("L106fs", "S102fs*16")
