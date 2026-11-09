"""Level 2 - structural features from PDB 3GXT (human alpha-Gal A + 1-deoxygalactonojirimycin, 2.7 A, chains A/B = biological dimer).

Residue numbering in 3GXT (auth_seq_id 32..422) equals UniProt P06280 / HGVS p. numbering (verified against the
CDS translation in `build_residue_table`), so no offset is applied.
RSA = SASA(dimer, protein atoms only) / maxASA. maxASA: Tien et al. 2013 (theoretical scale) - verify against the paper before submission.
"""
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, Selection
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.SeqUtils import seq1

from .common import RAW

MAX_ASA = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0,
           "I": 197.0, "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0,
           "Y": 263.0, "V": 174.0}
LIGAND = "NOJ"                      # 1-deoxygalactonojirimycin in 3GXT
CHAINS = ("A", "B")


def load_model(path=None):
    return MMCIFParser(QUIET=True).get_structure("3GXT", path or RAW / "3GXT.cif")[0]


def _protein_only(model, chains):
    """New model holding copies of the standard amino-acid residues of `chains` (no water, ligand, glycan)."""
    m = Model(0)
    for cid in chains:
        c = Chain(cid)
        for r in model[cid]:
            if r.id[0] == " " and is_aa(r, standard=True):
                c.add(r.copy())
        m.add(c)
    return m


def _sasa_per_residue(model, chains, **kw):
    sr = ShrakeRupley(probe_radius=1.4, n_points=200)
    sub = _protein_only(model, chains)
    sr.compute(sub, level="R")
    return {(c.id, r.id[1]): r.sasa for c in sub for r in c}


def build_residue_table(ref, model=None, active_cutoff=4.5, interface_dsasa_min=1.0):
    model = model or load_model()
    prot = ref.protein_gb
    rows = []
    # SASA: dimer, and each chain alone
    dimer = _sasa_per_residue(model, CHAINS)
    mono = {cid: _sasa_per_residue(model, (cid,)) for cid in CHAINS}
    lig_atoms = {cid: [a for r in model[cid] if r.resname == LIGAND for a in r] for cid in CHAINS}
    all_lig = [a for cid in CHAINS for a in lig_atoms[cid]]
    bad_numbering = []
    for cid in CHAINS:
        for r in model[cid]:
            if r.id[0] != " " or not is_aa(r, standard=True):
                continue
            pos, aa = r.id[1], seq1(r.resname)
            if prot[pos - 1] != aa:
                bad_numbering.append((cid, pos, aa, prot[pos - 1]))
            d_lig = min(float(np.linalg.norm(a.coord - b.coord)) for a in r for b in all_lig)
            rows.append({"protein_pos": pos, "aa": aa, "chain": cid,
                         "sasa_dimer": dimer[(cid, pos)], "sasa_monomer": mono[cid][(cid, pos)],
                         "dist_ligand": d_lig, "bfactor_ca": r["CA"].bfactor})
    assert not bad_numbering, f"structure numbering != CDS translation: {bad_numbering[:5]}"
    d = pd.DataFrame(rows)
    d["rsa_dimer"] = d.sasa_dimer / d.aa.map(MAX_ASA)
    d["rsa_monomer"] = d.sasa_monomer / d.aa.map(MAX_ASA)
    d["dsasa_interface"] = d.sasa_monomer - d.sasa_dimer
    # one row per position: chain-averaged values, chain A/B kept for QC
    g = d.groupby(["protein_pos", "aa"])
    out = g.agg(n_chains=("chain", "nunique"), rsa=("rsa_dimer", "mean"), rsa_monomer=("rsa_monomer", "mean"),
                dsasa_interface=("dsasa_interface", "mean"), dist_ligand=("dist_ligand", "min"),
                bfactor_ca=("bfactor_ca", "mean")).reset_index()
    ab = d.pivot(index="protein_pos", columns="chain", values="rsa_dimer")
    out["rsa_chainA"] = out.protein_pos.map(ab["A"]); out["rsa_chainB"] = out.protein_pos.map(ab["B"])
    out["active_site"] = out.dist_ligand <= active_cutoff
    out["interface"] = out.dsasa_interface >= interface_dsasa_min
    out["structure_resolved"] = True
    # complete 1..429 so unresolved positions are explicit (signal peptide 1-31, C-terminal tail 423-429)
    full = pd.DataFrame({"protein_pos": range(1, len(prot) + 1)}); full["aa"] = [prot[i - 1] for i in full.protein_pos]
    full = full.merge(out.drop(columns="aa"), on="protein_pos", how="left")
    full["structure_resolved"] = full.structure_resolved.fillna(False).astype(bool)
    full["signal_peptide"] = full.protein_pos <= 31
    full["glyco_site"] = full.protein_pos.isin(ref.glyco_codons)
    full["active_site"] = full.active_site.fillna(False).astype(bool)
    full["interface"] = full.interface.fillna(False).astype(bool)
    return full
