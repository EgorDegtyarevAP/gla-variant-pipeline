"""Stage 4 (level 2). IN: raw/3GXT.cif, raw/NM_000169.3.gb. OUT: processed/features_structural.parquet (+csv), one row per protein position 1..429.
Columns: rsa (dimer, chain mean), rsa_monomer, dsasa_interface, interface, dist_ligand (min heavy-atom distance to bound DGJ, A),
active_site (<= config distance), glyco_site (GenBank features), signal_peptide, structure_resolved, bfactor_ca."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from gla import reference
from gla.common import PROCESSED, RAW, log_source
from gla.structure import build_residue_table

if __name__ == "__main__":
    log_source(RAW / "3GXT.cif", "RCSB PDB 3GXT", "as downloaded by the user", "https://www.rcsb.org/structure/3GXT")
    log_source(RAW / "NM_000169.3.gb", "NCBI RefSeq NM_000169.3", "GenBank record dated 10-SEP-2026", "https://www.ncbi.nlm.nih.gov/nuccore/NM_000169.3")
    ref = reference.load()
    T = build_residue_table(ref)
    T.to_parquet(PROCESSED / "features_structural.parquet"); T.to_csv(PROCESSED / "features_structural.csv", index=False)
    r = T[T.structure_resolved]
    print("positions:", len(T), "resolved:", len(r), "active-site residues:", r[r.active_site].protein_pos.tolist())
    print("interface residues:", int(r.interface.sum()), " buried (RSA<0.2):", int((r.rsa < 0.2).sum()), " glyco:", T[T.glyco_site].protein_pos.tolist())
    print("catalytic D170/D231 dist to ligand:", r.set_index("protein_pos").loc[[170, 231], "dist_ligand"].round(2).tolist())
    print("chain A vs B RSA correlation:", round(r.rsa_chainA.corr(r.rsa_chainB), 3))
