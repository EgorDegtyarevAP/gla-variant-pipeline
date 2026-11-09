"""Stage 7b. Writes results/figures/fig6_structure.pml: colours 3GXT residues by the non-amenable share among Amicus-tested
protein changes at that position (blue = all amenable, orange = all non-amenable, grey = no S_amen data) and marks ClinVar P/LP positions
with sticks. NOT rendered in the build sandbox (no PyMOL): open the .pml in PyMOL with 3GXT.cif in data/raw and check the result by eye."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, pandas as pd
from gla.common import PROCESSED, RESULTS, RAW

LO, HI = np.array([0x2a, 0x78, 0xd6]), np.array([0xeb, 0x68, 0x34])
if __name__ == "__main__":
    M = pd.read_parquet(PROCESSED / "variants_master.parquet")
    pos = lambda s: s.str.extract(r"(\d+)")[0].astype(int)
    a = M[M.S_amen].drop_duplicates("protein").assign(pos=lambda d: pos(d.protein))
    share = a.groupby("pos").amenability_status.apply(lambda s: (s == "non_amenable").mean())
    plp = sorted(set(pos(M[M.S_path].protein)))
    L = ["load ../../data/raw/3GXT.cif, 3GXT", "hide everything", "show cartoon, 3GXT and polymer", "color grey70, 3GXT and polymer", "bg_color white",
         "set cartoon_transparency, 0.0", "show sticks, resn NOJ", "color black, resn NOJ"]
    for i in range(11):
        c = (LO + (HI - LO) * i / 10).astype(int)
        L.append(f"set_color gla_bin{i}, [{c[0]/255:.3f}, {c[1]/255:.3f}, {c[2]/255:.3f}]")
    for p, v in share.items():
        L.append(f"color gla_bin{int(round(v * 10))}, 3GXT and polymer and resi {p}")
    L.append("select plp_pos, 3GXT and polymer and resi " + "+".join(map(str, plp)))
    L += ["show sticks, plp_pos and not name N+C+O", "set stick_radius, 0.12", "orient", "# check by eye: colours = non-amenable share per position"]
    (RESULTS / "figures" / "fig6_structure.pml").write_text("\n".join(L) + "\n")
    print("written", len(share), "coloured positions,", len(plp), "P/LP positions")
