"""Figure 1. Map of GLA and alpha-Gal A.

Panels: (a) exon/domain track, (b) P/LP missense per residue by amenability
status, (c) positional conservation, (d) relative solvent accessibility.
"""
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

from gla.figures import style as S
from gla.figures.data import pos_of

NAME = "fig1_gene_protein_map"
FIGSIZE = (7.2, 6.4)
HEIGHT_RATIOS = [0.8, 2.3, 1.0, 1.0]
RSA_BURIED = 0.2
SIGNAL_PEPTIDE_END = 31


def _exon_track(ax, d):
    L = d.length
    ax.set_xlim(0, L + 1)
    ends = [(e - d.ref.utr5_len) for _, e in d.ref.exons_mrna[:-1]]
    bounds = [0] + [int(np.ceil(c / 3)) for c in ends] + [L]
    for i in range(len(bounds) - 1):
        ax.add_patch(Rectangle((bounds[i] + 0.5, 0.55), bounds[i + 1] - bounds[i], 0.35,
                               fc="#d9d8d3" if i % 2 else "#c3c2bb", ec=S.SURF, lw=0.8))
        ax.text((bounds[i] + bounds[i + 1]) / 2 + 0.5, 0.725, f"{i + 1}",
                ha="center", va="center", fontsize=7, color=S.INK)
    ax.add_patch(Rectangle((0.5, 0.2), SIGNAL_PEPTIDE_END, 0.2, fc="none", ec=S.INK2, lw=0.6, hatch="////"))
    ax.text(SIGNAL_PEPTIDE_END + 2, 0.3, "signal peptide", va="center", fontsize=7, color=S.INK2)
    st = d.structural
    ax.vlines(st.index[st.active_site], 0.0, 0.18, color=S.INK, lw=0.8)
    gl = st.index[st.glyco_site]
    ax.scatter(gl, np.full(len(gl), 0.09), marker="v", s=16, color=S.INK, zorder=3)
    ax.text(L + 3, 0.725, "exons", va="center", fontsize=7, color=S.INK2)
    ax.text(L + 3, 0.09, "active site (|) · N-glycan (▼)", va="center", fontsize=7, color=S.INK2)
    ax.set_ylim(-0.05, 1.0)
    ax.axis("off")


def _variant_density(ax, d):
    L = d.length
    sp = d.master[d.master.S_path].assign(pos=lambda x: pos_of(x.protein))
    cnt = {k: sp[sp.amenability_status.fillna("NA") == k].pos.value_counts().reindex(range(1, L + 1), fill_value=0)
           for k in ("amenable", "non_amenable", "NA")}
    x = np.arange(1, L + 1)
    bottom = np.zeros(L)
    for k, c, lab in (("amenable", S.BLUE, "amenable"),
                      ("non_amenable", S.ORANGE, "non-amenable"),
                      ("NA", S.GREY, "not in Amicus table")):
        ax.bar(x, cnt[k], bottom=bottom, width=0.9, color=c, edgecolor=S.SURF, linewidth=0.2, label=lab)
        bottom += cnt[k].to_numpy()
    reg = d.table("stats_E_regions")
    reg = reg[reg["sample"] == "P/LP missense (S_path)"]
    for r in reg.itertuples():
        ax.axvspan(r.start - 0.5, r.end + 0.5, color="#000000", alpha=0.05, lw=0)
        ax.text((r.start + r.end) / 2, bottom.max() * 1.03, f"{int(r.n_observed)} vs {r.n_expected:.0f} exp.",
                ha="center", va="bottom", fontsize=7, color=S.INK2)
    ax.set_ylabel("ClinVar P/LP missense\nvariants per residue")
    ax.set_ylim(0, bottom.max() * 1.32)
    S.grid_y(ax)
    ax.tick_params(labelbottom=False)
    ax.legend(loc="upper center", ncol=3, fontsize=7, handlelength=1, bbox_to_anchor=(0.5, 1.0))


def _conservation(ax, d):
    L = d.length
    cv = d.conservation
    cons = cv.conservation.where(cv.valid).reindex(range(1, L + 1))
    ax.plot(np.arange(1, L + 1), cons, color=S.INK2, lw=0.8)
    ax.set_ylim(0, 1.25)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("conservation", fontsize=7)
    ax.yaxis.set_label_coords(-0.075, 0.5)
    S.grid_y(ax)
    ax.tick_params(labelbottom=False)
    ax.text(L + 3, 0.5, "gaps: coverage\n< threshold", va="center", fontsize=6, color=S.INK2)


def _rsa(ax, d):
    L = d.length
    rsa = d.structural.rsa.reindex(range(1, L + 1))
    ax.plot(np.arange(1, L + 1), rsa, color=S.INK2, lw=0.8)
    ax.axhline(RSA_BURIED, color=S.INK2, lw=0.6, ls=(0, (3, 3)))
    ax.text(L + 3, RSA_BURIED, f"buried\n(RSA {RSA_BURIED})", va="center", fontsize=7, color=S.INK2)
    ax.set_ylabel("RSA (dimer)")
    ax.yaxis.set_label_coords(-0.075, 0.5)
    ax.set_xlabel("residue (UniProt P06280 / HGVS numbering)")
    ax.set_ylim(0, 1.25)
    S.grid_y(ax)
    ax.axvspan(0.5, SIGNAL_PEPTIDE_END + 0.5, color="#000000", alpha=0.05, lw=0)
    ax.text(SIGNAL_PEPTIDE_END / 2, 0.98, "signal peptide,\nno structure",
            ha="center", va="center", fontsize=6, color=S.INK2)


def build(d):
    fig = plt.figure(figsize=FIGSIZE)
    gs = fig.add_gridspec(4, 1, height_ratios=HEIGHT_RATIOS, hspace=0.16)
    a0, a1, a2, a3 = (fig.add_subplot(gs[i]) for i in range(4))
    for a in (a1, a2, a3):
        a.sharex(a0)
    _exon_track(a0, d)
    _variant_density(a1, d)
    _conservation(a2, d)
    _rsa(a3, d)
    a1.yaxis.set_label_coords(-0.075, 0.5)
    for ax, letter in ((a0, "a"), (a1, "b"), (a2, "c"), (a3, "d")):
        S.panel(ax, letter)
    return fig
