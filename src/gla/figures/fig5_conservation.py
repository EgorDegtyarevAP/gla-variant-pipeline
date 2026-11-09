"""Figure: positional conservation (level 3).

Panels: (a) ECDFs of conservation by group, (b) within-quartile shares,
(c) joint distribution of conservation and relative solvent accessibility.
"""
import numpy as np
import matplotlib.pyplot as plt

from gla.figures import style as S
from gla.figures.data import pos_of

NAME = "fig5_conservation"
FIGSIZE = (7.2, 2.8)
WIDTH_RATIOS = [1.1, 1, 1.1]


def _groups(d):
    M = d.master
    plp_pos = set(pos_of(M[M.S_path].protein))
    gn = M[(M.variant_class == "missense") & M.gnomad_pass.fillna(False).astype(bool)]
    gn_pos = set(pos_of(gn.protein))
    am = M[M.S_amen]
    am = am.assign(pos=pos_of(am.protein))
    by_status = {"amenable": set(am[am.amenability_status == "amenable"].pos),
                 "non_amenable": set(am[am.amenability_status == "non_amenable"].pos)}
    return plp_pos, gn_pos, by_status


def _ecdfs(ax, cv, gn_pos, by_status):
    for k, c, lab in (("amenable", S.BLUE, "amenable"), ("non_amenable", S.ORANGE, "non-amenable")):
        S.ecdf(ax, cv.loc[cv.index.isin(by_status[k]), "conservation"].to_numpy(), c, lab)
    S.ecdf(ax, cv.loc[cv.index.isin(gn_pos), "conservation"].to_numpy(), S.AQUA, "gnomAD missense")
    S.ecdf(ax, cv.conservation.to_numpy(), S.GREY, "all valid positions", lw=1.0)
    ax.set_xlabel("conservation")
    ax.set_ylabel("cumulative share of positions")
    S.grid_y(ax)
    ax.legend(loc="upper left", fontsize=6.5, handlelength=1.2)


def _quartiles(ax, Q):
    x = np.arange(4)
    ax.bar(x - 0.2, Q.frac_positions_with_PLP, width=0.38, color=S.BLUE, label="positions with P/LP")
    ax.bar(x + 0.2, Q.non_amenable_share_S_amen, width=0.38, color=S.ORANGE, label="non-amenable share")
    ax.set_xticks(x)
    ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4"])
    ax.set_xlabel("conservation quartile (Q1 least)")
    ax.set_ylim(0, 1.15)
    S.grid_y(ax)
    ax.legend(loc="upper left", fontsize=6.5, handlelength=1)


def _joint(ax, cv, structural, plp_pos):
    rsa = structural.reindex(cv.index).rsa
    carries = cv.index.isin(plp_pos)
    ax.scatter(rsa[~carries], cv.conservation[~carries], s=6, color=S.GREY, linewidth=0, label="no P/LP missense")
    ax.scatter(rsa[carries], cv.conservation[carries], s=7, color=S.BLUE, linewidth=0, label="P/LP missense")
    ax.set_xlabel("RSA (dimer)")
    ax.set_ylabel("conservation")
    ax.legend(loc="lower right", fontsize=6.5, handlelength=1, markerscale=1.5)


def build(d):
    cv = d.conservation[d.conservation.valid]
    plp_pos, gn_pos, by_status = _groups(d)
    fig, (d1, d2, d3) = plt.subplots(1, 3, figsize=FIGSIZE, gridspec_kw={"width_ratios": WIDTH_RATIOS})
    _ecdfs(d1, cv, gn_pos, by_status)
    _quartiles(d2, d.table("stats_F_conservation_quartiles"))
    _joint(d3, cv, d.structural, plp_pos)
    fig.tight_layout(w_pad=1.4)
    S.panel_fig(fig, [(d1, "a"), (d2, "b"), (d3, "c")])
    return fig
