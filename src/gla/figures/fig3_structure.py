"""Figure: structural context of amenable and non-amenable substitutions.

Panels: (a) ECDF of relative solvent accessibility, (b) ECDF of distance to the
bound ligand, (c) share of non-amenable variants by structural category with
Wilson 95% confidence intervals drawn as capped whiskers.
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from gla.figures import style as S
from gla.figures.data import pos_of

NAME = "fig3_structure_vs_amenability"
FIGSIZE = (7.2, 2.7)
WIDTH_RATIOS = [1, 1, 1.25]
RSA_BURIED = 0.2
CATEGORIES = [("active_site", "active\nsite"), ("dimer_interface", "dimer\ninterface"),
              ("buried_core", "buried\ncore"), ("surface", "surface")]


def _substitutions(d):
    a = d.master[d.master.S_amen].drop_duplicates("protein")[["protein", "amenability_status"]].copy()
    a["pos"] = pos_of(a.protein)
    st = d.structural
    a = a[a.pos.isin(st.index[st.structure_resolved])]
    a["rsa"] = a.pos.map(st.rsa)
    a["dist"] = a.pos.map(st.dist_ligand)
    return a


def _ecdf_panel(ax, a, stats, column, feature, xlabel):
    for k, c, lab in (("amenable", S.BLUE, "amenable"), ("non_amenable", S.ORANGE, "non-amenable")):
        S.ecdf(ax, a[a.amenability_status == k][column].to_numpy(), c, lab)
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, 1.02)
    S.grid_y(ax)
    r = stats.loc[feature]
    ax.text(0.97, 0.06, f"AUC {r.estimate:.2f} ({r.ci95_low:.2f}–{r.ci95_high:.2f})\nq = {r.q_BH:.1e}",
            transform=ax.transAxes, ha="right", fontsize=7, color=S.INK2)


def _categories(ax, cat):
    for i, (k, lab) in enumerate(CATEGORIES):
        r = cat.loc[k]
        ax.bar(i, r.frac_non_amenable, color=S.ORANGE, width=0.6)
        ax.bar(i, 1 - r.frac_non_amenable, bottom=r.frac_non_amenable, color=S.BLUE,
               width=0.6, edgecolor=S.SURF, linewidth=0.8)
        S.ci_y(ax, i, r.frac_non_amenable, r.ci95_low_nonamen, r.ci95_high_nonamen, color=S.INK)
        ax.text(i, 1.04, f"n={int(r.n_with_status)}", ha="center", fontsize=7, color=S.INK2)
    ax.set_xticks(range(len(CATEGORIES)))
    ax.set_xticklabels([l for _, l in CATEGORIES], fontsize=7.5)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("share of P/LP missense")
    S.grid_y(ax)
    ax.legend(handles=[Line2D([0], [0], marker="s", ls="", color=c, label=l, markersize=6)
                       for c, l in ((S.ORANGE, "non-amenable"), (S.BLUE, "amenable"))]
                      + [Line2D([0], [0], marker="|", ls="", color=S.INK, label="Wilson 95% CI",
                                markersize=9, markeredgewidth=1.2)],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=6.5,
              columnspacing=0.8, handletextpad=0.3)


def build(d):
    a = _substitutions(d)
    stats = d.table("stats_A_amenable_vs_nonamenable").set_index("feature")
    cat = d.table("stats_C_S_path_categories").set_index("category")
    fig, (c1, c2, c3) = plt.subplots(1, 3, figsize=FIGSIZE, gridspec_kw={"width_ratios": WIDTH_RATIOS})
    _ecdf_panel(c1, a, stats, "rsa", "rsa", "RSA of the substituted residue")
    _ecdf_panel(c2, a, stats, "dist", "dist_ligand", "distance to bound ligand (Å)")
    c1.set_ylabel("cumulative fraction")
    c1.legend(loc="center right", fontsize=7)
    c1.axvline(RSA_BURIED, color=S.INK2, lw=0.5, ls=(0, (3, 3)))
    _categories(c3, cat)
    fig.tight_layout(w_pad=1.4, rect=(0, 0, 1, 0.93))
    S.panel_fig(fig, [(c1, "a"), (c2, "b"), (c3, "c")])
    return fig
