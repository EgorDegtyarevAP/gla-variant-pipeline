"""Figure: CpG transitions among observed coding single-nucleotide variants.

Six samples (three population, three clinically ascertained) against the common
expectation from the enumerated space of possible coding SNVs. Exact
(Clopper-Pearson) 95% confidence intervals are drawn as capped whiskers.
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from gla.figures import style as S

NAME = "fig4_cpg"
FIGSIZE = (7.2, 2.6)
XMAX = 0.14
LABEL_X = 0.145


def build(d):
    E = d.table("cpg_enrichment").iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for i, r in E.iterrows():
        population = r["sample"].startswith("gnomAD")
        S.ci_x(ax, i, r.frac_observed, r.ci95_low, r.ci95_high)
        ax.scatter([r.frac_observed], [i], s=30, color=S.AQUA if population else S.BLUE,
                   zorder=3, edgecolor=S.SURF, linewidth=0.8)
        ax.scatter([r.frac_expected], [i], marker="|", s=110, color=S.INK, zorder=4, linewidth=1.4)
        ax.text(LABEL_X, i,
                f"{r.n_cpg_ts_observed}/{r.n_observed}   ×{r.ratio_obs_over_exp:.1f}   p = {r.p_binomial_two_sided:.2g}",
                va="center", fontsize=7.5, color=S.INK)
    ax.set_yticks(range(len(E)))
    ax.set_yticklabels(E["sample"], fontsize=7.5)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(-0.6, len(E) - 0.4)
    ax.set_xlabel("share of CpG transitions among observed CDS SNVs (exact 95% CI)")
    S.grid_x(ax)
    ax.text(LABEL_X, len(E) - 0.35, "k/n   enrichment   p", fontsize=7, color=S.INK2)
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color=S.AQUA, label="population (gnomAD PASS)", markersize=6),
                       Line2D([0], [0], marker="o", ls="", color=S.BLUE, label="patients (ClinVar P/LP)", markersize=6),
                       Line2D([0], [0], marker="|", ls="", color=S.INK, label="expected (possible-SNV space)",
                              markersize=9, markeredgewidth=1.4)],
              loc="upper center", ncol=3, fontsize=7, bbox_to_anchor=(0.5, -0.22))
    fig.tight_layout()
    return fig
