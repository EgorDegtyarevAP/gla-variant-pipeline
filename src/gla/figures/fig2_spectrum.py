"""Figure: variant spectrum and the migalastat coverage ceiling.

Panels: (a) ClinVar P/LP variants by class, missense split by amenability
status; (b) the four proportions that define the coverage ceiling. These are
descriptive proportions of a fixed table, so no confidence intervals are drawn
(see Methods, Section 2.7).
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


from gla.figures import style as S

NAME = "fig2_spectrum_and_coverage"
FIGSIZE = (7.2, 2.9)
WIDTH_RATIOS = [1.15, 1]

CLASS_ORDER = ["missense", "frameshift", "nonsense", "splice_site",
               "inframe_indel", "intronic", "other"]
CLASS_LABELS = {"missense": "missense", "frameshift": "frameshift", "nonsense": "nonsense",
                "splice_site": "canonical splice site", "inframe_indel": "in-frame indel",
                "intronic": "other intronic", "other": "other"}
COLLAPSE_TO_OTHER = ["mnv_missense", "start_lost", "stop_lost",
                     "splice_spanning_indel", "delins_unresolved"]


def _counts(d):
    plp = d.master[d.master.clinvar_class.isin(["P", "LP", "P/LP"])]
    cls = plp.variant_class.replace({k: "other" for k in COLLAPSE_TO_OTHER}).value_counts()
    order = [c for c in CLASS_ORDER if c in cls.index]
    mis = plp[plp.variant_class == "missense"]
    status = {"amenable": int((mis.amenability_status == "amenable").sum()),
              "non_amenable": int((mis.amenability_status == "non_amenable").sum()),
              "NA": int(mis.amenability_status.isna().sum())}
    return plp, cls, order, mis, status


def _classes(ax, cls, order, status, n_total):
    for i, c in enumerate(order):
        y = len(order) - 1 - i
        if c == "missense":
            left = 0
            for k, col in (("amenable", S.BLUE), ("non_amenable", S.ORANGE), ("NA", S.GREY)):
                ax.barh(y, status[k], left=left, color=col, height=0.62, edgecolor=S.SURF, linewidth=0.8)
                left += status[k]
        else:
            ax.barh(y, cls[c], color="#c3c2bb", height=0.62)
        ax.text(cls[c] + 8, y, f"{cls[c]} ({cls[c] / n_total:.0%})", va="center", fontsize=7.5, color=S.INK)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([CLASS_LABELS[c] for c in order][::-1])
    ax.set_xlim(0, cls.max() * 1.3)
    ax.set_xlabel(f"ClinVar P/LP variants (n = {n_total})")
    S.grid_x(ax)
    ax.legend(handles=[Line2D([0], [0], marker="s", ls="", color=c, label=l, markersize=6)
                       for c, l in ((S.BLUE, "amenable"), (S.ORANGE, "non-amenable"), (S.GREY, "not in table"))],
              loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=7, title="missense:",
              title_fontsize=7, handletextpad=0.3, alignment="left")


def _ceiling(ax, mis, status, n_total):
    am, na = status["amenable"], status["NA"]
    rows = [("missense / all P/LP", len(mis), n_total),
            ("amenable / P/LP missense", am, len(mis)),
            ("amenable / all P/LP: lower bound", am, n_total),
            ("amenable + not in table / all P/LP: upper bound", am + na, n_total)]
    # Descriptive proportions of a fixed table: no confidence intervals (Methods, 2.7).
    for i, (lab, k, n) in enumerate(rows):
        y = len(rows) - 1 - i
        p = k / n
        ax.barh(y, p, height=0.12, color=S.GREY, zorder=2)
        ax.scatter([p], [y], s=26, color=S.BLUE if i else S.INK2, zorder=3, edgecolor=S.SURF, linewidth=0.8)
        ax.text(p + 0.02, y, f"{p:.0%}  ({k}/{n})", va="center", fontsize=7.5, color=S.INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=7.5)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-0.5, len(rows) - 0.5)
    ax.set_xlabel("proportion of variants")
    S.grid_x(ax)


def build(d):
    plp, cls, order, mis, status = _counts(d)
    fig, (b1, b2) = plt.subplots(1, 2, figsize=FIGSIZE, gridspec_kw={"width_ratios": WIDTH_RATIOS})
    _classes(b1, cls, order, status, len(plp))
    _ceiling(b2, mis, status, len(plp))
    fig.tight_layout(w_pad=2.0, rect=(0, 0, 1, 0.95))
    S.panel_fig(fig, [(b1, "a"), (b2, "b")])
    return fig
