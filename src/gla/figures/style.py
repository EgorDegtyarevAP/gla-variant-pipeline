"""Shared figure style: palette, rcParams, saving, panel letters, error bars.

Edit this file to change the look of every figure at once. Edit one
``figN_*.py`` module to change a single figure.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gla.common import RESULTS

FIG = RESULTS / "figures"
TABLES = RESULTS / "tables"

# ---------------------------------------------------------------- palette ---
BLUE = "#2a78d6"      # slot 1: amenable / clinically ascertained
ORANGE = "#eb6834"    # slot 2: non-amenable
AQUA = "#1baf7a"      # slot 3: population sample (gnomAD)
GREY = "#a9a8a2"      # neutral: not in the amenability table / background
INK = "#0b0b0b"       # primary text and rules
INK2 = "#52514e"      # secondary text, axes, error bars
GRID = "#e6e5e1"
SURF = "#ffffff"

# --------------------------------------------------------------- rcParams ---
RC = {
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.edgecolor": INK2, "axes.linewidth": 0.6, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.facecolor": SURF, "figure.facecolor": SURF,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
}


def apply():
    plt.rcParams.update(RC)


# ----------------------------------------------------------------- output ---
FORMATS = ("pdf", "png", "svg")


def save(fig, name, formats=FORMATS, dpi=300):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(FIG / f"{name}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return name


# ----------------------------------------------------------- panel letters ---
def panel(ax, letter):
    ax.text(-0.125, 1.0, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", ha="left", va="top", color=INK)


def panel_fig(fig, axes_letters):
    """Panel letters at the top-left of each axes' tight bbox (no collisions)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    for ax, letter in axes_letters:
        bb = ax.get_tightbbox(r)
        fig.text(bb.x0 / W, min(bb.y1 / H + 0.01, 0.995), letter, fontsize=10,
                 fontweight="bold", ha="left", va="bottom", color=INK)


# -------------------------------------------------------------- error bars ---
# One place to change how every confidence interval is drawn. Caps are the
# short perpendicular ticks at the ends of each whisker.
CI_KW = dict(color=INK2, elinewidth=1.0, capsize=2.5, capthick=0.9)


def ci_x(ax, y, point, low, high, **kw):
    """Horizontal CI whisker with caps at (low, high) centred on ``point``."""
    k = {**CI_KW, **kw}
    ax.errorbar([point], [y], xerr=[[point - low], [high - point]], fmt="none",
                zorder=2, **k)


def ci_y(ax, x, point, low, high, **kw):
    """Vertical CI whisker with caps at (low, high) centred on ``point``."""
    k = {**CI_KW, **kw}
    ax.errorbar([x], [point], yerr=[[point - low], [high - point]], fmt="none",
                zorder=4, **k)


# ------------------------------------------------------------------ shapes ---
def ecdf(ax, x, color, label, lw=1.4):
    x = np.sort(x)
    y = np.arange(1, len(x) + 1) / len(x)
    ax.step(x, y, where="post", color=color, lw=lw, label=label)


def grid_x(ax):
    ax.xaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def grid_y(ax):
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
