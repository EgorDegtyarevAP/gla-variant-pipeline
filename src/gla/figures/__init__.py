"""Figure package: one module per figure, one shared style module.

Layout
------
style.py              palette, rcParams, save(), panel letters, CI whiskers
data.py               inputs loaded once (master table, features, results tables)
fig1_map.py           map of GLA and alpha-Gal A
fig2_spectrum.py      variant spectrum and the coverage ceiling
fig3_structure.py     structure versus amenability
fig4_cpg.py           CpG transitions
fig5_conservation.py  positional conservation

To change the appearance of every figure, edit ``style.py``. To change one
figure, edit its own module; nothing else imports it. Each module exposes
``NAME`` (output file stem) and ``build(inputs) -> matplotlib Figure``.
"""
from gla.figures import (fig1_map, fig2_spectrum, fig3_structure, fig4_cpg,
                         fig5_conservation, style)
from gla.figures.data import load

#: ordinal in the manuscript -> module
MODULES = {
    1: fig1_map,
    2: fig2_spectrum,
    3: fig3_structure,
    4: fig4_cpg,
    5: fig5_conservation,
}


def render(numbers=None, formats=style.FORMATS):
    """Build and save the selected figures; returns the file stems written."""
    style.apply()
    d = load()
    written = []
    for n in sorted(numbers or MODULES):
        module = MODULES[n]
        fig = module.build(d)
        written.append(style.save(fig, module.NAME, formats=formats))
    return written
