"""Stage 7. Figures.

The drawing code lives in ``src/gla/figures/``: one module per figure plus a
shared ``style`` module. This script only selects which figures to build.

    python scripts/11_figures.py            # all figures
    python scripts/11_figures.py 3 4        # only figures 3 and 4
    python scripts/11_figures.py 1 --png    # one figure, PNG only (fast preview)

Figures read processed tables only; no statistic is recomputed here except
ECDFs and the Wilson intervals of already-defined samples.
Output: results/figures/fig{1..5}.{pdf,png,svg}
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from gla.figures import MODULES, render
from gla.figures.style import FORMATS


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("numbers", nargs="*", type=int, choices=sorted(MODULES), metavar="N",
                   help=f"figure numbers to build ({', '.join(map(str, sorted(MODULES)))}); default: all")
    p.add_argument("--png", action="store_true", help="write PNG only instead of PDF, PNG and SVG")
    args = p.parse_args(argv)
    written = render(args.numbers or None, formats=("png",) if args.png else FORMATS)
    print("figures written:", ", ".join(written))


if __name__ == "__main__":
    main()
