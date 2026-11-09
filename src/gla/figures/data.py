"""Inputs for the figures: processed tables only, loaded once and cached.

No figure module recomputes a statistic. Everything except ECDFs and the
Wilson intervals of already-defined samples comes from results/tables/.
"""
from dataclasses import dataclass
from functools import lru_cache

import pandas as pd

from gla.common import PROCESSED
from gla import reference
from gla.figures.style import TABLES

PROTEIN_LENGTH = 429


def pos_of(s):
    """Residue number from an HGVS protein string (p.Asn215Ser -> 215)."""
    return s.str.extract(r"(\d+)")[0].astype(int)


@dataclass
class Inputs:
    ref: object                 # reference record (exon structure, UTR lengths)
    master: pd.DataFrame        # variants_master: one row per coding-level change
    structural: pd.DataFrame    # features_structural, indexed by protein_pos
    conservation: pd.DataFrame  # features_conservation, indexed by protein_pos
    length: int = PROTEIN_LENGTH

    def table(self, name):
        """A results/tables CSV by file stem."""
        return pd.read_csv(TABLES / f"{name}.csv")


@lru_cache(maxsize=1)
def load():
    return Inputs(
        ref=reference.load(),
        master=pd.read_parquet(PROCESSED / "variants_master.parquet"),
        structural=pd.read_parquet(PROCESSED / "features_structural.parquet").set_index("protein_pos"),
        conservation=pd.read_parquet(PROCESSED / "features_conservation.parquet").set_index("protein_pos"),
    )
