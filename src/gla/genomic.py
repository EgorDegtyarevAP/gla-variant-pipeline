"""Genomic (GRCh38, chrX, minus strand) <-> NM_000169 c. mapping.

No external reference is reachable offline, so the mapping is *derived from gnomAD*:
the canonical-transcript rows give (genomic position, c. name) pairs for exonic, intronic
and UTR variants. Exon offsets are inferred from them and validated by reproducing every
gnomAD SNV name (see tests/test_genomic.py and the audit report).
"""
import re
from collections import Counter

COMP = str.maketrans("ACGT", "TGCA")
_SNV = re.compile(r"^c\.(?P<b>[-*]?\d+)(?P<o>[+-]\d+)?(?P<ref>[ACGT])>(?P<alt>[ACGT])$")


class GlaMapper:
    """pos(c.N) = offset_k - N inside exon k (exon k spans c.starts[k]..c.ends[k])."""

    def __init__(self, exon_cdna, offsets):
        self.exon_cdna = exon_cdna      # [(c_first, c_last)] per exon (c.1 = A of ATG; UTR via sentinel)
        self.offsets = offsets          # per exon: genomic position + c_number
        self.pos_c1 = None

    @classmethod
    def from_gnomad(cls, rows, exon_cdna):
        """rows: iterable of (pos, 'c.xxx' transcript consequence) for SNVs of the canonical transcript."""
        votes = [Counter() for _ in exon_cdna]
        for pos, name in rows:
            m = _SNV.match(name)
            if not m:
                continue
            b, o = m["b"], int(m["o"]) if m["o"] else 0
            if b.startswith("*") or b.startswith("-"):
                continue
            n = int(b)
            # minus strand: c.N+k lies k nt below pos(c.N) -> pos(c.N) = pos + k ; c.M-k -> pos(c.M) = pos - k
            anchor_pos = pos + o
            for i, (a, z) in enumerate(exon_cdna):
                if a <= n <= z:
                    votes[i][anchor_pos + n] += 1
        offsets = []
        for v in votes:
            if not v:
                raise ValueError("exon without anchors")
            offsets.append(v.most_common(1)[0][0])
        mp = cls(exon_cdna, offsets)
        mp.vote_purity = [max(v.values()) / sum(v.values()) for v in votes]
        mp.pos_c1 = offsets[0] - 1
        mp.pos_cdna_end = offsets[-1] - exon_cdna[-1][1]      # genomic pos of c.1290 (last CDS nt)
        return mp

    def exon_bounds(self):
        """Genomic [hi, lo] for each exon using c. ranges (UTR extents given by observed variants)."""
        out = []
        for (a, z), off in zip(self.exon_cdna, self.offsets):
            out.append((off - a, off - z))
        return out

    def g_to_c(self, pos):
        """Genomic position -> c. coordinate string (no edit), nearest-exon HGVS rule."""
        pos = int(pos)
        c1, cend = self.pos_c1, self.pos_cdna_end
        if pos > c1:                       # 5' of ATG (minus strand: higher coordinate): 5'UTR or upstream flank
            return f"c.-{pos - c1}"
        if pos < cend:                     # downstream of last CDS nt
            return f"c.*{cend - pos}"
        for i, ((a, z), off) in enumerate(zip(self.exon_cdna, self.offsets)):
            hi, lo = off - a, off - z
            if lo <= pos <= hi:
                return f"c.{off - pos}"
        # intronic: between exon i (lower c.) and i+1
        for i in range(len(self.exon_cdna) - 1):
            (a, z), off = self.exon_cdna[i], self.offsets[i]
            (a2, _), off2 = self.exon_cdna[i + 1], self.offsets[i + 1]
            donor_pos = off - z            # genomic pos of last nt of exon i
            acc_pos = off2 - a2            # genomic pos of first nt of exon i+1
            if acc_pos < pos < donor_pos:
                d_don, d_acc = donor_pos - pos, pos - acc_pos
                return f"c.{z}+{d_don}" if d_don <= d_acc else f"c.{a2}-{d_acc}"
        raise ValueError(f"position {pos} cannot be mapped")

    def snv_to_c(self, pos, ref, alt):
        base = self.g_to_c(pos)
        return f"{base}{ref.translate(COMP)}>{alt.translate(COMP)}"
