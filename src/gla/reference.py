"""NM_000169.3 reference (GenBank) -> sequence checks, codon-level SNV consequences, CpG context."""
import re
from functools import lru_cache
from Bio import SeqIO
from Bio.Seq import Seq

from .common import RAW
from .hgvs import AA3

GB = "NM_000169.3.gb"
ONE = {v: k for k, v in AA3.items()}


class Reference:
    def __init__(self, path=None):
        rec = SeqIO.read(path or RAW / GB, "genbank")
        self.id = rec.id
        self.mrna = str(rec.seq).upper()
        cds = next(f for f in rec.features if f.type == "CDS")
        self.cds_start = int(cds.location.start)            # 0-based mRNA index of c.1
        self.cds = self.mrna[self.cds_start:int(cds.location.end)]
        self.utr5_len = self.cds_start
        self.utr3_len = len(self.mrna) - int(cds.location.end)
        self.protein = str(Seq(self.cds).translate(to_stop=False))      # includes terminal '*'
        self.protein_gb = cds.qualifiers["translation"][0]
        self.exons_mrna = [(int(f.location.start) + 1, int(f.location.end)) for f in rec.features if f.type == "exon"]
        self.sig_peptide = next((int(f.location.start) + 1, int(f.location.end)) for f in rec.features if f.type == "sig_peptide")
        self.glyco_codons = sorted({(int(f.location.start) - self.cds_start) // 3 + 1 for f in rec.features
                                    if f.type == "misc_feature" and "glycosylation" in f.qualifiers.get("note", [""])[0]})

    # ---- invariants
    def self_check(self):
        assert self.protein.rstrip("*") == self.protein_gb, "CDS translation != GenBank translation"
        assert len(self.cds) == 1290 and self.cds[:3] == "ATG"
        return {"mRNA_len": len(self.mrna), "cds_len": len(self.cds), "utr5": self.utr5_len, "utr3": self.utr3_len,
                "protein_len": len(self.protein_gb), "exon_ends_c": [e - self.utr5_len for _, e in self.exons_mrna[:-1]],
                "glyco_codons": self.glyco_codons}

    def exon_junction_c(self):
        """c. numbers (last base of each exon except the final one) - used to flag CpGs spanning a junction."""
        return {e - self.utr5_len for _, e in self.exons_mrna[:-1]}

    # ---- position helpers (exonic positions only)
    def base_at(self, c_token):
        """Reference base at an exonic c. token ('123', '-5', '*2'); None if outside the mRNA."""
        if c_token.startswith("-"):
            i = self.cds_start - int(c_token[1:])
        elif c_token.startswith("*"):
            i = self.cds_start + len(self.cds) + int(c_token[1:]) - 1
        else:
            i = self.cds_start + int(c_token) - 1
        return self.mrna[i] if 0 <= i < len(self.mrna) else None

    def snv_effect(self, n, ref, alt):
        """Codon-level effect of SNV at c.n (1..1290). Returns dict(codon, aa_ref, aa_alt, protein, ok_ref)."""
        i = n - 1
        cod_i = i // 3
        codon = self.cds[cod_i * 3: cod_i * 3 + 3]
        ok = self.cds[i] == ref
        mut = list(codon); mut[i % 3] = alt
        a0 = str(Seq(codon).translate()); a1 = str(Seq("".join(mut)).translate())
        return {"codon_no": cod_i + 1, "codon_pos": i % 3 + 1, "ref_ok": ok, "aa_ref": a0, "aa_alt": a1,
                "protein": f"{a0}{cod_i + 1}{a0 if a0 == a1 else a1}"}


@lru_cache(maxsize=1)
def load():
    ref = Reference()
    ref.self_check()
    return ref
