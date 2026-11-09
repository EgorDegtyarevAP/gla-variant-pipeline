"""Light-weight HGVS c. handling for NM_000169 (GLA).

Scope: syntactic normalisation only (no reference sequence is available offline).
Indel keys are therefore NOT 3'-shifted; ambiguous indel matches across sources are
flagged in the audit rather than silently merged. See docs/schema.md.
"""
import re

# Exon boundaries of NM_000169 (last c. position of exon n / first c. position of exon n+1).
# Verified against Amicus' own IVS -> c. conversions and gnomAD intron names (tests/test_hgvs.py).
EXON_ENDS = [194, 369, 547, 639, 801, 999]      # donor side of IVS1..IVS6
EXON_STARTS = [195, 370, 548, 640, 802, 1000]   # acceptor side of IVS1..IVS6
CDS_LEN = 1290
UTR5_LEN, UTR3_LEN = 22, 6     # NM_000169.3 (GenBank); asserted against the record in tests

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}

POS = r"(?:[-*]?\d+(?:[+-]\d+)?)"
_SNV = re.compile(rf"^c\.(?P<pos>{POS})(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
_RANGE = re.compile(rf"^c\.(?P<s>{POS})(?:_(?P<e>{POS}))?(?P<edit>.*)$")
_IVS = re.compile(r"^IVS(?P<n>\d+)(?P<sign>[+-])(?P<k>\d+)(?P<rest>.*)$")


def ivs_to_c(s):
    """'IVS4-1G>T' -> 'c.640-1G>T'. Returns input unchanged if not IVS notation."""
    m = _IVS.match(s.strip())
    if not m:
        return s
    n, k = int(m["n"]), int(m["k"])
    if not 1 <= n <= 6:
        raise ValueError(f"IVS{n} outside GLA intron range")
    anchor = EXON_ENDS[n - 1] if m["sign"] == "+" else EXON_STARTS[n - 1]
    return f"c.{anchor}{m['sign']}{k}{m['rest']}"


def norm_c(s):
    """Canonical c. string, or None if the string cannot be parsed to a single locus.

    - strips whitespace and trailing parenthetical comments ('(putative splicing site)')
    - converts IVS notation
    - drops explicit deleted/duplicated bases (delA -> del, dupCATG -> dup)
    - delinsX keeps inserted sequence only (delTGCACinsGCTCG -> delinsGCTCG)
    """
    if s is None or (isinstance(s, float) and s != s):
        return None
    s = re.sub(r"\s*\(.*?\)\s*$", "", str(s).strip())
    s = ivs_to_c(s)
    if not s.startswith("c."):
        return None
    if "?" in s or " or " in s or "[" in s or "/" in s or ";" in s:
        return None
    s = re.sub(r"del[ACGT]+ins([ACGT]+)$", r"delins\1", s)
    s = re.sub(r"(del|dup)[ACGT]+$", r"\1", s)
    if _SNV.match(s) or (_RANGE.match(s) and re.search(r"(del|dup|ins[ACGT]+|delins[ACGT]+)$", s)):
        return s
    return None


def split_alternatives(s):
    """'c.70T>C or c.70T>A' -> ['c.70T>C', 'c.70T>A'] (Amicus 'alternative' rows)."""
    return [p.strip() for p in re.split(r"\s+or\s+", str(s)) if p.strip()]


def _num(p):
    """Return (base_pos_signed, offset) for a c. position token."""
    m = re.match(r"^(?P<b>[-*]?\d+)(?P<o>[+-]\d+)?$", p)
    b = m["b"]
    return b, int(m["o"]) if m["o"] else 0


def region(c):
    """coarse region of a canonical c. string: utr5 / cds / utr3 / intron."""
    m = _RANGE.match(c)
    b, off = _num(m["s"])
    if off:
        return "intron"
    if b.startswith("-"):
        return "utr5"
    if b.startswith("*"):
        return "utr3"
    return "cds"


def splice_distance(c):
    """Absolute intronic offset of the first position (None if exonic)."""
    m = _RANGE.match(c)
    _, off = _num(m["s"])
    return abs(off) if off else None


def edit_length(c):
    """Number of nucleotides changed by a del/dup/ins; None if not a plain del/dup/ins."""
    m = _RANGE.match(c)
    e = m["edit"]
    if e in ("del", "dup"):
        s_b, s_o = _num(m["s"])
        if m["e"] is None:
            return 1
        e_b, e_o = _num(m["e"])
        if s_o or e_o:
            # intronic range; same base anchor -> offsets difference
            if s_b == e_b:
                return abs(e_o - s_o) + 1
            return None
        return abs(int(e_b.lstrip("*")) - int(s_b.lstrip("*"))) + 1
    m2 = re.match(r"^ins([ACGT]+)$", e)
    if m2:
        return len(m2[1])
    return None


def delins_net(c):
    """Net nucleotide change of a c.XdelinsY / c.X_YdelinsZ (inserted - deleted); None if not computable."""
    m = _RANGE.match(c)
    ins = re.match(r"^delins([ACGT]+)$", m["edit"])
    if not ins:
        return None
    s_b, s_o = _num(m["s"])
    if m["e"] is None:
        deleted = 1
    else:
        e_b, e_o = _num(m["e"])
        if s_o or e_o:
            return None
        deleted = abs(int(e_b.lstrip("*")) - int(s_b.lstrip("*"))) + 1
    return len(ins[1]) - deleted


def parse_protein(p):
    """Protein change -> canonical one-letter string.
    Accepts 3-letter HGVS ('p.Gly35Arg'), 1-letter ('G35R', 'W340X'), 'p.Asn336=' etc.
    Returns '' for unknown/uninformative values."""
    if p is None or (isinstance(p, float) and p != p):
        return ""
    p = str(p).strip()
    if p in ("", "-", "UNKNOWN") or p.startswith("UNKNOWN"):
        return ""
    p = p.replace("p.", "").replace("(", "").replace(")", "")
    m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)$", p)
    if m and m[1] in AA3:
        alt = AA3.get(m[3], "=") if m[3] != "=" else "="
        return f"{AA3[m[1]]}{m[2]}{AA3[m[1]] if alt == '=' else alt}"
    m = re.match(r"^([A-Z])(\d+)([A-Z*=X])$", p)
    if m:
        ref = "*" if m[1] == "X" else m[1]
        alt = "*" if m[3] == "X" else m[3]
        return f"{ref}{m[2]}{ref if alt == '=' else alt}"
    m = re.match(r"^([A-Z][a-z]{2})(\d+)(.*)$", p)   # fs, del, dup, ins, ext ...
    if m and m[1] in AA3:
        return f"{AA3[m[1]]}{m[2]}{m[3].replace('Ter', '*')}"
    return p


def protein_kind(p1):
    """Variant effect from canonical one-letter protein string."""
    if not p1:
        return None
    m = re.match(r"^([A-Z*])(\d+)(.*)$", p1)
    if not m:
        return None
    ref, pos, rest = m[1], int(m[2]), m[3]
    if ref == "*":
        return "stop_lost"
    if "fs" in rest:
        return "frameshift"
    if any(t in rest for t in ("del", "dup", "ins")):
        return "inframe_indel"
    if rest == ref:
        return "synonymous"
    if rest == "*":
        return "nonsense"
    if pos == 1 and ref == "M":
        return "start_lost"
    if len(rest) == 1 and rest.isalpha():
        return "missense"
    if rest.startswith("ext"):
        return "stop_lost"
    return "other"


def classify(c, p1):
    """Variant class for the master table (single rule set, applied to every row).

    Non-coding regions first (utr5/utr3/splice_site/intronic); then
    - SNVs: protein effect from the best available protein string (else 'cds_snv_unannotated',
      resolved at stage 3 from the CDS);
    - del/dup/ins: nucleotide length rule (multiple of 3 -> inframe_indel, else frameshift),
      NOT the protein string, because ClinVar/Amicus protein notations for indels are inconsistent;
    - ranges that cross an exon/intron boundary -> splice_spanning_indel; delins -> delins_unresolved.
    """
    reg = region(c)
    if reg in ("utr5", "utr3"):
        n = int(re.match(r"^c\.[-*](\d+)", c)[1])
        if n > (UTR5_LEN if reg == "utr5" else UTR3_LEN):
            return "upstream_flank" if reg == "utr5" else "downstream_flank"
        return reg
    if reg == "intron":
        d = splice_distance(c)
        return "splice_site" if d is not None and d <= 2 else "intronic"
    m = _RANGE.match(c)
    if _SNV.match(c):
        return protein_kind(p1) or "cds_snv_unannotated"
    if m["e"] is not None and re.search(r"[+-]\d+$", m["e"]):
        return "splice_spanning_indel"
    if m["edit"].startswith("delins"):
        net = delins_net(c)
        if net is None:
            return "delins_unresolved"
        if net % 3:
            return "frameshift"
        return "mnv_missense" if protein_kind(p1) == "missense" else "inframe_indel"
    n = edit_length(c)
    if n is None:
        return "other"
    if p1 and "ext" in p1:
        return "stop_lost"
    return "inframe_indel" if n % 3 == 0 else "frameshift"


def proteins_compatible(a, b):
    """True if two canonical protein strings describe the same event at the precision the notations allow.
    Substitutions must match exactly; fs/indel strings only need the same first position (+-3),
    because tools shift repeats differently (3' vs 5')."""
    if not a or not b:
        return True
    ka, kb = protein_kind(a), protein_kind(b)
    sub = ("missense", "synonymous", "nonsense", "start_lost")
    if ka in sub and kb in sub:
        return a == b
    pa, pb = re.search(r"\d+", a), re.search(r"\d+", b)
    return bool(pa and pb and abs(int(pa[0]) - int(pb[0])) <= 3)
