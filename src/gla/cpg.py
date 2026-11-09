"""Level 1 - mutational context (CpG) in the NM_000169.3 CDS.

Definitions (fixed here, used everywhere):
- CpG site: base is the C of a 'CG' or the G of a 'CG' on the coding strand (strand-symmetric, so the
  template strand gives the same sites). A dinucleotide spanning an exon-exon junction is NOT counted
  (the genomic neighbour is intronic and unknown): such positions are flagged `junction_ambiguous`.
- CpG transition: C>T at the C of CpG, or G>A at the G of CpG (the deamination-driven hypermutable class).
- Possible-SNV space = all 3 x 1290 substitutions of the CDS; the expected share of a class is its share
  among possible SNVs with the same consequence set as the observed sample.
"""
import pandas as pd
from scipy.stats import binomtest

from .hgvs import protein_kind

BASES = "ACGT"


def possible_snvs(ref):
    junction = ref.exon_junction_c()                    # last c. base of exons 1..6
    cds = ref.cds
    rows = []
    for n in range(1, len(cds) + 1):
        b = cds[n - 1]
        nxt = cds[n] if n < len(cds) else ""
        prv = cds[n - 2] if n > 1 else ""
        span_next = n in junction          # pair (n, n+1) crosses a junction
        span_prev = (n - 1) in junction    # pair (n-1, n) crosses a junction
        cpg_c = b == "C" and nxt == "G" and not span_next
        cpg_g = b == "G" and prv == "C" and not span_prev
        amb = (b == "C" and nxt == "G" and span_next) or (b == "G" and prv == "C" and span_prev)
        for alt in BASES:
            if alt == b:
                continue
            e = ref.snv_effect(n, b, alt)
            rows.append({
                "c_pos": n, "ref": b, "alt": alt, "hgvs_c": f"c.{n}{b}>{alt}",
                "consequence": protein_kind(e["protein"]) or "other", "protein": e["protein"],
                "cpg_site": cpg_c or cpg_g, "junction_ambiguous": amb,
                "cpg_transition": (cpg_c and alt == "T") or (cpg_g and alt == "A"),
                "transition": {b, alt} in ({"A", "G"}, {"C", "T"}),
            })
    return pd.DataFrame(rows)


def enrichment(obs_keys, possible, consequences, label):
    """Share of CpG transitions among observed SNVs vs among possible SNVs of the same consequence set."""
    space = possible[possible.consequence.isin(consequences) & ~possible.junction_ambiguous]
    obs = space[space.hgvs_c.isin(set(obs_keys))]
    n_obs, k_obs = len(obs), int(obs.cpg_transition.sum())
    n_exp, k_exp = len(space), int(space.cpg_transition.sum())
    p_exp = k_exp / n_exp
    t = binomtest(k_obs, n_obs, p_exp) if n_obs else None
    ci = t.proportion_ci(method="exact") if n_obs else (float("nan"),) * 2
    # Expected COUNT must be derived from the same p_exp as the test. Reporting it from a
    # global CpG proportion instead of the class-specific one is the error that reached the
    # manuscript text (expected 2.2/1.0/5.8/1.1 instead of 2.0/1.2/5.2/2.4); the assertion below
    # makes any future divergence between the printed expectation and the tested null fail loudly.
    n_cpg_ts_expected = n_obs * p_exp
    ratio = (k_obs / n_obs) / p_exp if n_obs else float("nan")
    if n_obs:
        assert abs(ratio - k_obs / n_cpg_ts_expected) < 1e-9, (
            f"{label}: printed expectation {n_cpg_ts_expected:.3f} is inconsistent with the tested null {p_exp:.6f}")
    out = {"sample": label, "consequences": "+".join(sorted(consequences)), "n_observed": n_obs, "n_cpg_ts_observed": k_obs,
           "frac_observed": k_obs / n_obs if n_obs else float("nan"), "ci95_low": ci[0], "ci95_high": ci[1],
           "n_possible": n_exp, "n_cpg_ts_possible": k_exp, "frac_expected": p_exp,
           "n_cpg_ts_expected": n_cpg_ts_expected,
           "ratio_obs_over_exp": ratio,
           "p_binomial_two_sided": t.pvalue if t else float("nan")}
    if n_obs:
        out["mde_fold_80pct_power"] = mde_fold(n_obs, p_exp)
        out["power_at_observed_fold"] = power_fold(n_obs, p_exp, ratio)
    return out


def _crit(n, p0, alpha=0.05):
    """Smallest upper-tail count whose two-sided exact binomial p falls below alpha."""
    from math import floor
    for k in range(floor(n * p0), n + 1):
        if binomtest(k, n, p0).pvalue < alpha:
            return k
    return None


def power_fold(n, p0, fold, alpha=0.05):
    """Exact power of the two-sided binomial test against a true rate of fold * p0."""
    from scipy.stats import binom
    k = _crit(n, p0, alpha)
    if k is None:
        return float("nan")
    return float(1 - binom.cdf(k - 1, n, min(p0 * fold, 1.0)))


def mde_fold(n, p0, target=0.80, alpha=0.05, hi=40.0, step=0.01):
    """Minimum enrichment detectable with `target` power; nan if unreachable below `hi`."""
    f = 1.0
    while f < hi:
        if power_fold(n, p0, f, alpha) >= target:
            return round(f, 2)
        f += step
    return float("nan")
