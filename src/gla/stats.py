"""Statistical helpers. Unit of analysis is stated for every test (see docs/stats_plan.md):
- variant-level tests with a position-clustered permutation null (features are constant within a position),
- position-level tests (positions are independent units),
- descriptive proportions with Wilson intervals."""
import numpy as np
import pandas as pd
from scipy.stats import rankdata, mannwhitneyu, fisher_exact, false_discovery_control
from scipy.stats.contingency import odds_ratio


def wilson(k, n, z=1.959964):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def auc(x, y):
    """P(x_pos > x_neg) + 0.5 ties, y boolean."""
    r = rankdata(x); n1 = y.sum(); n0 = len(y) - n1
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def _prop_diff(x, y):
    return x[y].mean() - x[~y].mean()


def cluster_permutation(pos_idx, y, feat_by_pos, stat, n_perm, rng):
    """Permute the position->feature assignment (keeps within-position label structure): two-sided p."""
    obs = stat(feat_by_pos[pos_idx], y)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = stat(rng.permutation(feat_by_pos)[pos_idx], y)
    c = null.mean()
    p = (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (n_perm + 1)
    return obs, p


def cluster_bootstrap_ci(pos_idx, y, feat_by_pos, stat, n_boot, rng, alpha=0.05):
    """Resample positions (clusters) with replacement; percentile CI of the statistic."""
    groups = [np.flatnonzero(pos_idx == k) for k in range(len(feat_by_pos))]
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        idx = np.concatenate([groups[k] for k in pick if len(groups[k])])
        yy = y[idx]
        if yy.all() or (~yy).all():
            continue
        vals.append(stat(feat_by_pos[pos_idx[idx]], yy))
    return tuple(np.quantile(vals, [alpha / 2, 1 - alpha / 2]))


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    out[ok] = false_discovery_control(p[ok], method="bh")
    return out


def fisher_or(a, b, c, d):
    """2x2 [[a,b],[c,d]] -> OR, CI, p (exact)."""
    t = [[a, b], [c, d]]
    o, p = fisher_exact(t)
    ci = odds_ratio(t).confidence_interval(0.95)
    return o, ci.low, ci.high, p


def window_scan(pos_obs, pos_possible_counts, width, n_perm, rng):
    """Window count of observed variants vs expectation from mutational opportunity.
    pos_obs: protein positions of observed variants (one entry per variant);
    pos_possible_counts: array len L+1 (index = position) of possible missense SNVs.
    Null: the same number of variants drawn without replacement from the possible-SNV pool.
    Returns per-window table with observed, expected, local p (upper tail) and family-wise p (max-statistic)."""
    L = len(pos_possible_counts) - 1
    pool = np.repeat(np.arange(L + 1), pos_possible_counts.astype(int))
    n = len(pos_obs)
    starts = np.arange(1, L - width + 2)
    cs = lambda counts: np.concatenate([[0], np.cumsum(counts)])
    def win(counts):
        c = cs(counts); return c[starts + width - 1] - c[starts - 1]
    obs = win(np.bincount(pos_obs, minlength=L + 1))
    exp = n * win(pos_possible_counts) / pos_possible_counts.sum()
    null_max = np.empty(n_perm); ge = np.zeros(len(starts))
    for i in range(n_perm):
        draw = rng.choice(pool, size=n, replace=False)
        w = win(np.bincount(draw, minlength=L + 1))
        null_max[i] = w.max(); ge += (w >= obs)
    return pd.DataFrame({"start": starts, "end": starts + width - 1, "observed": obs, "expected": exp,
                         "p_local": (ge + 1) / (n_perm + 1),
                         "p_fwer": [(1 + np.sum(null_max >= o)) / (n_perm + 1) for o in obs]})


def within_cluster_permutation(pos_idx, y, x, stat, n_perm, rng):
    """For features that vary WITHIN a position (e.g. BLOSUM62): permute the labels among the variants of the same position.
    Null = no association between label and feature given the position."""
    groups = [np.flatnonzero(pos_idx == k) for k in np.unique(pos_idx)]
    groups = [g for g in groups if len(g) > 1]
    obs = stat(x, y)
    null = np.empty(n_perm)
    for i in range(n_perm):
        yy = y.copy()
        for g in groups:
            yy[g] = y[rng.permutation(g)]
        null[i] = stat(x, yy)
    c = null.mean()
    return obs, (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (n_perm + 1)
