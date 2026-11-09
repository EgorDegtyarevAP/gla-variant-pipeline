import numpy as np
from gla import stats as st

def test_auc_extremes_and_ties():
    y = np.array([True, True, False, False])
    assert st.auc(np.array([3, 4, 1, 2.]), y) == 1.0
    assert st.auc(np.array([1, 2, 3, 4.]), y) == 0.0
    assert st.auc(np.ones(4), y) == 0.5

def test_wilson_known_values():
    lo, hi = st.wilson(0, 10); assert lo == 0 or lo < 1e-12
    assert abs(hi - 0.2775) < 2e-3
    lo, hi = st.wilson(50, 100); assert abs(lo - 0.4038) < 2e-3 and abs(hi - 0.5962) < 2e-3

def test_cluster_permutation_detects_and_ignores():
    rng = np.random.default_rng(1)
    n_pos = 60; reps = 8
    pos_idx = np.repeat(np.arange(n_pos), reps)
    feat = rng.normal(size=n_pos)
    y_signal = (feat[pos_idx] + rng.normal(scale=.5, size=len(pos_idx))) > 0
    _, p1 = st.cluster_permutation(pos_idx, y_signal, feat, st.auc, 500, np.random.default_rng(2))
    y_noise = rng.random(len(pos_idx)) > .5
    _, p2 = st.cluster_permutation(pos_idx, y_noise, feat, st.auc, 500, np.random.default_rng(3))
    assert p1 < 0.01 and p2 > 0.05

def test_window_scan_null_is_not_significant_and_hotspot_is():
    L = 200; opp = np.r_[0, np.full(L, 3.)]
    pool = np.repeat(np.arange(L + 1), opp.astype(int))
    rng = np.random.default_rng(0)
    null_obs = rng.choice(pool, 80, replace=False)
    w0 = st.window_scan(null_obs, opp, 20, 300, np.random.default_rng(1))
    assert w0.p_fwer.min() > 0.01
    hot = np.r_[rng.choice(pool, 40, replace=False), np.full(25, 100)]
    w1 = st.window_scan(hot, opp, 20, 300, np.random.default_rng(1))
    assert w1.p_fwer.min() < 0.01 and 81 <= w1.loc[w1.p_fwer.idxmin(), "start"] + 20 <= 121
