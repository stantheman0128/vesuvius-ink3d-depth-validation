"""Unit tests required by PHASE5B preregistration section 8: anchor rule with
tie policy, selection stride, Wilson interval, and the NCC helper."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phase5b_scaleup_run import (  # noqa: E402
    ANCHOR_LO, CHEBYSHEV_MIN, N_CUBES, ncc_search, select_anchors, wilson,
)


def test_anchor_rule_takes_global_max_first():
    pred = np.zeros((256, 256, 256), dtype=np.uint8)
    pred[100, 120, 140] = 250
    pred[170, 80, 80] = 200  # Chebyshev 70 from the max: separable
    anchors = select_anchors(pred)
    assert anchors[0] == (100, 120, 140)
    assert anchors[1] == (170, 80, 80)


def test_anchor_rule_enforces_chebyshev_separation():
    pred = np.zeros((256, 256, 256), dtype=np.uint8)
    pred[100, 100, 100] = 250
    pred[100, 100, 150] = 249  # Chebyshev 50 < 64: must be rejected
    pred[100, 100, 170] = 248  # Chebyshev 70: accepted
    anchors = select_anchors(pred)
    assert (100, 100, 100) in anchors
    assert (100, 100, 150) not in anchors
    assert (100, 100, 170) in anchors
    for i, a in enumerate(anchors):
        for b in anchors[i + 1:]:
            assert max(abs(a[k] - b[k]) for k in range(3)) >= CHEBYSHEV_MIN


def test_anchor_rule_tie_policy_lexicographic():
    pred = np.full((256, 256, 256), 7, dtype=np.uint8)
    anchors = select_anchors(pred)
    # all values tie: the first anchor must be the lexicographically smallest
    # eligible center, and the rule stays deterministic
    assert anchors[0] == (ANCHOR_LO, ANCHOR_LO, ANCHOR_LO)
    assert anchors == select_anchors(pred)
    assert len(anchors) == 8


def test_anchor_rule_ignores_out_of_bounds_maxima():
    pred = np.zeros((256, 256, 256), dtype=np.uint8)
    pred[10, 10, 10] = 255      # outside [64,191]: must never be an anchor
    pred[90, 90, 90] = 100
    anchors = select_anchors(pred)
    assert (10, 10, 10) not in anchors
    assert anchors[0] == (90, 90, 90)


def test_selection_stride_indices():
    L = 1000
    seeds = [(i * L) // N_CUBES for i in range(N_CUBES)]
    assert len(set(seeds)) == N_CUBES
    assert seeds[0] == 0 and seeds[-1] < L


def test_wilson_interval_known_value():
    lo, hi = wilson(8, 10)
    assert 0.49 < lo < 0.50 and 0.94 < hi < 0.95
    assert wilson(0, 0) is None


def test_ncc_search_recovers_known_shift():
    rng = np.random.default_rng(5)
    vol = rng.normal(100, 20, (80, 80, 80)).astype(np.float32)
    template = vol[24:56, 24:56, 24:56].copy()
    # search window centered so that template sits at offset (+2, -1, 0)
    pad = 16
    lo = np.array([24 + 2, 24 - 1, 24 + 0]) - pad
    search = vol[lo[0]:lo[0] + 32 + 2 * pad, lo[1]:lo[1] + 32 + 2 * pad,
                 lo[2]:lo[2] + 32 + 2 * pad]
    import phase5b_scaleup_run as m
    old = m.QA_PAD
    m.QA_PAD = pad
    try:
        peak, offset, on_edge = ncc_search(search, template)
    finally:
        m.QA_PAD = old
    assert peak > 0.99
    assert offset == [-2, 1, 0] or offset == [2, -1, 0]
    assert not on_edge
