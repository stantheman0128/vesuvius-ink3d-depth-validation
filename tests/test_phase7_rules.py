"""Unit tests required by PHASE7 preregistration section 9."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phase7_w00_depth_labels_run import (  # noqa: E402
    BAND_HALF, CLAMP, L_Z, band_block, correction_and_status, node_normal,
    signed_surface_offset, stride_subsample, valid_mask,
)


def test_valid_mask_rejects_minus_one_and_nan():
    x = np.array([[1.0, -1.0], [np.nan, 5.0]])
    y = np.ones((2, 2))
    z = np.ones((2, 2))
    m = valid_mask(x, y, z)
    assert m.tolist() == [[True, False], [False, True]]


def test_node_normal_on_synthetic_plane():
    # plane: x = col, y = row, z = 0.5*col  (tilted around the y axis)
    rr, cc = np.mgrid[0:9, 0:9].astype(float)
    xm, ym, zm = cc, rr, 0.5 * cc
    n = node_normal(xm, ym, zm, 4, 4)
    expect = np.array([1.0, 0.0, -0.5])  # (z,y,x), up to sign
    expect = expect / np.linalg.norm(expect)
    assert abs(abs(float(np.dot(n, expect))) - 1.0) < 1e-9


def test_signed_surface_offset_sign_and_position():
    s = np.arange(-24, 25, dtype=float)
    prof = 60 + 120 / (1 + np.exp(-(s - 6.0) / 1.5))  # step at s=+6
    got = signed_surface_offset(prof)
    assert got is not None and abs(got - 6.0) <= 0.5
    prof2 = 60 + 120 / (1 + np.exp(-(s + 6.0) / 1.5))  # step at s=-6
    got2 = signed_surface_offset(prof2)
    assert got2 is not None and abs(got2 + 6.0) <= 0.5


def test_stride_subsample_cap_and_determinism():
    items = [(i, 0) for i in range(1000)]
    out = stride_subsample(items, 800)
    assert len(out) == 800
    assert out == stride_subsample(items, 800)
    assert stride_subsample(items[:10], 800) == items[:10]


def test_correction_field_clamp_and_unvalidated_radius():
    nodes = [(10, 10), (10, 40)]
    vals = [20.0, -20.0]  # beyond clamp: must clamp to +/-CLAMP
    field, status = correction_and_status(nodes, vals, (60, 60), radius=24.0)
    assert field[10, 10] == CLAMP and field[10, 40] == -CLAMP
    assert status[10, 10] == 2 and status[10, 40] == 2
    assert status[10, 25] == 1            # between nodes: interpolated
    assert status[59, 59] == 0            # far away: unvalidated
    assert field[59, 59] == 0.0


def test_band_block_edges_and_ink_gate():
    ink = np.zeros((2, 2), np.uint8)
    ink[0, 0] = 1
    centers = np.array([[0.0, 0.0], [0.0, 0.0]]) + np.array([[1, 0], [0, 0]])
    # pixel (0,0): center 1 -> band clamped to [0, 1+BAND_HALF]
    blk = band_block(ink, centers)
    assert blk.shape == (L_Z, 2, 2)
    assert blk[0, 0, 0] == 255 and blk[1 + BAND_HALF, 0, 0] == 255
    assert blk[2 + BAND_HALF, 0, 0] == 0
    assert blk[:, 0, 1].sum() == 0        # no ink: no label
    hi = np.full((1, 1), float(L_Z + 10))
    blk2 = band_block(np.ones((1, 1), np.uint8), hi)
    assert blk2[L_Z - 1, 0, 0] == 255     # clamped at the top edge
