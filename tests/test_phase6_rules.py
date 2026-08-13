"""Unit tests required by PHASE6 preregistration section 7."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phase6_dual_energy_run import NCC_MIN, apply_offset, to137  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_chain_direction_on_official_landmarks():
    """canonical -> 137keV chain must map 78keV moving landmarks near the
    137keV moving landmarks that share the same fixed landmark, using the
    transform files themselves (round trip through fixed space)."""
    d78 = json.load(open(ROOT / "data_crossscan_meta_2026-08-07/transform_20260411134726.json"))
    d137 = json.load(open(ROOT / "data_crossscan_meta_2026-08-07/transform_20260323153942.json"))
    M137 = np.array(d137["transformation_matrix"], float)
    A137, t137 = M137[:, :3], M137[:, 3]
    fn = to137()
    # take a 137keV moving landmark, push it to fixed with the official
    # transform, pull it back to canonical with the 78keV inverse, then the
    # chain must return (about) the original 137keV landmark.
    M78 = np.array(d78["transformation_matrix"], float)
    A78, t78 = M78[:, :3], M78[:, 3]
    m137 = np.array(d137["moving_landmarks"], float)[0]  # (x,y,z)
    fixed = A137 @ m137 + t137
    canonical_xyz = np.linalg.inv(A78) @ (fixed - t78)
    p_zyx = canonical_xyz[::-1][None]
    out_zyx = fn(p_zyx)[0]
    assert np.allclose(out_zyx[::-1], m137, atol=1e-6)


def test_apply_offset_is_input_shift():
    fn = to137()
    p = np.array([[40000.0, 12000.0, 11000.0]])
    off = np.array([2.0, -3.0, 1.0])
    got = apply_offset(fn, off)(p)
    expect = fn(p + off)
    assert np.allclose(got, expect)
    assert not np.allclose(got, fn(p))


def test_ncc_gate_threshold_constant():
    assert NCC_MIN == 0.85
