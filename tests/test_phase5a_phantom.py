"""Phantom tests required by PHASE5A preregistration section 11.

(a) a known synthetic planar surface is recovered within 0.5 voxel;
(b) a known synthetic FWHM is recovered within 1 voxel;
(c) a profile that never falls below half height on one side is censored.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phase5a_depth_attribution_run import (  # noqa: E402
    FWHM_CENSOR, WIN, analyze_candidate, fwhm_from_profile,
    surface_D_from_profile,
)


def make_plane_phantom(normal, d0, steepness=1.5):
    """69^3 volume: smooth step across a plane at signed distance d0 from center."""
    n = 2 * WIN + 1
    zz, yy, xx = np.meshgrid(*(np.arange(n) - WIN,) * 3, indexing="ij")
    s = zz * normal[0] + yy * normal[1] + xx * normal[2]
    return (60 + 120 / (1 + np.exp(-(s - d0) / steepness))).astype(np.float32)


def test_planar_surface_recovered_within_half_voxel():
    normal = np.array([0.3, 0.8, 0.52])
    normal = normal / np.linalg.norm(normal)
    d0 = 5.0
    vol = make_plane_phantom(normal, d0)
    pred = np.zeros((256, 256, 256), dtype=np.uint8)
    rec = analyze_candidate(vol, pred, (128, 128, 128))
    assert rec["anisotropy_ok"]
    got = np.array(rec["normal_zyx"])
    assert abs(float(np.dot(got, normal))) > 0.99
    assert rec["D_voxels"] is not None
    assert abs(rec["D_voxels"] - d0) <= 0.5


def test_planar_surface_at_center():
    normal = np.array([0.0, 0.0, 1.0])
    vol = make_plane_phantom(normal, 0.0)
    pred = np.zeros((256, 256, 256), dtype=np.uint8)
    rec = analyze_candidate(vol, pred, (128, 128, 128))
    assert rec["D_voxels"] is not None
    assert rec["D_voxels"] <= 0.5


def test_fwhm_recovered_within_one_voxel():
    s = np.arange(-24, 25, dtype=float)
    sigma = 3.0
    prof = 10 + 200 * np.exp(-(s ** 2) / (2 * sigma ** 2))
    expected = 2.354820045 * sigma
    fwhm, censored = fwhm_from_profile(prof)
    assert not censored
    assert abs(fwhm - expected) <= 1.0


def test_fwhm_censored_when_no_crossing():
    s = np.arange(-24, 25, dtype=float)
    prof = 100 + 0.5 * s  # monotone: no crossing on the peak's right side
    fwhm, censored = fwhm_from_profile(prof)
    assert censored
    assert fwhm == FWHM_CENSOR


def test_flat_profile_yields_no_surface():
    prof = np.full(49, 77.0)
    assert surface_D_from_profile(prof) is None


def test_gate_tables_match_preregistration():
    from phase5a_depth_attribution_run import G2_PRIMARY, G2_SECONDARY
    assert G2_PRIMARY == {4: 4, 5: 4, 6: 5, 7: 5, 8: 5}
    assert G2_SECONDARY == {4: 3, 5: 3, 6: 3, 7: 4, 8: 4}
