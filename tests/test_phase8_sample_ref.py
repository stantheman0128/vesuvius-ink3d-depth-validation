"""Phase 8 sample_ref must trilinear-interpolate across chunk faces.

Independent expected value: assemble adjacent 128^3 blocks into one array and
map_coordinates on that assembly. Sampling each block alone returns 0 for any
local coordinate > 127 (scipy constant-mode OOB), which is not trilinear.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phase5a_depth_attribution_run import CHUNK  # noqa: E402
from phase8_khj1222_band_score import sample_ref  # noqa: E402


def test_sample_ref_crosses_chunk_face_with_known_values():
    cz = cy = 10
    cx0, cx1 = 5, 6
    blk0 = np.full((CHUNK, CHUNK, CHUNK), 100, np.uint8)
    blk1 = np.full((CHUNK, CHUNK, CHUNK), 200, np.uint8)
    store = {(cz, cy, cx0): blk0, (cz, cy, cx1): blk1}

    gz = cz * CHUNK
    gy = cy * CHUNK
    # Midway between last voxel of cx0 and first voxel of cx1.
    gx = cx0 * CHUNK + 127.5
    pt = np.array([[gz, gy, gx]], float)

    sub = np.zeros((CHUNK, CHUNK, 2 * CHUNK), np.float32)
    sub[:, :, :CHUNK] = blk0
    sub[:, :, CHUNK:] = blk1
    sub_origin = np.array([cz * CHUNK, cy * CHUNK, cx0 * CHUNK], float)
    expect = float(map_coordinates(sub, (pt - sub_origin).T, order=1,
                                   mode="constant", cval=0.0)[0])
    assert expect == 150.0

    got = sample_ref(pt, store)[0]
    assert abs(got - expect) < 1e-6


def test_sample_ref_interior_unchanged():
    cz = cy = cx = 3
    blk = np.arange(CHUNK ** 3, dtype=np.uint8).reshape(CHUNK, CHUNK, CHUNK)
    store = {(cz, cy, cx): blk}
    local = np.array([40.25, 50.5, 60.75])
    pt = (np.array([cz, cy, cx]) * CHUNK + local)[None].astype(float)
    expect = float(map_coordinates(blk.astype(np.float32), local[:, None],
                                   order=1, mode="constant", cval=0.0)[0])
    got = sample_ref(pt, store)[0]
    assert abs(got - expect) < 1e-6
