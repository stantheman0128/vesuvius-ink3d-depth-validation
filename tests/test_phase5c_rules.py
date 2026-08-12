"""Unit tests required by PHASE5C preregistration: the prescreen predicate."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phase5c_scaleup_run import CONTENT_MIN_FRACTION, content_ok  # noqa: E402


def test_content_ok_accepts_full_chunk():
    assert content_ok(np.full((128, 128, 128), 90, np.uint8))


def test_content_ok_rejects_masked_air():
    assert not content_ok(np.zeros((128, 128, 128), np.uint8))


def test_content_ok_boundary_is_inclusive():
    arr = np.zeros((100,), np.uint8)
    arr[:50] = 7  # exactly the 0.5 fraction
    assert CONTENT_MIN_FRACTION == 0.5
    assert content_ok(arr)
    arr[49] = 0  # 0.49: below threshold
    assert not content_ok(arr)
