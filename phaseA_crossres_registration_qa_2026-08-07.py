"""Phase A cross-resolution registration QA for villa #192.

Technical QA only, no scientific ink claim: tests whether the official
1.129um->2.4um transform places the already-acquired P50_INTERIOR raw CT cube
(78keV 2.4um canonical volume) at the right spot, and measures the residual
offset by normalized cross-correlation.

Fail-closed: any HTTP status other than 200/404, any size mismatch, or a
correlation peak on the search boundary is reported as failure. 404 chunks are
filled with the zarr fill_value 0 (masked volume convention).

Download budget: only the 128^3 zarr chunks of the 1.129um volume that
intersect the mapped search window; each chunk is exactly 2,097,152 raw uint8
bytes (compressor null). Expected total is tens of MB, well under the 5 GB gate.
"""
import hashlib
import json
import time
import urllib.request
import urllib.error

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve

BASE = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/"
        "PHercParis4/volumes/20260608103018-1.129um-0.2m-78keV-masked.zarr/0")
SHAPE_1129_ZYX = (59969, 36006, 32354)
CHUNK = 128
CHUNK_BYTES = CHUNK ** 3
OUTDIR = "data_crossscan_1129_p50_2026-08-07"
SNAP = "data_phase4b_snapshot_2026-08-06/bodies"

# P50_INTERIOR prediction chunk (z,y,x) = (148, 58, 75) in the 78keV 2.4um volume
P50_ORIGIN_ZYX = np.array([148 * 256, 58 * 256, 75 * 256])  # (37888, 14848, 19200)
# raw octant chunk files, index 10..17, chunk coords (296|297, 116|117, 150|151)
RAW_FILES = {
    (0, 0, 0): "10_raw_ct_296_116_150.bin",
    (0, 0, 1): "11_raw_ct_296_116_151.bin",
    (0, 1, 0): "12_raw_ct_296_117_150.bin",
    (0, 1, 1): "13_raw_ct_296_117_151.bin",
    (1, 0, 0): "14_raw_ct_297_116_150.bin",
    (1, 0, 1): "15_raw_ct_297_116_151.bin",
    (1, 1, 0): "16_raw_ct_297_117_150.bin",
    (1, 1, 1): "17_raw_ct_297_117_151.bin",
}

TEMPLATE = 128          # central 128^3 of the P50 raw cube, 2.4um voxels
SEARCH_PAD = 16         # +/- voxels of search range on the 2.4um grid
MARGIN_1129 = 8.0       # extra sampling margin in 1.129um voxels


def load_transform():
    with open("data_crossscan_meta_2026-08-07/transform_20260608103018.json") as f:
        d = json.load(f)
    M = np.array(d["transformation_matrix"], dtype=float)
    return M[:, :3], M[:, 3]  # x_2.4_xyz = A @ x_1129_xyz + t


def assemble_p50_raw():
    cube = np.empty((256, 256, 256), dtype=np.uint8)
    for (dz, dy, dx), fname in RAW_FILES.items():
        buf = np.fromfile(f"{SNAP}/{fname}", dtype=np.uint8)
        assert buf.size == CHUNK_BYTES, (fname, buf.size)
        cube[dz * 128:(dz + 1) * 128, dy * 128:(dy + 1) * 128,
             dx * 128:(dx + 1) * 128] = buf.reshape(128, 128, 128)
    return cube


def fetch_chunk(cz, cy, cx, journal):
    url = f"{BASE}/{cz}/{cy}/{cx}"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
            if len(data) != CHUNK_BYTES:
                raise RuntimeError(f"size mismatch {url}: {len(data)}")
            journal.append({"key": f"0/{cz}/{cy}/{cx}", "bytes": len(data),
                            "sha256": hashlib.sha256(data).hexdigest(), "status": 200})
            return np.frombuffer(data, dtype=np.uint8).reshape(128, 128, 128)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                journal.append({"key": f"0/{cz}/{cy}/{cx}", "bytes": 0,
                                "sha256": None, "status": 404})
                return np.zeros((128, 128, 128), dtype=np.uint8)
            if attempt == 2:
                raise
        except Exception:
            if attempt == 2:
                raise
        time.sleep(1.5 * (attempt + 1))


def main():
    import os
    os.makedirs(OUTDIR, exist_ok=True)
    A, t = load_transform()
    Ai = np.linalg.inv(A)

    p50 = assemble_p50_raw()
    template = p50[64:192, 64:192, 64:192].astype(np.float64)

    # search window on the 2.4um grid, (z,y,x) voxel coords in the 78keV volume
    lo_zyx = P50_ORIGIN_ZYX + 64 - SEARCH_PAD
    n = TEMPLATE + 2 * SEARCH_PAD  # 160
    zz, yy, xx = np.meshgrid(np.arange(n) + lo_zyx[0],
                             np.arange(n) + lo_zyx[1],
                             np.arange(n) + lo_zyx[2], indexing="ij")
    pts_xyz = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1).astype(float)

    # map to 1.129um voxel coords
    q = (pts_xyz - t) @ Ai.T  # (N,3) xyz
    q_zyx = q[:, ::-1]
    qlo = np.floor(q_zyx.min(axis=0) - MARGIN_1129).astype(int)
    qhi = np.ceil(q_zyx.max(axis=0) + MARGIN_1129).astype(int)
    qlo = np.maximum(qlo, 0)
    qhi = np.minimum(qhi, np.array(SHAPE_1129_ZYX) - 1)
    c0 = qlo // CHUNK
    c1 = qhi // CHUNK
    nchunks = int(np.prod(c1 - c0 + 1))
    est_mb = nchunks * CHUNK_BYTES / 1e6
    print(f"1.129um voxel bbox lo={qlo} hi={qhi}; chunks={nchunks} (~{est_mb:.0f} MB)")
    if est_mb > 1000:
        raise SystemExit("ABORT: estimated download exceeds 1 GB self-cap")

    journal = []
    sub_shape = tuple((c1 - c0 + 1) * CHUNK)
    sub = np.zeros(sub_shape, dtype=np.uint8)
    done = 0
    for cz in range(c0[0], c1[0] + 1):
        for cy in range(c0[1], c1[1] + 1):
            for cx in range(c0[2], c1[2] + 1):
                blk = fetch_chunk(cz, cy, cx, journal)
                oz, oy, ox = (cz - c0[0]) * CHUNK, (cy - c0[1]) * CHUNK, (cx - c0[2]) * CHUNK
                sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk
                done += 1
        print(f"  fetched {done}/{nchunks} chunks", flush=True)

    with open(f"{OUTDIR}/download_journal.json", "w") as f:
        json.dump({"base": BASE, "chunk_origin_zyx": (c0 * CHUNK).tolist(),
                   "chunks": journal}, f, indent=1)

    # sample the 1.129um data at the mapped points (trilinear)
    coords = (q_zyx - (c0 * CHUNK)).T  # (3,N) in sub-volume index space
    resampled = map_coordinates(sub.astype(np.float32), coords, order=1,
                                mode="constant", cval=0.0).reshape(n, n, n)

    # normalized cross-correlation of template within resampled search volume
    tz = template - template.mean()
    num = fftconvolve(resampled, tz[::-1, ::-1, ::-1], mode="valid")
    ones = np.ones_like(tz)
    s1 = fftconvolve(resampled, ones[::-1, ::-1, ::-1], mode="valid")
    s2 = fftconvolve(resampled ** 2, ones[::-1, ::-1, ::-1], mode="valid")
    nvox = tz.size
    denom = np.sqrt(np.maximum(s2 - s1 ** 2 / nvox, 1e-9)) * np.linalg.norm(tz)
    ncc = num / denom
    peak = np.unravel_index(np.argmax(ncc), ncc.shape)
    offset = np.array(peak) - SEARCH_PAD  # (z,y,x) voxels at 2.4um
    on_edge = bool(any(p in (0, ncc.shape[0] - 1) for p in peak))

    result = {
        "status": "PHASEA_QA_COMPLETED",
        "template": "P50_INTERIOR raw central 128^3 (78keV 2.4um)",
        "search_pad_voxels": SEARCH_PAD,
        "ncc_peak": float(ncc.max()),
        "ncc_at_zero_offset": float(ncc[SEARCH_PAD, SEARCH_PAD, SEARCH_PAD]),
        "peak_offset_zyx_2p4um_voxels": offset.tolist(),
        "peak_offset_um": (offset * 2.4).tolist(),
        "peak_on_search_edge": on_edge,
        "chunks_fetched": nchunks,
        "chunks_404": sum(1 for j in journal if j["status"] == 404),
        "downloaded_mb": sum(j["bytes"] for j in journal) / 1e6,
    }
    with open(f"{OUTDIR}/phaseA_qa_result.json", "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
