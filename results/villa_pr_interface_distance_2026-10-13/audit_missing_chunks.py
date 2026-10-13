"""Which reference chunks do the w00 runs need, and which are absent from the local mirror?

Reads the per-point reference centres from the sigma-0 points run, derives the
read box each point needs (window 34 point-voxels * scale, padded), lists the
128^3 chunk keys it touches, and compares with the chunk files present in
local_1129_level0. Also estimates the smaller set needed if only the profile
tube (49 samples along the supplied normal, plus smoothing support) were read.

Output: JSON with per-point missing counts and totals, printed and written next
to this file. Pure arithmetic, no network, no zarr.
"""
import json
import os
from collections import defaultdict

import numpy as np

STG = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(STG, "w00_interface_distance_sigma0.json")
MIRROR = os.path.join(STG, "local_1129_level0")
CSV = (r"C:\Users\stans\Projects\vesuvius-ink3d-depth-validation"
       r"\data_phase8_khj1222_2026-08-25\w00_20231016151002_depth_anchors.csv")
CHUNK = 128
WINDOW = 34
PROFILE_HALF = 24
SMOOTH_SUPPORT = 3 * 2.0  # 3 sigma at sigma 2, in point voxels

with open(RUN) as f:
    run = json.load(f)
scale = run["reference_voxels_per_point_voxel"]
M = np.array(run["chain_matrix_xyz"], dtype=float)
shape = np.array(run["reference_shape"])

present = {k for k in os.listdir(MIRROR) if k.count(".") == 2}
present_nonempty = {k for k in present if os.path.getsize(os.path.join(MIRROR, k)) > 0}

normals = {}
with open(CSV) as f:
    header = f.readline().strip().split(",")
    idx = {h: i for i, h in enumerate(header)}
    for line in f:
        p = line.strip().split(",")
        if p[idx["sampled_from"]] != "center":
            continue
        pid = f"{p[idx['cell_row']]}_{p[idx['cell_col']]}"
        n = np.array([float(p[idx["normal_z"]]), float(p[idx["normal_y"]]), float(p[idx["normal_x"]])])
        normals[pid] = n


def keys_for_box(lo, hi):
    lo = np.maximum(np.floor(lo).astype(int), 0)
    hi = np.minimum(np.ceil(hi).astype(int), shape - 1)
    if np.any(hi < lo):
        return set()
    c0 = lo // CHUNK
    c1 = hi // CHUNK
    return {f"{z}.{y}.{x}"
            for z in range(c0[0], c1[0] + 1)
            for y in range(c0[1], c1[1] + 1)
            for x in range(c0[2], c1[2] + 1)}


def keys_for_tube(centre, n_ref, half_ref, radius_ref):
    keys = set()
    for s in np.arange(-half_ref, half_ref + 1e-9, CHUNK / 4):
        c = centre + s * n_ref
        keys |= keys_for_box(c - radius_ref, c + radius_ref)
    return keys


rows = []
for rec in run["records"]:
    if "outside_reference" in rec["non_evaluable_reasons"]:
        continue
    centre = np.array(rec["reference_zyx"])
    half_box = WINDOW * scale * 1.03  # axis-aligned bound of the mapped cube; off-diagonals are ~1.3% of scale
    box_keys = keys_for_box(centre - half_box - 2, centre + half_box + 3)
    n = normals.get(rec["id"])
    tube_keys = None
    if n is not None and np.all(np.isfinite(n)) and np.linalg.norm(n) > 0:
        # map the normal direction through the linear part (ZYX in, ZYX out)
        A = M[:, :3][::-1, ::-1]
        n_ref = A @ (n / np.linalg.norm(n))
        n_ref /= np.linalg.norm(n_ref)
        tube_keys = keys_for_tube(centre, n_ref, PROFILE_HALF * scale, SMOOTH_SUPPORT * scale + 2)
    rows.append({
        "id": rec["id"],
        "evaluable": rec["evaluable"],
        "D": rec["D_voxels"],
        "box_keys": len(box_keys),
        "box_missing": sorted(box_keys - present_nonempty),
        "tube_keys": None if tube_keys is None else len(tube_keys),
        "tube_missing": None if tube_keys is None else sorted(tube_keys - present_nonempty),
    })

all_box_missing = set()
all_tube_missing = set()
for r in rows:
    all_box_missing |= set(r["box_missing"])
    if r["tube_missing"] is not None:
        all_tube_missing |= set(r["tube_missing"])

summary = {
    "points_in_footprint": len(rows),
    "mirror_chunks_present_nonempty": len(present_nonempty),
    "mirror_chunks_present_total_files": len(present),
    "points_with_any_box_chunk_missing": sum(1 for r in rows if r["box_missing"]),
    "points_with_any_tube_chunk_missing": sum(1 for r in rows if r["tube_missing"]),
    "evaluable_points_with_tube_chunk_missing": sum(1 for r in rows if r["evaluable"] and r["tube_missing"]),
    "distinct_box_chunks_missing": len(all_box_missing),
    "distinct_tube_chunks_missing": len(all_tube_missing),
    "est_download_MB_box_at_1MB_per_chunk": len(all_box_missing),
    "est_download_MB_tube_at_1MB_per_chunk": len(all_tube_missing),
}
out = {"summary": summary, "rows": rows}
with open(os.path.join(STG, "audit_missing_chunks.json"), "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps(summary, indent=1))
