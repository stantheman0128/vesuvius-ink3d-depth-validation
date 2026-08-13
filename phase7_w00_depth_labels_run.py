"""Phase 7 run: depth-validated 3D ink labels for w00_20231016151002.

Per PHASE7_PREREGISTRATION_W00_DEPTH_VALIDATED_LABELS_2026-08-13.md. Imports
frozen 5A/5B functions unmodified. Artifact + QC only; no confirmatory gate.

Usage:
    python phase7_w00_depth_labels_run.py --plan   # node/chunk budget only
    python phase7_w00_depth_labels_run.py          # full run (resumable)
"""
import hashlib
import json
import os
import sys

import numpy as np
import tifffile
import zarr
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import find_peaks
from scipy.spatial import cKDTree

from phase5a_depth_attribution_run import (
    CHUNK, REF_SHAPE_ZYX, canonical_to_ref_zyx, load_chain,
)
from phase5b_scaleup_run import CHUNKDIR as SHARED_CHUNKDIR
from phase5b_scaleup_run import REF_BASE, cached_chunk

D7 = "data_phase7_w00_2026-08-13"
INK_TIF = ("data_phase3/official_w00_2d_maps_2026-07-26/ink/phercparis4/"
           "w00_20231016151002/w00_20231016151002_inklabels.tif")
INK_SHA = "31810cd7fd8948618a29232cd0839cf249f0f3906a8218d3408156f315c2d014"
CANON_SHAPE_ZYX = (75784, 32693, 32693)
L_Z, L_MID = 65, 32
FULL_SHAPE = (32249, 51380)

STRIDE = 8
INK_FRAC_MIN = 0.01
NODE_CAP = 800
WIN = 34
PROFILE_HALF = 24
REF_PAD = 8
SIGMA_SMOOTH = 2.0
PROMINENCE_FRAC = 0.25
ZERO_FRAC_CAP = 0.10
CLAMP = 8.0
UNVALIDATED_RADIUS = 24.0
BAND_HALF = 4
NEW_CHUNK_CAP = 2200
RESULT_JSON = "phase7_result_2026-08-13.json"


# ---------------------------------------------------------------- pure rules
def valid_mask(xm, ym, zm):
    m = np.isfinite(xm) & np.isfinite(ym) & np.isfinite(zm)
    return m & (xm >= 0) & (ym >= 0) & (zm >= 0)


def node_normal(xm, ym, zm, r, c):
    """Unit surface normal in (z,y,x) order from map central differences."""
    tc = np.array([(xm[r, c + 1] - xm[r, c - 1]) / 2,
                   (ym[r, c + 1] - ym[r, c - 1]) / 2,
                   (zm[r, c + 1] - zm[r, c - 1]) / 2], float)
    tr = np.array([(xm[r + 1, c] - xm[r - 1, c]) / 2,
                   (ym[r + 1, c] - ym[r - 1, c]) / 2,
                   (zm[r + 1, c] - zm[r - 1, c]) / 2], float)
    n = np.cross(tc, tr)
    norm = np.linalg.norm(n)
    if not np.isfinite(norm) or norm < 1e-9:
        return None
    n = n / norm
    return n[::-1]  # (x,y,z) -> (z,y,x)


def signed_surface_offset(ref_profile):
    """Signed s of the qualifying |gradient| peak nearest to s=0, or None."""
    g = np.abs(np.gradient(np.asarray(ref_profile, float)))
    if g.max() <= 0:
        return None
    peaks, _ = find_peaks(g, prominence=PROMINENCE_FRAC * g.max())
    if len(peaks) == 0:
        return None
    s = peaks - PROFILE_HALF
    return float(s[np.argmin(np.abs(s))])


def stride_subsample(items, cap):
    L = len(items)
    if L <= cap:
        return list(items)
    idx = sorted({(i * L) // cap for i in range(cap)})
    return [items[i] for i in idx]


def correction_and_status(nodes_rc, svals, shape, radius=UNVALIDATED_RADIUS,
                          clamp=CLAMP):
    """Map-resolution correction field and status (0 unval, 1 interp, 2 node)."""
    field = np.zeros(shape, np.float32)
    status = np.zeros(shape, np.uint8)
    if len(nodes_rc) == 0:
        return field, status
    pts = np.asarray(nodes_rc, float)
    vals = np.clip(np.asarray(svals, float), -clamp, clamp)
    rr, cc = np.mgrid[0:shape[0], 0:shape[1]]
    grid = np.stack([rr.ravel(), cc.ravel()], axis=1)
    near = griddata(pts, vals, grid, method="nearest")
    if len(pts) >= 3:
        try:
            interp = griddata(pts, vals, grid, method="linear")
        except Exception:
            interp = near.copy()
    else:
        interp = near.copy()
    interp = np.where(np.isnan(interp), near, interp).reshape(shape)
    dist, _ = cKDTree(pts).query(grid, k=1)
    dist = dist.reshape(shape)
    inside = dist <= radius
    field[inside] = interp[inside].astype(np.float32)
    status[inside] = 1
    for (r, c) in nodes_rc:
        r, c = int(r), int(c)
        if 0 <= r < shape[0] and 0 <= c < shape[1]:
            status[r, c] = 2
    return field, status


def band_block(ink_tile, center_tile, l_z=L_Z, half=BAND_HALF):
    """uint8 (l_z, h, w): 255 in [center-half, center+half] where ink>0."""
    h, w = ink_tile.shape
    out = np.zeros((l_z, h, w), np.uint8)
    centers = np.rint(center_tile).astype(int)
    lo = np.clip(centers - half, 0, l_z - 1)
    hi = np.clip(centers + half, 0, l_z - 1)
    layers = np.arange(l_z)[:, None, None]
    mask = (ink_tile > 0)[None, :, :] & (layers >= lo[None]) & (layers <= hi[None])
    out[mask] = 255
    return out


# ------------------------------------------------------------------ pipeline
def load_maps():
    xm = tifffile.imread(f"{D7}/x.tif")
    ym = tifffile.imread(f"{D7}/y.tif")
    zm = tifffile.imread(f"{D7}/z.tif")
    return xm, ym, zm


def frame_check(xm, ym, zm):
    v = valid_mask(xm, ym, zm)
    inb = (v & (xm < CANON_SHAPE_ZYX[2]) & (ym < CANON_SHAPE_ZYX[1])
           & (zm < CANON_SHAPE_ZYX[0]))
    frac = inb.sum() / max(v.sum(), 1)
    if frac < 0.99:
        raise SystemExit(f"TECHNICAL_FAIL_FRAME: in-bounds fraction {frac:.4f}")
    return v


def ink_cell_fractions(map_shape):
    data = open(INK_TIF, "rb").read()
    if hashlib.sha256(data).hexdigest() != INK_SHA:
        raise SystemExit("TECHNICAL_FAIL_IDENTITY: inklabels sha mismatch")
    ink = tifffile.imread(INK_TIF)
    assert ink.shape == FULL_SHAPE, ink.shape
    H, W = map_shape
    sy, sx = FULL_SHAPE[0] / H, FULL_SHAPE[1] / W
    rows = np.round(np.arange(0, H + 1, STRIDE) * sy).astype(int)
    cols = np.round(np.arange(0, W + 1, STRIDE) * sx).astype(int)
    frac = {}
    for i in range(len(rows) - 1):
        y0, y1 = rows[i], rows[i + 1]
        band = (ink[y0:y1] > 0)
        colsum = np.add.reduceat(band.sum(axis=0, dtype=np.int64), cols[:-1])
        for j in range(len(cols) - 1):
            area = (y1 - y0) * (cols[j + 1] - cols[j])
            frac[(i * STRIDE, j * STRIDE)] = colsum[j] / max(area, 1)
    return ink, frac


def node_in_ref(p_zyx, Ai, t):
    corners = np.array([[p_zyx[0] + a, p_zyx[1] + b, p_zyx[2] + d]
                        for a in (-WIN, WIN) for b in (-WIN, WIN)
                        for d in (-WIN, WIN)], float)
    q = canonical_to_ref_zyx(corners, Ai, t)
    return bool(np.all(q - REF_PAD >= 0)
                and np.all(q + REF_PAD <= np.array(REF_SHAPE_ZYX) - 1))


def node_chunks(p_zyx, Ai, t):
    corners = np.array([[p_zyx[0] + a, p_zyx[1] + b, p_zyx[2] + d]
                        for a in (-WIN, WIN) for b in (-WIN, WIN)
                        for d in (-WIN, WIN)], float)
    q = canonical_to_ref_zyx(corners, Ai, t)
    lo = np.maximum(np.floor(q.min(axis=0) - REF_PAD).astype(int), 0)
    hi = np.minimum(np.ceil(q.max(axis=0) + REF_PAD).astype(int),
                    np.array(REF_SHAPE_ZYX) - 1)
    k0, k1 = lo // CHUNK, hi // CHUNK
    return [(cz, cy, cx)
            for cz in range(k0[0], k1[0] + 1)
            for cy in range(k0[1], k1[1] + 1)
            for cx in range(k0[2], k1[2] + 1)]


def measure_node(p_zyx, normal_zyx, Ai, t, chunk_store, journal):
    n = 2 * WIN + 1
    off = np.stack(np.meshgrid(*(np.arange(n) - WIN,) * 3, indexing="ij"),
                   axis=-1).reshape(-1, 3)
    q = canonical_to_ref_zyx(p_zyx + off, Ai, t)
    ks = node_chunks(p_zyx, Ai, t)
    arr = np.array(ks)
    k0 = arr.min(axis=0)
    k1 = arr.max(axis=0)
    sub = np.zeros(tuple((k1 - k0 + 1) * CHUNK), np.uint8)
    for (cz, cy, cx) in ks:
        blk = chunk_store.get((cz, cy, cx))
        if blk is None:
            blk = cached_chunk(REF_BASE, SHARED_CHUNKDIR, cz, cy, cx, journal)
            chunk_store[(cz, cy, cx)] = blk
            if len(chunk_store) > 600:
                chunk_store.pop(next(iter(chunk_store)))
        oz, oy, ox = (np.array([cz, cy, cx]) - k0) * CHUNK
        sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk
    window = map_coordinates(sub.astype(np.float32),
                             (q - k0 * CHUNK).T, order=1,
                             mode="constant", cval=0.0).reshape(n, n, n)
    zero_frac = float(np.mean(window < 0.5))
    R = gaussian_filter(window, SIGMA_SMOOTH)
    s = np.arange(-PROFILE_HALF, PROFILE_HALF + 1, dtype=float)
    pts = WIN + s[:, None] * np.asarray(normal_zyx)[None, :]
    prof = map_coordinates(R, pts.T, order=1, mode="constant", cval=0.0)
    s_star = signed_surface_offset(prof)
    ok = zero_frac <= ZERO_FRAC_CAP and s_star is not None
    return {"zero_frac": round(zero_frac, 4), "s_star": s_star,
            "evaluable": bool(ok)}


def main():
    plan_only = "--plan" in sys.argv
    os.makedirs(D7, exist_ok=True)
    A, Ai, t = load_chain()
    xm, ym, zm = load_maps()
    v = frame_check(xm, ym, zm)
    H, W = xm.shape

    sel_path = f"{D7}/node_selection.json"
    ink = None
    if os.path.exists(sel_path):
        sel = json.load(open(sel_path))
        nodes = [tuple(x) for x in sel["nodes"]]
    else:
        ink, frac = ink_cell_fractions((H, W))
        eligible = []
        for r in range(0, H - 1, STRIDE):
            for c in range(0, W - 1, STRIDE):
                if r == 0 or c == 0:
                    continue
                if frac.get((r, c), 0.0) < INK_FRAC_MIN:
                    continue
                if not (v[r, c] and v[r - 1, c] and v[r + 1, c]
                        and v[r, c - 1] and v[r, c + 1]):
                    continue
                p = np.array([zm[r, c], ym[r, c], xm[r, c]], float)
                if not node_in_ref(p, Ai, t):
                    continue
                eligible.append((r, c))
        nodes = stride_subsample(eligible, NODE_CAP)
        json.dump({"eligible_count": len(eligible),
                   "nodes": [list(n) for n in nodes]}, open(sel_path, "w"))
    print(f"nodes: {len(nodes)}")

    all_chunks = set()
    for (r, c) in nodes:
        p = np.array([zm[r, c], ym[r, c], xm[r, c]], float)
        all_chunks.update(node_chunks(p, Ai, t))
    new = sum(1 for (cz, cy, cx) in all_chunks if not os.path.exists(
        f"{SHARED_CHUNKDIR}/chunk_{cz}_{cy}_{cx}.bin"))
    print(f"chunk plan: {len(all_chunks)} total, {new} new "
          f"(~{new * CHUNK**3 / 1e6:.0f} MB)")
    if new > NEW_CHUNK_CAP:
        raise SystemExit("TECHNICAL_FAIL_BUDGET")
    if plan_only:
        return

    meas_path = f"{D7}/node_measurements.json"
    if os.path.exists(meas_path):
        meas = json.load(open(meas_path))
    else:
        journal, store, meas = [], {}, {}
        for i, (r, c) in enumerate(nodes, 1):
            p = np.array([zm[r, c], ym[r, c], xm[r, c]], float)
            nrm = node_normal(xm, ym, zm, r, c)
            if nrm is None:
                meas[f"{r}_{c}"] = {"evaluable": False, "s_star": None,
                                    "zero_frac": None, "reason": "degenerate_normal"}
                continue
            meas[f"{r}_{c}"] = measure_node(p, nrm, Ai, t, store, journal)
            if i % 100 == 0 or i == len(nodes):
                print(f"  measured {i}/{len(nodes)}", flush=True)
        json.dump(meas, open(meas_path, "w"))
        json.dump(journal, open(f"{D7}/ref_chunk_journal.json", "w"))

    ev_nodes, ev_vals = [], []
    for key, m in meas.items():
        if m.get("evaluable"):
            r, c = map(int, key.split("_"))
            ev_nodes.append((r, c))
            ev_vals.append(m["s_star"])
    n_eval = len(ev_nodes)
    print(f"evaluable nodes: {n_eval}/{len(nodes)}")

    field, status = correction_and_status(ev_nodes, ev_vals, (H, W))

    if ink is None:
        ink = tifffile.imread(INK_TIF)

    lab_path = f"{D7}/w00_20231016151002_inklabels_v5.zarr"
    qc_path = f"{D7}/w00_20231016151002_inklabels_v5_qc.zarr"
    g = zarr.open_group(qc_path, mode="w")
    g.create_array("s_star", shape=field.shape, dtype="f4")[:] = field
    g.create_array("status", shape=status.shape, dtype="u1")[:] = status
    dmap = np.full((H, W), np.nan, np.float32)
    for (r, c), val in zip(ev_nodes, ev_vals):
        dmap[r, c] = abs(val)
    g.create_array("D_node", shape=dmap.shape, dtype="f4")[:] = dmap

    za = zarr.open_array(lab_path, mode="w", shape=(L_Z,) + FULL_SHAPE,
                         chunks=(L_Z, 256, 256), dtype="u1", fill_value=0)
    sy, sx = FULL_SHAPE[0] / H, FULL_SHAPE[1] / W
    tiles_written = 0
    for y0 in range(0, FULL_SHAPE[0], 256):
        for x0 in range(0, FULL_SHAPE[1], 256):
            y1, x1 = min(y0 + 256, FULL_SHAPE[0]), min(x0 + 256, FULL_SHAPE[1])
            tile = ink[y0:y1, x0:x1]
            if not tile.any():
                continue
            yy = np.clip((np.arange(y0, y1) / sy).astype(int), 0, H - 1)
            xx = np.clip((np.arange(x0, x1) / sx).astype(int), 0, W - 1)
            centers = L_MID + field[np.ix_(yy, xx)]
            za[:, y0:y1, x0:x1] = band_block(tile, centers)
            tiles_written += 1
    print(f"label tiles written: {tiles_written}")

    sv = np.array(ev_vals, float)
    result = {
        "preregistration":
            "PHASE7_PREREGISTRATION_W00_DEPTH_VALIDATED_LABELS_2026-08-13.md",
        "n_nodes": len(nodes), "n_eval": n_eval,
        "s_star_quartiles": ([float(x) for x in np.percentile(sv, [25, 50, 75])]
                             if n_eval else None),
        "D_le_3_fraction": (float(np.mean(np.abs(sv) <= 3)) if n_eval else None),
        "tiles_written": tiles_written,
        "label_zarr": lab_path, "qc_zarr": qc_path,
    }
    # validated ink fraction, computed directly
    stat_full = status[np.ix_(
        np.clip((np.arange(FULL_SHAPE[0]) / sy).astype(int), 0, H - 1),
        np.clip((np.arange(FULL_SHAPE[1]) / sx).astype(int), 0, W - 1))]
    result["validated_ink_fraction"] = float((stat_full[ink > 0] > 0).mean())
    with open(RESULT_JSON, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
