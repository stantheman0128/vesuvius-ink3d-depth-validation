"""Phase 8: D/FWHM of khj1222's w00 measured band against the 1.129um scan.

Per PHASE8_PREREGISTRATION_KHJ1222_W00_BAND_2026-08-25.md.
Uses the frozen 5A chain and 5B chunk cache helper. Does not flip normals.

Usage:
    python phase8_khj1222_band_score.py --plan
    python phase8_khj1222_band_score.py
"""
import csv
import hashlib
import json
import os
import sys
import urllib.request

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.signal import find_peaks

from phase5a_depth_attribution_run import (
    CHUNK, CHUNK_BYTES, FWHM_CENSOR, PROFILE_HALF, PROMINENCE_FRAC,
    REF_SHAPE_ZYX, ZERO_FRAC_CAP, canonical_to_ref_zyx, load_chain,
)
from phase5b_scaleup_run import CHUNKDIR as SHARED_CHUNKDIR
from phase5b_scaleup_run import REF_BASE, cached_chunk

CSV_URL = ("https://raw.githubusercontent.com/khj1222/vesuvius-challenge/"
           "main/submission/depth_anchors/w00_20231016151002_depth_anchors.csv")
JSON_URL = ("https://raw.githubusercontent.com/khj1222/vesuvius-challenge/"
            "main/submission/depth_anchors/w00_20231016151002_depth_anchors.json")
OUTDIR = "data_phase8_khj1222_2026-08-25"
RESULT_JSON = "results/phase8_result_2026-08-25.json"
ANCHOR_CAP = 400
CHUNK_CAP = 2200
D_THRESHOLD = 3.0


def download(url, dest):
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    if os.path.exists(dest):
        data = open(dest, "rb").read()
    else:
        data = urllib.request.urlopen(url, timeout=120).read()
        with open(dest, "wb") as f:
            f.write(data)
    return data


def load_center_rows(csv_path):
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            if rec.get("sampled_from") != "center":
                continue
            rows.append({
                "cell_row": int(rec["cell_row"]),
                "cell_col": int(rec["cell_col"]),
                "region": rec["region"],
                "center_layer": float(rec["center_layer"]),
                "offset_from_plane": float(rec["offset_from_plane"]),
                "half_width_layers": float(rec["half_width_layers"]),
                "base": np.array([float(rec["base_x"]),
                                  float(rec["base_y"]),
                                  float(rec["base_z"])], float),
                "normal": np.array([float(rec["normal_x"]),
                                    float(rec["normal_y"]),
                                    float(rec["normal_z"])], float),
            })
    rows.sort(key=lambda r: (r["cell_row"], r["cell_col"]))
    return rows


def stride_subsample(items, cap):
    L = len(items)
    if L <= cap:
        return list(items)
    idx = sorted({(i * L) // cap for i in range(cap)})
    return [items[i] for i in idx]


def unit(v):
    n = np.linalg.norm(v)
    if not np.isfinite(n) or n < 1e-9:
        return None
    return v / n


def band_center_xyz(row):
    n = unit(row["normal"])
    if n is None:
        return None, None
    return row["base"] + row["offset_from_plane"] * n, n


def profile_canonical_zyx(center_xyz, normal_xyz):
    s = np.arange(-PROFILE_HALF, PROFILE_HALF + 1, dtype=float)
    xyz = center_xyz[None, :] + s[:, None] * normal_xyz[None, :]
    return xyz[:, ::-1], s


def fwhm_abs_grad(g, peak_idx):
    p = np.asarray(g, float)
    half = (p.max() + p.min()) / 2.0
    left = right = None
    for i in range(peak_idx, 0, -1):
        if p[i - 1] < half <= p[i]:
            left = (i - 1) + (half - p[i - 1]) / (p[i] - p[i - 1])
            break
    for i in range(peak_idx, len(p) - 1):
        if p[i + 1] < half <= p[i]:
            right = i + (p[i] - half) / (p[i] - p[i + 1])
            break
    if left is None or right is None:
        return FWHM_CENSOR, True
    return float(right - left), False


def measure_profile(ref_prof):
    rec = {}
    zero_frac = float(np.mean(np.asarray(ref_prof) < 0.5))
    rec["masked_zero_fraction"] = round(zero_frac, 4)
    g = np.abs(np.gradient(np.asarray(ref_prof, float)))
    rec["evaluable"] = False
    rec["D_voxels"] = None
    rec["fwhm_voxels"] = None
    rec["fwhm_censored"] = True
    rec["non_evaluable_reasons"] = []
    if zero_frac > ZERO_FRAC_CAP:
        rec["non_evaluable_reasons"].append("masked_zero_fraction_exceeds_cap")
    if g.max() <= 0:
        rec["non_evaluable_reasons"].append("no_qualifying_surface")
        return rec
    peaks, _ = find_peaks(g, prominence=PROMINENCE_FRAC * g.max())
    if len(peaks) == 0:
        rec["non_evaluable_reasons"].append("no_qualifying_surface")
        return rec
    s = peaks - PROFILE_HALF
    nearest = int(np.argmin(np.abs(s)))
    rec["D_voxels"] = round(float(np.abs(s[nearest])), 2)
    fwhm, censored = fwhm_abs_grad(g, int(peaks[nearest]))
    rec["fwhm_voxels"] = round(fwhm, 2)
    rec["fwhm_censored"] = censored
    if zero_frac <= ZERO_FRAC_CAP:
        rec["evaluable"] = True
        rec["non_evaluable_reasons"] = []
    return rec


def chunk_keys_for_points(pts_zyx, Ai, t):
    """Chunks covering floor and ceil of each inside sample (trilinear stencil)."""
    q = canonical_to_ref_zyx(pts_zyx, Ai, t)
    shape = np.array(REF_SHAPE_ZYX)
    inside = np.all(q >= 0, axis=1) & np.all(q <= (shape - 1), axis=1)
    if not np.any(inside):
        return set(), q, inside
    q_in = q[inside]
    flo = np.floor(q_in).astype(int)
    cei = np.minimum(np.ceil(q_in).astype(int), shape - 1)
    keys = set()
    for pts in (flo, cei):
        qc = pts // CHUNK
        keys.update(tuple(int(x) for x in cc) for cc in np.unique(qc, axis=0))
    return keys, q, inside


def assemble_store(chunks, journal):
    store = {}
    for n, (cz, cy, cx) in enumerate(sorted(chunks), 1):
        store[(cz, cy, cx)] = cached_chunk(REF_BASE, SHARED_CHUNKDIR, cz, cy, cx,
                                           journal)
        if n % 25 == 0 or n == len(chunks):
            print(f"  fetched {n}/{len(chunks)}", flush=True)
    return store


def sample_ref(q, store):
    """Trilinear sample of reference voxels at q (N,3) zyx.

    Assembles a sub-volume from cached chunks so order-1 interpolation can
    cross chunk faces. Sampling each 128^3 block alone returns cval=0 for any
    local coordinate > 127 (fractional samples on the high face).
    """
    q = np.asarray(q, float)
    vals = np.zeros(len(q), float)
    shape = np.array(REF_SHAPE_ZYX)
    inside = np.all(q >= 0, axis=1) & np.all(q <= (shape - 1), axis=1)
    if not np.any(inside):
        return vals
    q_in = q[inside]
    flo = np.floor(q_in.min(axis=0)).astype(int)
    cei = np.minimum(np.ceil(q_in.max(axis=0)).astype(int), shape - 1)
    k0 = flo // CHUNK
    k1 = cei // CHUNK
    sub = np.zeros(tuple((k1 - k0 + 1) * CHUNK), np.float32)
    for cz in range(int(k0[0]), int(k1[0]) + 1):
        for cy in range(int(k0[1]), int(k1[1]) + 1):
            for cx in range(int(k0[2]), int(k1[2]) + 1):
                blk = store.get((cz, cy, cx))
                if blk is None:
                    continue
                oz, oy, ox = (np.array([cz, cy, cx]) - k0) * CHUNK
                sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk
    sub_origin = k0 * CHUNK
    vals[inside] = map_coordinates(sub, (q_in - sub_origin).T, order=1,
                                   mode="constant", cval=0.0)
    return vals


def main():
    plan_only = "--plan" in sys.argv
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs("results", exist_ok=True)
    csv_path = f"{OUTDIR}/w00_20231016151002_depth_anchors.csv"
    json_path = f"{OUTDIR}/w00_20231016151002_depth_anchors.json"
    csv_bytes = download(CSV_URL, csv_path)
    json_bytes = download(JSON_URL, json_path)
    sidecar = json.loads(json_bytes.decode("utf-8"))
    if sidecar.get("assumptions", {}).get("layer_step_voxels") != 1.0:
        raise SystemExit("TECHNICAL_FAIL_SIDECAR: layer_step_voxels changed")

    rows = load_center_rows(csv_path)
    A, Ai, t = load_chain()
    eligible = []
    skipped_normal = 0
    for row in rows:
        center, normal = band_center_xyz(row)
        if center is None:
            skipped_normal += 1
            continue
        pts, s = profile_canonical_zyx(center, normal)
        keys, q, inside = chunk_keys_for_points(pts, Ai, t)
        if not np.any(inside):
            continue
        eligible.append((row, center, normal, pts, q, inside, keys))
    selected = stride_subsample(eligible, ANCHOR_CAP)
    chunks = set()
    prepared = []
    for item in selected:
        prepared.append(item)
        chunks |= item[-1]

    est_mb = len(chunks) * CHUNK_BYTES / 1e6
    print(f"center rows: {len(rows)}; eligible in 1.129um footprint: "
          f"{len(eligible)}; selected: {len(selected)}; "
          f"bad normals: {skipped_normal}")
    print(f"chunk plan: {len(chunks)} chunks (~{est_mb:.0f} MB)")
    print(f"csv sha256: {hashlib.sha256(csv_bytes).hexdigest()}")
    print(f"json sha256: {hashlib.sha256(json_bytes).hexdigest()}")
    if len(chunks) > CHUNK_CAP:
        raise SystemExit("TECHNICAL_FAIL_BUDGET: chunk plan exceeds self-cap")
    if plan_only:
        with open(f"{OUTDIR}/plan.json", "w") as f:
            json.dump({
                "n_center": len(rows), "n_eligible": len(eligible),
                "n_selected": len(selected),
                "n_chunks": len(chunks), "est_mb": round(est_mb, 1),
                "csv_sha256": hashlib.sha256(csv_bytes).hexdigest(),
                "json_sha256": hashlib.sha256(json_bytes).hexdigest(),
            }, f, indent=1)
        return

    journal = []
    store = assemble_store(chunks, journal)
    records = []
    for row, center, normal, pts, q, inside, keys in prepared:
        rec = {
            "cell_row": row["cell_row"], "cell_col": row["cell_col"],
            "region": row["region"],
            "center_layer": row["center_layer"],
            "half_width_layers": row["half_width_layers"],
            "offset_from_plane": row["offset_from_plane"],
        }
        if not np.any(inside):
            rec["evaluable"] = False
            rec["non_evaluable_reasons"] = ["outside_reference_footprint"]
            rec["D_voxels"] = None
            rec["fwhm_voxels"] = None
            records.append(rec)
            continue
        prof = sample_ref(q, store)
        rec.update(measure_profile(prof))
        records.append(rec)

    ev = [r for r in records if r.get("evaluable")]
    ds = [r["D_voxels"] for r in ev if r["D_voxels"] is not None]
    k = sum(1 for d in ds if d <= D_THRESHOLD)
    result = {
        "preregistration": "PHASE8_PREREGISTRATION_KHJ1222_W00_BAND_2026-08-25.md",
        "n_center": len(rows), "n_selected": len(selected),
        "n_records": len(records), "n_eval": len(ev),
        "K_D_le_3": k,
        "median_D": None if not ds else round(float(np.median(ds)), 2),
        "iqr_D": None if len(ds) < 2 else round(float(np.percentile(ds, 75) -
                                                     np.percentile(ds, 25)), 2),
        "median_fwhm": None if not ev else round(float(np.median(
            [r["fwhm_voxels"] for r in ev])), 2),
        "csv_sha256": hashlib.sha256(csv_bytes).hexdigest(),
        "json_sha256": hashlib.sha256(json_bytes).hexdigest(),
        "records": records,
    }
    with open(RESULT_JSON, "w") as f:
        json.dump(result, f, indent=1)
    with open(f"{OUTDIR}/acquisition_journal.json", "w") as f:
        json.dump({"base": REF_BASE, "chunks": journal}, f, indent=1)
    print(f"n_eval={len(ev)} median_D={result['median_D']} "
          f"K_D_le_3={k}/{len(ev)} iqr_D={result['iqr_D']}")


if __name__ == "__main__":
    main()
