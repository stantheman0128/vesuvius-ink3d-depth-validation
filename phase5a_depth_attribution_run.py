"""Phase 5A depth-attribution pilot, per the frozen preregistration
PHASE5A_PREREGISTRATION_DEPTH_ATTRIBUTION_1129_2026-08-07.md.

Single deterministic run. Fail-closed: any acquisition anomaly, size mismatch,
or out-of-coverage sample aborts before a scientific statement is emitted.

Usage:
    python phase5a_depth_attribution_run.py            # acquire + analyze
    python phase5a_depth_attribution_run.py --dry-run  # print chunk plan only
"""
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import find_peaks

# ---------------------------------------------------------------- frozen data
REF_BASE = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/"
            "PHercParis4/volumes/20260608103018-1.129um-0.2m-78keV-masked.zarr/0")
REF_SHAPE_ZYX = (59969, 36006, 32354)
CHUNK = 128
CHUNK_BYTES = CHUNK ** 3
TRANSFORM_JSON = "data_crossscan_meta_2026-08-07/transform_20260608103018.json"
PHASEA_OFFSET_ZYX = np.array([0.0, -1.0, 0.0])  # frozen Phase A NCC correction

P50_ORIGIN_ZYX = np.array([37888, 14848, 19200])
ANCHORS_LOCAL_ZYX = {
    "P50:1": (118, 136, 173), "P50:2": (136, 186, 82),
    "P50:3": (182, 156, 172), "P50:4": (129, 64, 144),
    "P50:5": (102, 74, 80), "P50:6": (64, 64, 144),
    "P50:7": (64, 138, 64), "P50:8": (166, 64, 65),
}
PRED_BODY = "data_phase4b_snapshot_2026-08-06/bodies/09_prediction_148_58_75.bin"

WIN = 34            # canonical half-window per candidate
PROFILE_HALF = 24   # profile s = -24..24
REF_PAD = 8         # reference-voxel padding of each mapped bbox
SIGMA_SMOOTH = 2.0
SIGMA_GRAD = 2.0
SIGMA_TENSOR = 4.0
ANISO_RATIO = 2.0
PROMINENCE_FRAC = 0.25
D_THRESHOLD = 3.0
FWHM_CENSOR = 48.0
FWHM_GATE_MEDIAN = 12.0
ZERO_FRAC_CAP = 0.10
CHUNK_SELF_CAP = 500
G2_PRIMARY = {4: 4, 5: 4, 6: 5, 7: 5, 8: 5}       # p0 = 2/7, alpha 0.05
G2_SECONDARY = {4: 3, 5: 3, 6: 3, 7: 4, 8: 4}     # p0 = 7/49, descriptive
OUTDIR = "data_phase5a_1129_chunks_2026-08-07"
RESULT_JSON = "phase5a_result_2026-08-07.json"


# ------------------------------------------------------------------- geometry
def load_chain():
    with open(TRANSFORM_JSON) as f:
        d = json.load(f)
    M = np.array(d["transformation_matrix"], dtype=float)
    A, t = M[:, :3], M[:, 3]
    Ai = np.linalg.inv(A)
    return A, Ai, t


def canonical_to_ref_zyx(pts_zyx, Ai, t):
    """Frozen chain: canonical (z,y,x) voxels -> reference (z,y,x) voxels."""
    p = np.atleast_2d(np.asarray(pts_zyx, dtype=float)) + PHASEA_OFFSET_ZYX
    q_xyz = (p[:, ::-1] - t) @ Ai.T
    return q_xyz[:, ::-1]


# ---------------------------------------------------------------- acquisition
def chunk_plan(Ai, t):
    chunks = set()
    boxes = {}
    for name, local in ANCHORS_LOCAL_ZYX.items():
        c = P50_ORIGIN_ZYX + np.array(local)
        corners = np.array([[c[0] + dz, c[1] + dy, c[2] + dx]
                            for dz in (-WIN, WIN) for dy in (-WIN, WIN)
                            for dx in (-WIN, WIN)], float)
        q = canonical_to_ref_zyx(corners, Ai, t)
        lo = np.floor(q.min(axis=0) - REF_PAD).astype(int)
        hi = np.ceil(q.max(axis=0) + REF_PAD).astype(int)
        lo = np.maximum(lo, 0)
        hi = np.minimum(hi, np.array(REF_SHAPE_ZYX) - 1)
        boxes[name] = (lo.tolist(), hi.tolist())
        c0, c1 = lo // CHUNK, hi // CHUNK
        for cz in range(c0[0], c1[0] + 1):
            for cy in range(c0[1], c1[1] + 1):
                for cx in range(c0[2], c1[2] + 1):
                    chunks.add((cz, cy, cx))
    return sorted(chunks), boxes


def fetch_chunk(cz, cy, cx, journal):
    path = f"{OUTDIR}/chunk_{cz}_{cy}_{cx}.bin"
    url = f"{REF_BASE}/{cz}/{cy}/{cx}"
    if os.path.exists(path):
        data = open(path, "rb").read()
        if len(data) not in (0, CHUNK_BYTES):
            raise RuntimeError(f"local chunk corrupt: {path}")
        journal.append({"key": f"0/{cz}/{cy}/{cx}", "bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest() if data else None,
                        "status": "local_reuse"})
        return (np.frombuffer(data, dtype=np.uint8).reshape(CHUNK, CHUNK, CHUNK)
                if data else np.zeros((CHUNK,) * 3, np.uint8))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
            if len(data) != CHUNK_BYTES:
                raise RuntimeError(f"size mismatch {url}: {len(data)}")
            with open(path, "wb") as f:
                f.write(data)
            journal.append({"key": f"0/{cz}/{cy}/{cx}", "bytes": len(data),
                            "sha256": hashlib.sha256(data).hexdigest(),
                            "status": 200})
            return np.frombuffer(data, dtype=np.uint8).reshape(CHUNK, CHUNK, CHUNK)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                open(path, "wb").close()  # zero-length marker for masked chunk
                journal.append({"key": f"0/{cz}/{cy}/{cx}", "bytes": 0,
                                "sha256": None, "status": 404})
                return np.zeros((CHUNK,) * 3, np.uint8)
            if attempt == 2:
                raise
        except Exception:
            if attempt == 2:
                raise
        time.sleep(1.5 * (attempt + 1))


# ------------------------------------------------------------------- analysis
def surface_D_from_profile(ref_profile):
    """Distances of qualifying |gradient| peaks from s=0; None if no peak."""
    g = np.abs(np.gradient(np.asarray(ref_profile, float)))
    if g.max() <= 0:
        return None
    peaks, _ = find_peaks(g, prominence=PROMINENCE_FRAC * g.max())
    if len(peaks) == 0:
        return None
    s = peaks - PROFILE_HALF
    return float(np.min(np.abs(s)))


def fwhm_from_profile(prof):
    """Linear-interpolated FWHM; FWHM_CENSOR if either side never crosses."""
    p = np.asarray(prof, float)
    ipk = int(np.argmax(p))
    half = (p.max() + p.min()) / 2.0
    left = right = None
    for i in range(ipk, 0, -1):
        if p[i - 1] < half <= p[i]:
            left = (i - 1) + (half - p[i - 1]) / (p[i] - p[i - 1])
            break
    for i in range(ipk, len(p) - 1):
        if p[i + 1] < half <= p[i]:
            right = i + (p[i] - half) / (p[i] - p[i + 1])
            break
    if left is None or right is None:
        return FWHM_CENSOR, True
    return float(right - left), False


def analyze_candidate(window_ref, pred_cube, local_center):
    """window_ref: float32 (69,69,69) resampled reference, anchor at center."""
    rec = {}
    zero_frac = float(np.mean(window_ref < 0.5))
    rec["masked_zero_fraction"] = round(zero_frac, 4)
    R_s = gaussian_filter(window_ref.astype(np.float32), SIGMA_SMOOTH)

    g = [gaussian_filter(R_s, SIGMA_GRAD, order=tuple(int(a == i) for a in range(3)))
         for i in range(3)]
    c = WIN
    J = np.empty((3, 3))
    for i in range(3):
        for j in range(i, 3):
            J[i, j] = J[j, i] = gaussian_filter(g[i] * g[j], SIGMA_TENSOR)[c, c, c]
    w, v = np.linalg.eigh(J)
    lam = w[::-1]
    normal = v[:, ::-1][:, 0]  # (z,y,x) unit vector
    rec["lambda1"] = float(lam[0])
    rec["lambda2"] = float(lam[1])
    aniso_ok = bool(lam[0] >= ANISO_RATIO * lam[1] and lam[0] > 0)
    rec["anisotropy_ok"] = aniso_ok
    rec["normal_zyx"] = [round(float(x), 4) for x in normal]

    s = np.arange(-PROFILE_HALF, PROFILE_HALF + 1, dtype=float)
    pts = c + s[:, None] * normal[None, :]
    ref_prof = map_coordinates(R_s, pts.T, order=1, mode="constant", cval=0.0)
    ppts = np.array(local_center, float) + s[:, None] * normal[None, :]
    pred_prof = map_coordinates(pred_cube.astype(np.float32), ppts.T, order=1,
                                mode="constant", cval=0.0)
    rec["ref_profile"] = [round(float(x), 2) for x in ref_prof]
    rec["pred_profile"] = [round(float(x), 2) for x in pred_prof]

    D = surface_D_from_profile(ref_prof)
    rec["D_voxels"] = D
    fwhm, censored = fwhm_from_profile(pred_prof)
    rec["fwhm_voxels"] = round(fwhm, 2)
    rec["fwhm_censored"] = censored

    rec["evaluable"] = bool(zero_frac <= ZERO_FRAC_CAP and aniso_ok and D is not None)
    reasons = []
    if zero_frac > ZERO_FRAC_CAP:
        reasons.append("masked_zero_fraction_exceeds_cap")
    if not aniso_ok:
        reasons.append("anisotropy_fail")
    if D is None:
        reasons.append("no_qualifying_surface")
    rec["non_evaluable_reasons"] = reasons
    return rec


def decode_prediction():
    import numcodecs
    raw = open(PRED_BODY, "rb").read()
    dec = numcodecs.Blosc().decode(raw)
    arr = np.frombuffer(dec, dtype=np.uint8)
    if arr.size != 256 ** 3:
        raise RuntimeError(f"prediction decode size {arr.size}")
    return arr.reshape(256, 256, 256)


def main():
    dry = "--dry-run" in sys.argv
    A, Ai, t = load_chain()
    chunks, boxes = chunk_plan(Ai, t)
    est_mb = len(chunks) * CHUNK_BYTES / 1e6
    print(f"chunk plan: {len(chunks)} chunks (~{est_mb:.0f} MB)")
    if len(chunks) > CHUNK_SELF_CAP:
        raise SystemExit("TECHNICAL_FAIL_BUDGET: chunk plan exceeds self-cap")
    if dry:
        print(json.dumps(boxes, indent=1))
        return

    os.makedirs(OUTDIR, exist_ok=True)
    journal = []
    store = {}
    for n, (cz, cy, cx) in enumerate(chunks, 1):
        store[(cz, cy, cx)] = fetch_chunk(cz, cy, cx, journal)
        if n % 25 == 0 or n == len(chunks):
            print(f"  fetched {n}/{len(chunks)}", flush=True)
    with open(f"{OUTDIR}/acquisition_journal.json", "w") as f:
        json.dump({"base": REF_BASE, "chunks": journal}, f, indent=1)

    all_c = np.array(chunks)
    c_lo, c_hi = all_c.min(axis=0), all_c.max(axis=0)
    sub = np.zeros(tuple((c_hi - c_lo + 1) * CHUNK), dtype=np.uint8)
    acquired = set(chunks)
    for (cz, cy, cx), blk in store.items():
        oz, oy, ox = (np.array([cz, cy, cx]) - c_lo) * CHUNK
        sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk
    sub_origin = c_lo * CHUNK

    pred_cube = decode_prediction()
    records = {}
    n = 2 * WIN + 1
    off = np.stack(np.meshgrid(*(np.arange(n) - WIN,) * 3, indexing="ij"),
                   axis=-1).reshape(-1, 3)
    for name, local in ANCHORS_LOCAL_ZYX.items():
        center = P50_ORIGIN_ZYX + np.array(local)
        q = canonical_to_ref_zyx(center + off, Ai, t)
        qc = (np.floor(q).astype(int)) // CHUNK
        missing = {tuple(cc) for cc in np.unique(qc, axis=0)} - acquired
        if missing:
            raise SystemExit(f"TECHNICAL_FAIL_COVERAGE: {name} needs {missing}")
        window = map_coordinates(sub.astype(np.float32), (q - sub_origin).T,
                                 order=1, mode="constant", cval=0.0
                                 ).reshape(n, n, n)
        records[name] = analyze_candidate(window, pred_cube, local)
        print(f"  {name}: evaluable={records[name]['evaluable']} "
              f"D={records[name]['D_voxels']} fwhm={records[name]['fwhm_voxels']}")

    ev = [r for r in records.values() if r["evaluable"]]
    n_eval = len(ev)
    K = sum(1 for r in ev if r["D_voxels"] is not None and r["D_voxels"] <= D_THRESHOLD)
    fwhms = [r["fwhm_voxels"] for r in ev]
    med_fwhm = float(np.median(fwhms)) if fwhms else None
    g1 = n_eval >= 4
    g2 = bool(g1 and K >= G2_PRIMARY[n_eval])
    g2_desc = bool(g1 and K >= G2_SECONDARY[n_eval])
    g3 = bool(g1 and med_fwhm is not None and med_fwhm <= FWHM_GATE_MEDIAN)
    if not g1:
        outcome = "O4_INCONCLUSIVE_INSUFFICIENT_EVALUABLE"
    elif g2 and g3:
        outcome = "O1_GO_DEPTH_ATTRIBUTION_SUPPORTED"
    elif g2:
        outcome = "O2_DEPTH_SMEARED"
    else:
        outcome = "O3_NO_SURFACE_PREFERENCE"

    result = {
        "preregistration": "PHASE5A_PREREGISTRATION_DEPTH_ATTRIBUTION_1129_2026-08-07.md",
        "n_eval": n_eval, "K_D_le_3": K, "median_fwhm_voxels": med_fwhm,
        "G1_evaluability": g1, "G2_primary_pass": g2,
        "G2_secondary_descriptive_pass": g2_desc, "G3_thinness": g3,
        "outcome": outcome,
        "chunks_fetched": len(chunks),
        "downloaded_mb": sum(j["bytes"] for j in journal) / 1e6,
        "records": records,
    }
    with open(RESULT_JSON, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=1))


if __name__ == "__main__":
    main()
