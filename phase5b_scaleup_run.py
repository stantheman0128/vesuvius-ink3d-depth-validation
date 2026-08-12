"""Phase 5B scale-up run, per PHASE5B_PREREGISTRATION_SCALEUP_2026-08-07.md.

Imports the frozen Phase 5A analysis functions unmodified. Resumable: per-cube
result files and on-disk chunks are reused; a completed cube is never
re-analyzed.

Usage:
    python phase5b_scaleup_run.py --plan      # cube selection + budget only
    python phase5b_scaleup_run.py             # full run (selection cached)
"""
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve
from scipy.stats import binomtest

from phase5a_depth_attribution_run import (
    CHUNK, CHUNK_BYTES, REF_SHAPE_ZYX, WIN, D_THRESHOLD, FWHM_GATE_MEDIAN,
    analyze_candidate, canonical_to_ref_zyx, load_chain,
)

S3 = "https://vesuvius-challenge-open-data.s3.amazonaws.com/"
PRED_BASE = (S3 + "PHercParis4/representations/predictions/ink-3d/"
             "20260411134726-ink3d-20260428123845-v3-78k-fullsup.zarr/0")
RAW_BASE = S3 + "PHercParis4/volumes/20260411134726-2.400um-0.2m-78keV-masked.zarr/0"
REF_BASE = S3 + "PHercParis4/volumes/20260608103018-1.129um-0.2m-78keV-masked.zarr/0"
E137_BASE = S3 + "PHercParis4/volumes/20260323153942-2.400um-0.2m-137keV-masked.zarr/0"
E137_SHAPE_ZYX = (6625, 8431, 8431)
CANON_SHAPE_ZYX = (75784, 32693, 32693)

N_CUBES = 24
ANCHOR_LO, ANCHOR_HI = 64, 191          # inclusive eligibility per axis
ENVELOPE_LO, ENVELOPE_HI = 30, 225      # anchor box +/- WIN
CHEBYSHEV_MIN = 64
ANCHORS_PER_CUBE = 8
REF_PAD = 8
CHUNK_SELF_CAP = 2200
EVAL_FLOOR = 96
P0 = 2.0 / 7.0
QA_PAD = 16  # NCC search pad, Phase A construct

OUTROOT = "data_phase5b_2026-08-07"
CHUNKDIR = f"{OUTROOT}/ref_chunks"
PREDDIR = f"{OUTROOT}/pred_chunks"
QADIR = f"{OUTROOT}/dual_energy_qa"
CUBEDIR = f"{OUTROOT}/cube_results"
SELECTION_JSON = f"{OUTROOT}/cube_selection.json"
RESULT_JSON = "phase5b_result_2026-08-07.json"

TRANSFORM_137 = "data_crossscan_meta_2026-08-07/transform_20260323153942.json"
TRANSFORM_78 = "data_crossscan_meta_2026-08-07/transform_20260411134726.json"


# ------------------------------------------------------------------ utilities
def http_get(url, expect_bytes=None, retries=3, timeout=90):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                data = r.read()
            if expect_bytes is not None and len(data) != expect_bytes:
                raise RuntimeError(f"size mismatch {url}: {len(data)}")
            return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == retries - 1:
                raise
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(1.5 * (attempt + 1))


def http_head_exists(url):
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60):
            return True
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return False
        raise


def cached_chunk(base, cdir, cz, cy, cx, journal):
    """Fetch/cache one uncompressed 128^3 uint8 chunk; None-body 404 -> zeros."""
    os.makedirs(cdir, exist_ok=True)
    path = f"{cdir}/chunk_{cz}_{cy}_{cx}.bin"
    if os.path.exists(path):
        data = open(path, "rb").read()
        if len(data) not in (0, CHUNK_BYTES):
            raise RuntimeError(f"corrupt cache {path}")
        journal.append({"key": f"{cz}/{cy}/{cx}", "bytes": len(data),
                        "status": "local_reuse",
                        "sha256": hashlib.sha256(data).hexdigest() if data else None})
    else:
        data = http_get(f"{base}/{cz}/{cy}/{cx}", expect_bytes=CHUNK_BYTES)
        if data is None:
            open(path, "wb").close()
            journal.append({"key": f"{cz}/{cy}/{cx}", "bytes": 0,
                            "status": 404, "sha256": None})
            data = b""
        else:
            with open(path, "wb") as f:
                f.write(data)
            journal.append({"key": f"{cz}/{cy}/{cx}", "bytes": len(data),
                            "status": 200,
                            "sha256": hashlib.sha256(data).hexdigest()})
    return (np.frombuffer(data, dtype=np.uint8).reshape(CHUNK, CHUNK, CHUNK)
            if data else np.zeros((CHUNK,) * 3, np.uint8))


def wilson(k, n, z=1.959964):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return [round(float((c - h) / d), 4), round(float((c + h) / d), 4)]


# --------------------------------------------------------- selection geometry
def envelope_corners_global(cube_zyx):
    o = np.array(cube_zyx) * 256
    return np.array([[o[0] + a, o[1] + b, o[2] + c]
                     for a in (ENVELOPE_LO, ENVELOPE_HI)
                     for b in (ENVELOPE_LO, ENVELOPE_HI)
                     for c in (ENVELOPE_LO, ENVELOPE_HI)], float)


def cube_inside_ref(cube_zyx, Ai, t):
    q = canonical_to_ref_zyx(envelope_corners_global(cube_zyx), Ai, t)
    return bool(np.all(q - REF_PAD >= 0)
                and np.all(q + REF_PAD <= np.array(REF_SHAPE_ZYX) - 1))


def plan_selection():
    A, Ai, t = load_chain()
    # canonical bbox of the reference volume: x_can = A @ x_ref + t (xyz)
    ref_corners = np.array([[z, y, x]
                            for z in (0, REF_SHAPE_ZYX[0] - 1)
                            for y in (0, REF_SHAPE_ZYX[1] - 1)
                            for x in (0, REF_SHAPE_ZYX[2] - 1)], float)
    can_xyz = ref_corners[:, ::-1] @ A.T + t
    can_zyx = can_xyz[:, ::-1]
    lo = np.floor(can_zyx.min(axis=0)).astype(int)
    hi = np.ceil(can_zyx.max(axis=0)).astype(int)
    c0 = np.maximum(lo // 256, 0)
    c1 = np.minimum(hi // 256, (np.array(CANON_SHAPE_ZYX) - 1) // 256)
    eligible = []
    for cz in range(c0[0], c1[0] + 1):
        for cy in range(c0[1], c1[1] + 1):
            for cx in range(c0[2], c1[2] + 1):
                if cube_inside_ref((cz, cy, cx), Ai, t):
                    eligible.append((cz, cy, cx))
    eligible.sort()
    L = len(eligible)
    head_journal = []
    selected = []
    taken = set()
    for i in range(N_CUBES):
        j = (i * L) // N_CUBES
        while True:
            cand = eligible[j % L]
            j += 1
            if cand in taken:
                continue
            ok = http_head_exists(f"{PRED_BASE}/{cand[0]}/{cand[1]}/{cand[2]}")
            head_journal.append({"cube": list(cand), "head_200": ok})
            if ok:
                selected.append(cand)
                taken.add(cand)
                break
            taken.add(cand)
            if len(taken) >= L:
                break
    return A, Ai, t, eligible, selected, head_journal


def ref_chunk_plan_for_cube(cube_zyx, anchors_local, Ai, t):
    chunks = set()
    o = np.array(cube_zyx) * 256
    for local in anchors_local:
        c = o + np.array(local)
        corners = np.array([[c[0] + a, c[1] + b, c[2] + d]
                            for a in (-WIN, WIN) for b in (-WIN, WIN)
                            for d in (-WIN, WIN)], float)
        q = canonical_to_ref_zyx(corners, Ai, t)
        qlo = np.maximum(np.floor(q.min(axis=0) - REF_PAD).astype(int), 0)
        qhi = np.minimum(np.ceil(q.max(axis=0) + REF_PAD).astype(int),
                         np.array(REF_SHAPE_ZYX) - 1)
        k0, k1 = qlo // CHUNK, qhi // CHUNK
        for cz in range(k0[0], k1[0] + 1):
            for cy in range(k0[1], k1[1] + 1):
                for cx in range(k0[2], k1[2] + 1):
                    chunks.add((cz, cy, cx))
    return sorted(chunks)


# -------------------------------------------------------------- anchor sampler
def select_anchors(pred_cube):
    """Frozen Phase 4B rule with the frozen lexicographic tie policy."""
    sub = pred_cube[ANCHOR_LO:ANCHOR_HI + 1,
                    ANCHOR_LO:ANCHOR_HI + 1,
                    ANCHOR_LO:ANCHOR_HI + 1]
    n = sub.shape[0]
    flat = sub.reshape(-1)
    zz, yy, xx = np.unravel_index(np.arange(flat.size), sub.shape)
    order = np.lexsort((xx, yy, zz, -flat.astype(np.int16)))
    anchors = []
    for idx in order:
        v = (int(zz[idx]) + ANCHOR_LO, int(yy[idx]) + ANCHOR_LO,
             int(xx[idx]) + ANCHOR_LO)
        if all(max(abs(v[0] - a[0]), abs(v[1] - a[1]), abs(v[2] - a[2]))
               >= CHEBYSHEV_MIN for a in anchors):
            anchors.append(v)
            if len(anchors) == ANCHORS_PER_CUBE:
                break
    return anchors


def decode_pred_bytes(raw):
    import numcodecs
    arr = np.frombuffer(numcodecs.Blosc().decode(raw), dtype=np.uint8)
    if arr.size != 256 ** 3:
        raise RuntimeError(f"prediction decode size {arr.size}")
    return arr.reshape(256, 256, 256)


# ------------------------------------------------------------- dual-energy QA
def chain_can_to_137():
    with open(TRANSFORM_78) as f:
        d78 = json.load(f)
    with open(TRANSFORM_137) as f:
        d137 = json.load(f)
    M78 = np.array(d78["transformation_matrix"], float)
    M137 = np.array(d137["transformation_matrix"], float)
    A78, t78 = M78[:, :3], M78[:, 3]
    A137, t137 = M137[:, :3], M137[:, 3]
    A137i = np.linalg.inv(A137)

    def to137_zyx(pts_zyx):
        p = np.atleast_2d(np.asarray(pts_zyx, float))
        q = ((p[:, ::-1] @ A78.T) + t78 - t137) @ A137i.T
        return q[:, ::-1]
    return to137_zyx


def cube_inside_137(cube_zyx, to137):
    q = to137(envelope_corners_global(cube_zyx))
    return bool(np.all(q - REF_PAD >= 0)
                and np.all(q + REF_PAD <= np.array(E137_SHAPE_ZYX) - 1))


def ncc_search(search, template):
    tz = template.astype(np.float64) - template.mean()
    num = fftconvolve(search.astype(np.float64), tz[::-1, ::-1, ::-1], mode="valid")
    ones = np.ones_like(tz)
    s1 = fftconvolve(search.astype(np.float64), ones, mode="valid")
    s2 = fftconvolve(search.astype(np.float64) ** 2, ones, mode="valid")
    denom = np.sqrt(np.maximum(s2 - s1 ** 2 / tz.size, 1e-9)) * np.linalg.norm(tz)
    ncc = num / denom
    peak = np.unravel_index(np.argmax(ncc), ncc.shape)
    return float(ncc.max()), [int(p) - QA_PAD for p in peak], bool(
        any(p in (0, s - 1) for p, s in zip(peak, ncc.shape)))


def run_dual_energy_qa(selected, Ai, t):
    to137 = chain_can_to_137()
    qa_cubes = [c for c in selected if cube_inside_137(c, to137)]
    if not qa_cubes:
        return {"status": "DUAL_ENERGY_QA_NOT_AVAILABLE"}
    cube = sorted(qa_cubes)[0]
    os.makedirs(QADIR, exist_ok=True)
    journal = []
    raw = np.empty((256, 256, 256), np.uint8)
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                blk = cached_chunk(RAW_BASE, QADIR,
                                   2 * cube[0] + dz, 2 * cube[1] + dy,
                                   2 * cube[2] + dx, journal)
                raw[dz * 128:(dz + 1) * 128, dy * 128:(dy + 1) * 128,
                    dx * 128:(dx + 1) * 128] = blk
    template = raw[64:192, 64:192, 64:192]

    o = np.array(cube) * 256
    n = 128 + 2 * QA_PAD
    lo = o + 64 - QA_PAD
    zz, yy, xx = np.meshgrid(np.arange(n) + lo[0], np.arange(n) + lo[1],
                             np.arange(n) + lo[2], indexing="ij")
    pts = np.stack([zz.ravel(), yy.ravel(), xx.ravel()], axis=1).astype(float)
    q = to137(pts)
    qlo = np.maximum(np.floor(q.min(axis=0) - REF_PAD).astype(int), 0)
    qhi = np.minimum(np.ceil(q.max(axis=0) + REF_PAD).astype(int),
                     np.array(E137_SHAPE_ZYX) - 1)
    k0, k1 = qlo // CHUNK, qhi // CHUNK
    sub = np.zeros(tuple((k1 - k0 + 1) * CHUNK), np.uint8)
    for cz in range(k0[0], k1[0] + 1):
        for cy in range(k0[1], k1[1] + 1):
            for cx in range(k0[2], k1[2] + 1):
                blk = cached_chunk(E137_BASE, QADIR, cz, cy, cx, journal)
                oz, oy, ox = (np.array([cz, cy, cx]) - k0) * CHUNK
                sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk
    search = map_coordinates(sub.astype(np.float32), (q - k0 * CHUNK).T,
                             order=1, mode="constant", cval=0.0).reshape(n, n, n)
    peak, offset, on_edge = ncc_search(search, template)
    result = {"status": "DUAL_ENERGY_QA_DONE", "cube": list(cube),
              "ncc_peak": round(peak, 4), "offset_zyx": offset,
              "peak_on_edge": on_edge,
              "downloaded_mb": sum(j["bytes"] for j in journal) / 1e6}
    with open(f"{QADIR}/qa_journal.json", "w") as f:
        json.dump(journal, f, indent=1)
    return result


# ------------------------------------------------------------------- main run
def main():
    plan_only = "--plan" in sys.argv
    os.makedirs(OUTROOT, exist_ok=True)
    if os.path.exists(SELECTION_JSON):
        sel = json.load(open(SELECTION_JSON))
        A, Ai, t = load_chain()
        eligible = [tuple(c) for c in sel["eligible"]]
        selected = [tuple(c) for c in sel["selected"]]
    else:
        A, Ai, t, eligible, selected, head_journal = plan_selection()
        with open(SELECTION_JSON, "w") as f:
            json.dump({"eligible_count": len(eligible),
                       "eligible": [list(c) for c in eligible],
                       "selected": [list(c) for c in selected],
                       "head_journal": head_journal}, f)
    print(f"eligible cubes: {len(eligible)}; selected: {len(selected)}")
    print(f"selected: {selected}")

    # budget check across all cubes using provisional 8 center-spread anchors
    # (exact plans are computed per cube after anchors are known; the cap is
    #  enforced cumulatively during the run as well)
    if plan_only:
        return

    os.makedirs(CUBEDIR, exist_ok=True)
    total_chunks = 0
    for cube in selected:
        tag = f"{cube[0]}_{cube[1]}_{cube[2]}"
        out = f"{CUBEDIR}/cube_{tag}.json"
        if os.path.exists(out):
            total_chunks += json.load(open(out)).get("ref_chunks", 0)
            print(f"cube {tag}: cached")
            continue
        os.makedirs(PREDDIR, exist_ok=True)
        ppath = f"{PREDDIR}/pred_{tag}.bin"
        if os.path.exists(ppath):
            praw = open(ppath, "rb").read()
        else:
            praw = http_get(f"{PRED_BASE}/{cube[0]}/{cube[1]}/{cube[2]}")
            if praw is None:
                raise SystemExit(f"TECHNICAL_FAIL_ACQUISITION: pred 404 {tag}")
            with open(ppath, "wb") as f:
                f.write(praw)
        pred_cube = decode_pred_bytes(praw)
        anchors = select_anchors(pred_cube)
        chunks = ref_chunk_plan_for_cube(cube, anchors, Ai, t)
        total_chunks += len(chunks)
        if total_chunks > CHUNK_SELF_CAP:
            raise SystemExit("TECHNICAL_FAIL_BUDGET: cumulative chunk cap")
        journal = []
        store = {}
        for (cz, cy, cx) in chunks:
            store[(cz, cy, cx)] = cached_chunk(REF_BASE, CHUNKDIR, cz, cy, cx,
                                               journal)
        arr = np.array(chunks)
        k0 = arr.min(axis=0)
        k1 = arr.max(axis=0)
        sub = np.zeros(tuple((k1 - k0 + 1) * CHUNK), np.uint8)
        for (cz, cy, cx), blk in store.items():
            oz, oy, ox = (np.array([cz, cy, cx]) - k0) * CHUNK
            sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk
        sub_origin = k0 * CHUNK
        acquired = set(chunks)

        o = np.array(cube) * 256
        nwin = 2 * WIN + 1
        off = np.stack(np.meshgrid(*(np.arange(nwin) - WIN,) * 3,
                                   indexing="ij"), axis=-1).reshape(-1, 3)
        records = {}
        for ai, local in enumerate(anchors, 1):
            center = o + np.array(local)
            q = canonical_to_ref_zyx(center + off, Ai, t)
            qc = {tuple(cc) for cc in np.unique(np.floor(q).astype(int) // CHUNK,
                                                axis=0)}
            if qc - acquired:
                raise SystemExit(f"TECHNICAL_FAIL_COVERAGE: {tag} anchor {ai}")
            window = map_coordinates(sub.astype(np.float32),
                                     (q - sub_origin).T, order=1,
                                     mode="constant", cval=0.0
                                     ).reshape(nwin, nwin, nwin)
            records[f"{tag}:{ai}"] = dict(
                analyze_candidate(window, pred_cube, local),
                local_anchor_zyx=list(local))
        with open(out, "w") as f:
            json.dump({"cube": list(cube), "anchors": [list(a) for a in anchors],
                       "anchor_shortfall": ANCHORS_PER_CUBE - len(anchors),
                       "ref_chunks": len(chunks),
                       "downloaded_mb": sum(j["bytes"] for j in journal) / 1e6,
                       "records": records}, f, indent=1)
        print(f"cube {tag}: anchors={len(anchors)} chunks={len(chunks)} "
              f"total_chunks={total_chunks}", flush=True)

    # ---------------------------------------------------------------- aggregate
    all_records = {}
    per_cube = {}
    for cube in selected:
        tag = f"{cube[0]}_{cube[1]}_{cube[2]}"
        d = json.load(open(f"{CUBEDIR}/cube_{tag}.json"))
        all_records.update(d["records"])
        ev = [r for r in d["records"].values() if r["evaluable"]]
        per_cube[tag] = {
            "n_anchors": len(d["records"]), "n_eval": len(ev),
            "median_D": (float(np.median([r["D_voxels"] for r in ev]))
                         if ev else None),
            "median_fwhm": (float(np.median([r["fwhm_voxels"] for r in ev]))
                            if ev else None)}

    ev = [r for r in all_records.values() if r["evaluable"]]
    n_eval = len(ev)
    n_total = len(all_records)
    K = sum(1 for r in ev if r["D_voxels"] <= D_THRESHOLD)
    c1 = sum(1 for r in ev if r["D_voxels"] <= D_THRESHOLD
             and r["fwhm_voxels"] <= FWHM_GATE_MEDIAN)
    c2 = K - c1
    c3 = n_eval - K
    c4 = n_total - n_eval
    if n_eval < EVAL_FLOOR:
        outcome = "O4b_INCONCLUSIVE_INSUFFICIENT_EVALUABLE_5B"
        pval = None
        g2s = None
    else:
        pval = float(binomtest(K, n_eval, P0, alternative="greater").pvalue)
        g2s = bool(pval <= 0.05)
        outcome = ("O1b_LOCALIZATION_REPLICATED" if g2s
                   else "O3b_LOCALIZATION_NOT_REPLICATED")

    qa = run_dual_energy_qa(selected, Ai, t)

    result = {
        "preregistration": "PHASE5B_PREREGISTRATION_SCALEUP_2026-08-07.md",
        "n_cubes": len(selected), "n_anchors": n_total, "n_eval": n_eval,
        "K_D_le_3": K, "binom_p_vs_2_7": pval, "G2s_pass": g2s,
        "outcome": outcome,
        "classes": {
            "C1_localized_thin": {"n": c1, "wilson": wilson(c1, n_eval)},
            "C2_localized_smeared": {"n": c2, "wilson": wilson(c2, n_eval)},
            "C3_not_localized": {"n": c3, "wilson": wilson(c3, n_eval)},
            "C4_non_evaluable": {"n": c4, "wilson": wilson(c4, n_total)},
        },
        "per_cube": per_cube,
        "dual_energy_qa": qa,
        "records": all_records,
    }
    with open(RESULT_JSON, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("records", "per_cube")}, indent=1))


if __name__ == "__main__":
    main()
