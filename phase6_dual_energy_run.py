"""Phase 6 run: depth attribution against the 137keV energy-independent
reference, per PHASE6_PREREGISTRATION_DUAL_ENERGY_2026-08-13.md.

Imports frozen 5A/5B functions unmodified. New code: the 137keV chain, the
per-cube NCC correction gate, and the 137keV-specific selection.

Usage:
    python phase6_dual_energy_run.py --plan
    python phase6_dual_energy_run.py
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.stats import binomtest

from phase5a_depth_attribution_run import (
    CHUNK, WIN, D_THRESHOLD, FWHM_GATE_MEDIAN, analyze_candidate,
)
from phase5b_scaleup_run import (
    ANCHORS_PER_CUBE, PRED_BASE, RAW_BASE, E137_BASE, E137_SHAPE_ZYX,
    QA_PAD, cached_chunk, chain_can_to_137, decode_pred_bytes,
    envelope_corners_global, http_get, http_head_exists, ncc_search,
    select_anchors, wilson,
)

N_CUBES = 12
REF_PAD = 8
NCC_MIN = 0.85
NEW_CHUNK_CAP = 1100
EVAL_FLOOR = 48
P0 = 2.0 / 7.0
CONTENT_MIN_FRACTION = 0.5
G2_TABLE_ALPHA = 0.05

OUTROOT = "data_phase6_2026-08-13"
CHUNKDIR = f"{OUTROOT}/e137_chunks"
RAWDIR = f"{OUTROOT}/raw_chunks"
PREDDIR = f"{OUTROOT}/pred_chunks"
CUBEDIR = f"{OUTROOT}/cube_results"
SELECTION_JSON = f"{OUTROOT}/cube_selection.json"
RESULT_JSON = "phase6_result_2026-08-13.json"
SEEN_5B = "data_phase5b_2026-08-07/cube_selection.json"
SEEN_5C = "data_phase5c_2026-08-07/cube_selection.json"


def to137():
    return chain_can_to_137()


def cube_inside_137(cube, fn, pad=REF_PAD):
    q = fn(envelope_corners_global(cube))
    return bool(np.all(q - pad >= 0)
                and np.all(q + pad <= np.array(E137_SHAPE_ZYX) - 1))


def apply_offset(fn, off_zyx):
    """Corrected chain fn'(v) = fn(v + offset), the 5A input-shift convention:
    the NCC peak says canonical content at v is found at fn(v + offset)."""
    o = np.asarray(off_zyx, float)

    def g(pts):
        return fn(np.asarray(pts, float) + o)
    return g


def cube_ncc_correction(cube, fn, journal):
    """Phase A construct: raw central 128^3 vs 137keV search window."""
    raw = np.empty((256, 256, 256), np.uint8)
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                blk = cached_chunk(RAW_BASE, RAWDIR, 2 * cube[0] + dz,
                                   2 * cube[1] + dy, 2 * cube[2] + dx, journal)
                raw[dz * 128:(dz + 1) * 128, dy * 128:(dy + 1) * 128,
                    dx * 128:(dx + 1) * 128] = blk
    template = raw[64:192, 64:192, 64:192]
    o = np.array(cube) * 256
    n = 128 + 2 * QA_PAD
    lo = o + 64 - QA_PAD
    zz, yy, xx = np.meshgrid(np.arange(n) + lo[0], np.arange(n) + lo[1],
                             np.arange(n) + lo[2], indexing="ij")
    pts = np.stack([zz.ravel(), yy.ravel(), xx.ravel()], axis=1).astype(float)
    q = fn(pts)
    qlo = np.maximum(np.floor(q.min(axis=0) - REF_PAD).astype(int), 0)
    qhi = np.minimum(np.ceil(q.max(axis=0) + REF_PAD).astype(int),
                     np.array(E137_SHAPE_ZYX) - 1)
    k0, k1 = qlo // CHUNK, qhi // CHUNK
    sub = np.zeros(tuple((k1 - k0 + 1) * CHUNK), np.uint8)
    for cz in range(k0[0], k1[0] + 1):
        for cy in range(k0[1], k1[1] + 1):
            for cx in range(k0[2], k1[2] + 1):
                blk = cached_chunk(E137_BASE, CHUNKDIR, cz, cy, cx, journal)
                oz, oy, ox = (np.array([cz, cy, cx]) - k0) * CHUNK
                sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk
    search = map_coordinates(sub.astype(np.float32), (q - k0 * CHUNK).T,
                             order=1, mode="constant", cval=0.0).reshape(n, n, n)
    peak, offset, on_edge = ncc_search(search, template)
    ok = bool(peak >= NCC_MIN and not on_edge)
    # search grid index i corresponds to canonical lo+i; template sits at
    # canonical o+64. Peak at index p means content of canonical v is found
    # at fn(v + offset) with offset = p - QA_PAD, i.e. the corrected chain is
    # fn'(v) = fn(v + offset).
    return {"ncc": round(peak, 4), "offset_zyx": offset, "on_edge": on_edge,
            "pass": ok}


def plan_selection(fn):
    lim = np.array(E137_SHAPE_ZYX)
    # canonical bbox of the 137keV volume via inverse mapping of its corners:
    # walk the known overlap band instead (chunk ranges measured in the
    # crossscan analysis: z 164..191, y 36..71, x 36..70) and corner-test.
    eligible = []
    for cz in range(164, 192):
        for cy in range(36, 72):
            for cx in range(36, 71):
                if cube_inside_137((cz, cy, cx), fn):
                    eligible.append((cz, cy, cx))
    eligible.sort()
    seen = set()
    for p in (SEEN_5B, SEEN_5C):
        if os.path.exists(p):
            seen.update(tuple(c) for c in json.load(open(p))["selected"])
    L = len(eligible)
    selected, taken, walk = [], set(seen), []
    probe_journal = []
    for i in range(N_CUBES):
        j = (i * L) // N_CUBES
        while True:
            cand = eligible[j % L]
            j += 1
            if cand in taken:
                continue
            taken.add(cand)
            ok = http_head_exists(f"{PRED_BASE}/{cand[0]}/{cand[1]}/{cand[2]}")
            e = {"cube": list(cand), "head_200": ok}
            if ok:
                center = np.array(cand) * 256 + 128
                q = fn(center[None].astype(float))[0]
                ck = tuple((np.floor(q).astype(int) // CHUNK).tolist())
                blk = cached_chunk(E137_BASE, CHUNKDIR, *ck, probe_journal)
                frac = float(np.mean(blk > 0))
                e["content_fraction"] = round(frac, 4)
                e["content_ok"] = bool(frac >= CONTENT_MIN_FRACTION)
            walk.append(e)
            if ok and e.get("content_ok"):
                selected.append(cand)
                break
            if len(taken) >= L + len(seen):
                break
    return eligible, selected, walk, probe_journal


def ref_chunks_for_anchors(cube, anchors, fn):
    chunks = set()
    o = np.array(cube) * 256
    for local in anchors:
        c = o + np.array(local)
        corners = np.array([[c[0] + a, c[1] + b, c[2] + d]
                            for a in (-WIN, WIN) for b in (-WIN, WIN)
                            for d in (-WIN, WIN)], float)
        q = fn(corners)
        lo = np.maximum(np.floor(q.min(axis=0) - REF_PAD).astype(int), 0)
        hi = np.minimum(np.ceil(q.max(axis=0) + REF_PAD).astype(int),
                        np.array(E137_SHAPE_ZYX) - 1)
        k0, k1 = lo // CHUNK, hi // CHUNK
        for cz in range(k0[0], k1[0] + 1):
            for cy in range(k0[1], k1[1] + 1):
                for cx in range(k0[2], k1[2] + 1):
                    chunks.add((cz, cy, cx))
    return sorted(chunks)


def main():
    plan_only = "--plan" in sys.argv
    os.makedirs(OUTROOT, exist_ok=True)
    base_fn = to137()

    if os.path.exists(SELECTION_JSON):
        sel = json.load(open(SELECTION_JSON))
        selected = [tuple(c) for c in sel["selected"]]
    else:
        eligible, selected, walk, probe_journal = plan_selection(base_fn)
        json.dump({"eligible_count": len(eligible),
                   "selected": [list(c) for c in selected],
                   "walk_journal": walk, "probe_journal": probe_journal},
                  open(SELECTION_JSON, "w"))
    print(f"selected {len(selected)} cubes: {selected}")
    if plan_only:
        return

    os.makedirs(CUBEDIR, exist_ok=True)
    os.makedirs(PREDDIR, exist_ok=True)
    total_new = 0
    for cube in selected:
        tag = f"{cube[0]}_{cube[1]}_{cube[2]}"
        out = f"{CUBEDIR}/cube_{tag}.json"
        if os.path.exists(out):
            continue
        journal = []
        corr = cube_ncc_correction(cube, base_fn, journal)
        rec = {"cube": list(cube), "ncc_correction": corr}
        if not corr["pass"]:
            rec["records"] = {}
            rec["cube_chain_fail"] = True
            json.dump(rec, open(out, "w"), indent=1)
            print(f"cube {tag}: CUBE_CHAIN_FAIL ncc={corr['ncc']}")
            continue
        fn = apply_offset(base_fn, corr["offset_zyx"])

        ppath = f"{PREDDIR}/pred_{tag}.bin"
        if os.path.exists(ppath):
            praw = open(ppath, "rb").read()
        else:
            praw = http_get(f"{PRED_BASE}/{cube[0]}/{cube[1]}/{cube[2]}")
            if praw is None:
                raise SystemExit(f"TECHNICAL_FAIL_ACQUISITION: pred 404 {tag}")
            open(ppath, "wb").write(praw)
        pred_cube = decode_pred_bytes(praw)
        anchors = select_anchors(pred_cube)
        chunks = ref_chunks_for_anchors(cube, anchors, fn)
        new = sum(1 for (cz, cy, cx) in chunks if not os.path.exists(
            f"{CHUNKDIR}/chunk_{cz}_{cy}_{cx}.bin"))
        total_new += new
        if total_new > NEW_CHUNK_CAP:
            raise SystemExit("TECHNICAL_FAIL_BUDGET")
        store = {}
        for k in chunks:
            store[k] = cached_chunk(E137_BASE, CHUNKDIR, *k, journal)
        arr = np.array(chunks)
        k0, k1 = arr.min(axis=0), arr.max(axis=0)
        sub = np.zeros(tuple((k1 - k0 + 1) * CHUNK), np.uint8)
        for k, blk in store.items():
            oz, oy, ox = (np.array(k) - k0) * CHUNK
            sub[oz:oz + CHUNK, oy:oy + CHUNK, ox:ox + CHUNK] = blk

        o = np.array(cube) * 256
        nwin = 2 * WIN + 1
        off = np.stack(np.meshgrid(*(np.arange(nwin) - WIN,) * 3,
                                   indexing="ij"), axis=-1).reshape(-1, 3)
        acquired = set(chunks)
        records = {}
        for ai, local in enumerate(anchors, 1):
            center = o + np.array(local)
            q = fn(center + off)
            qc = {tuple(cc) for cc in np.unique(
                np.floor(q).astype(int) // CHUNK, axis=0)}
            if qc - acquired:
                raise SystemExit(f"TECHNICAL_FAIL_COVERAGE: {tag}:{ai}")
            window = map_coordinates(sub.astype(np.float32),
                                     (q - k0 * CHUNK).T, order=1,
                                     mode="constant", cval=0.0
                                     ).reshape(nwin, nwin, nwin)
            records[f"{tag}:{ai}"] = dict(
                analyze_candidate(window, pred_cube, local),
                local_anchor_zyx=list(local))
        rec["records"] = records
        rec["ref_chunks"] = len(chunks)
        json.dump(rec, open(out, "w"), indent=1)
        json.dump(journal, open(f"{CUBEDIR}/journal_{tag}.json", "w"))
        print(f"cube {tag}: ncc={corr['ncc']} off={corr['offset_zyx']} "
              f"anchors={len(anchors)} chunks={len(chunks)}", flush=True)

    all_records = {}
    chain_fails = 0
    ncc_summary = {}
    for cube in selected:
        tag = f"{cube[0]}_{cube[1]}_{cube[2]}"
        d = json.load(open(f"{CUBEDIR}/cube_{tag}.json"))
        ncc_summary[tag] = d["ncc_correction"]
        if d.get("cube_chain_fail"):
            chain_fails += 1
        all_records.update(d.get("records", {}))

    ev = [r for r in all_records.values() if r["evaluable"]]
    n_eval, n_total = len(ev), len(all_records)
    K = sum(1 for r in ev if r["D_voxels"] <= D_THRESHOLD)
    c1 = sum(1 for r in ev if r["D_voxels"] <= D_THRESHOLD
             and r["fwhm_voxels"] <= FWHM_GATE_MEDIAN)
    if n_eval < EVAL_FLOOR:
        outcome, pval, g2e = "O4e_INCONCLUSIVE_INSUFFICIENT_EVALUABLE_6", None, None
    else:
        pval = float(binomtest(K, n_eval, P0, alternative="greater").pvalue)
        g2e = bool(pval <= G2_TABLE_ALPHA)
        outcome = ("O1e_LOCALIZATION_CONFIRMED_DUAL_ENERGY" if g2e
                   else "O3e_NOT_CONFIRMED")

    result = {
        "preregistration": "PHASE6_PREREGISTRATION_DUAL_ENERGY_2026-08-13.md",
        "n_cubes": len(selected), "cube_chain_fails": chain_fails,
        "n_anchors": n_total, "n_eval": n_eval, "K_D_le_3": K,
        "binom_p_vs_2_7": pval, "G2e_pass": g2e, "outcome": outcome,
        "classes": {
            "C1_localized_thin": {"n": c1, "wilson": wilson(c1, n_eval)},
            "C2_localized_smeared": {"n": K - c1, "wilson": wilson(K - c1, n_eval)},
            "C3_not_localized": {"n": n_eval - K, "wilson": wilson(n_eval - K, n_eval)},
            "C4_non_evaluable": {"n": n_total - n_eval, "wilson": wilson(n_total - n_eval, n_total)},
        },
        "ncc_per_cube": ncc_summary,
        "records": all_records,
    }
    json.dump(result, open(RESULT_JSON, "w"), indent=1)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("records", "ncc_per_cube")}, indent=1))


if __name__ == "__main__":
    main()
