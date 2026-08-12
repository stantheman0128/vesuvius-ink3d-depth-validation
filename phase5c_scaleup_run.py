"""Phase 5C run, per PHASE5C_PREREGISTRATION_CONTENT_PRESCREEN_2026-08-07.md.

Identical to the frozen 5B flow except: (1) masked-content prescreen in the
selection walk; (2) the 24 Phase 5B cubes are excluded; (3) no dual-energy QA.
Imports frozen 5A/5B functions unmodified.
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.stats import binomtest

from phase5a_depth_attribution_run import (
    CHUNK, WIN, D_THRESHOLD, FWHM_GATE_MEDIAN, analyze_candidate,
    canonical_to_ref_zyx, load_chain,
)
from phase5b_scaleup_run import (
    ANCHORS_PER_CUBE, CHUNK_SELF_CAP, EVAL_FLOOR, N_CUBES, P0, PRED_BASE,
    REF_BASE, cached_chunk, cube_inside_ref, decode_pred_bytes,
    http_get, http_head_exists, ref_chunk_plan_for_cube, select_anchors,
    wilson,
)
from phase5b_scaleup_run import CHUNKDIR as CHUNKDIR_5B
from phase5b_scaleup_run import SELECTION_JSON as SELECTION_5B

CONTENT_MIN_FRACTION = 0.5

OUTROOT = "data_phase5c_2026-08-07"
CHUNKDIR = CHUNKDIR_5B  # shared journaled cache, per preregistration
PREDDIR = f"{OUTROOT}/pred_chunks"
CUBEDIR = f"{OUTROOT}/cube_results"
SELECTION_JSON = f"{OUTROOT}/cube_selection.json"
RESULT_JSON = "phase5c_result_2026-08-07.json"


def content_ok(chunk_arr):
    """Prescreen predicate: masked-content fraction of the probe chunk."""
    return bool(np.mean(chunk_arr > 0) >= CONTENT_MIN_FRACTION)


def plan_selection_5c(Ai, t, eligible, excluded):
    L = len(eligible)
    journal = []
    selected = []
    taken = set(excluded)
    probe_journal = []
    for i in range(N_CUBES):
        j = (i * L) // N_CUBES
        while True:
            cand = eligible[j % L]
            j += 1
            if cand in taken:
                continue
            taken.add(cand)
            ok_head = http_head_exists(
                f"{PRED_BASE}/{cand[0]}/{cand[1]}/{cand[2]}")
            entry = {"cube": list(cand), "head_200": ok_head}
            if ok_head:
                center = np.array(cand) * 256 + 128
                q = canonical_to_ref_zyx(center.astype(float), Ai, t)[0]
                ck = tuple((np.floor(q).astype(int) // CHUNK).tolist())
                blk = cached_chunk(REF_BASE, CHUNKDIR, *ck, probe_journal)
                frac = float(np.mean(blk > 0))
                entry["probe_chunk"] = list(ck)
                entry["content_fraction"] = round(frac, 4)
                entry["content_ok"] = content_ok(blk)
            journal.append(entry)
            if ok_head and entry.get("content_ok"):
                selected.append(cand)
                break
            if len(taken) >= L:
                break
    return selected, journal, probe_journal


def main():
    plan_only = "--plan" in sys.argv
    os.makedirs(OUTROOT, exist_ok=True)
    A, Ai, t = load_chain()
    sel5b = json.load(open(SELECTION_5B))
    eligible = [tuple(c) for c in sel5b["eligible"]]
    excluded = {tuple(c) for c in sel5b["selected"]}

    if os.path.exists(SELECTION_JSON):
        sel = json.load(open(SELECTION_JSON))
        selected = [tuple(c) for c in sel["selected"]]
    else:
        selected, walk_journal, probe_journal = plan_selection_5c(
            Ai, t, eligible, excluded)
        with open(SELECTION_JSON, "w") as f:
            json.dump({"selected": [list(c) for c in selected],
                       "excluded_5b": [list(c) for c in sorted(excluded)],
                       "walk_journal": walk_journal,
                       "probe_journal": probe_journal}, f, indent=1)
    print(f"selected {len(selected)} cubes: {selected}", flush=True)
    if plan_only:
        return

    os.makedirs(CUBEDIR, exist_ok=True)
    os.makedirs(PREDDIR, exist_ok=True)
    total_new_chunks = 0
    for cube in selected:
        tag = f"{cube[0]}_{cube[1]}_{cube[2]}"
        out = f"{CUBEDIR}/cube_{tag}.json"
        if os.path.exists(out):
            print(f"cube {tag}: cached")
            continue
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
        new = sum(1 for (cz, cy, cx) in chunks if not os.path.exists(
            f"{CHUNKDIR}/chunk_{cz}_{cy}_{cx}.bin"))
        total_new_chunks += new
        if total_new_chunks > CHUNK_SELF_CAP:
            raise SystemExit("TECHNICAL_FAIL_BUDGET: cumulative chunk cap")
        journal = []
        store = {}
        for (cz, cy, cx) in chunks:
            store[(cz, cy, cx)] = cached_chunk(REF_BASE, CHUNKDIR, cz, cy, cx,
                                               journal)
        arr = np.array(chunks)
        k0, k1 = arr.min(axis=0), arr.max(axis=0)
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
            qc = {tuple(cc) for cc in np.unique(
                np.floor(q).astype(int) // CHUNK, axis=0)}
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
            json.dump({"cube": list(cube),
                       "anchors": [list(a) for a in anchors],
                       "anchor_shortfall": ANCHORS_PER_CUBE - len(anchors),
                       "ref_chunks": len(chunks), "new_chunks": new,
                       "records": records}, f, indent=1)
        print(f"cube {tag}: anchors={len(anchors)} chunks={len(chunks)} "
              f"new={new} total_new={total_new_chunks}", flush=True)

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
    n_eval, n_total = len(ev), len(all_records)
    K = sum(1 for r in ev if r["D_voxels"] <= D_THRESHOLD)
    c1 = sum(1 for r in ev if r["D_voxels"] <= D_THRESHOLD
             and r["fwhm_voxels"] <= FWHM_GATE_MEDIAN)
    c2, c3, c4 = K - c1, n_eval - K, n_total - n_eval
    if n_eval < EVAL_FLOOR:
        outcome, pval, g2s = "O4c_INCONCLUSIVE_INSUFFICIENT_EVALUABLE", None, None
    else:
        pval = float(binomtest(K, n_eval, P0, alternative="greater").pvalue)
        g2s = bool(pval <= 0.05)
        outcome = ("O1c_LOCALIZATION_CONFIRMED" if g2s
                   else "O3c_LOCALIZATION_NOT_CONFIRMED")

    result = {
        "preregistration":
            "PHASE5C_PREREGISTRATION_CONTENT_PRESCREEN_2026-08-07.md",
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
        "records": all_records,
    }
    with open(RESULT_JSON, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("records", "per_cube")}, indent=1))


if __name__ == "__main__":
    main()
