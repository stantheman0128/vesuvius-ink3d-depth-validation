# Phase 5B preregistration: depth-attribution scale-up in the 1.129um footprint

Date frozen: 2026-08-07, Asia/Taipei.
Status: approved continuation. Stan approved the option-3 plan (scale-up, then
packaging) in session on 2026-08-07 after reviewing the Phase 5A O2 outcome.
Execution still requires the code-freeze addendum of section 8 before any body
acquisition.

Inherits every analysis definition, constant, registration chain, and
fail-closed rule of `PHASE5A_PREREGISTRATION_DEPTH_ATTRIBUTION_1129_2026-08-07.md`
section 4, 5, and 6 unchanged, including the composed Phase A offset. Only the
deltas are specified here.

## 1. Purpose

Phase 5A (n = 8 anchors, one cube) ended `O2_DEPTH_SMEARED`: surface
localization passed at threshold, depth thinness failed, with strong per-anchor
heterogeneity. Phase 5B estimates, at roughly 20x the sample, (a) whether
surface localization replicates and (b) the proportions of depth-attribution
quality classes across the 1.129um footprint. Claim boundary identical to 5A:
geometric depth attribution only, never ink identity, never issue completion.

## 2. Cube selection (deterministic, frozen)

1. Candidate cubes: all prediction 256^3 chunks (z,y,x) of the ink-3d zarr.
2. Geometric eligibility: the cube's anchor envelope, local box [30, 225] per
   axis in cube coordinates (anchor eligibility [64, 191] plus the +/-34
   analysis window), mapped through the frozen chain at its 8 corners with the
   8-reference-voxel pad, must lie inside the reference array bounds.
3. From the lexicographically sorted (z,y,x) geometrically eligible list of
   length L, take seed indices floor(i * L / 24) for i = 0..23. For each seed,
   HTTP HEAD the prediction chunk object; on 404 walk forward through the
   sorted list (wrapping at the end) to the next cube not already selected
   until a 200 is found. Every HEAD result is journaled. This lazy-HEAD walk
   is fully deterministic and exists only because HEAD-probing the entire
   eligible list would cost tens of thousands of requests.
4. The result is N <= 24 distinct cubes. No other selection freedom exists.
   P50_INTERIOR (148,58,75) is analyzed again only if selected by this rule;
   its Phase 5A result stays reported either way.

## 3. Anchor rule per cube (frozen, from Phase 4B section 7.1)

Decode the prediction chunk (blosc, manifest settings). Eligible centers have
every coordinate in [64, 191]. Rank voxels by (-prediction value, z, y, x)
ascending (the lexicographic tie policy is frozen here). Greedily retain a
voxel when its Chebyshev distance to every retained anchor is >= 64. Stop at 8.
Cubes yielding fewer than 8 anchors keep what they have; the shortfall is
recorded, and the cube is not replaced. Target sample: up to 192 anchors.

## 4. Statistics

Primary confirmatory gate (the only pass/fail claim):

- **G2s (surface localization):** K = anchors with D <= 3 among evaluable
  anchors, tested by exact one-sided binomial against p0 = 2/7
  (scipy.stats.binomtest, alternative "greater"). Pass iff p <= 0.05.

Everything else is estimation, reported with Wilson 95% intervals:

- **C1** D <= 3 and FWHM <= 12 (localized and thin)
- **C2** D <= 3 and FWHM > 12 (localized, depth-smeared)
- **C3** D > 3 (not localized)
- **C4** non-evaluable (with reasons)

plus per-cube medians of D and FWHM, and the same per-anchor records as 5A.
No thinness pass/fail gate exists in 5B; the pilot already showed
heterogeneity, and 5B quantifies it instead of gating on it.

Evaluability floor: if fewer than 96 evaluable anchors result, the confirmatory
gate is not computed and the run reports
`INCONCLUSIVE_INSUFFICIENT_EVALUABLE_5B` with all records.

## 5. Outcomes

- **O1b LOCALIZATION_REPLICATED:** G2s passes. Deliverable: proportions with
  intervals, feeding the packaging step.
- **O3b LOCALIZATION_NOT_REPLICATED:** G2s fails. Honest negative; Phase 5A's
  threshold-level pass is then treated as fragile and reported as such.
- **O4b INCONCLUSIVE_INSUFFICIENT_EVALUABLE_5B**, **O5b TECHNICAL_FAIL_*:** as
  in 5A.

## 6. Dual-energy chain QA (descriptive side-branch, no scientific gate)

Among selected cubes whose anchor envelope also lies fully inside the 137keV
volume (frozen chain: canonical -> fixed 20230205180739 -> 137keV, both
official transforms, corner test with 8-voxel pad), pick the lexicographically
first. For that cube: acquire its 8 raw-CT 128^3 octant chunks and the mapped
137keV region, then run the Phase A NCC construct (template = raw central
128^3, search pad 16). Report NCC peak, offset, and edge flag. If no selected
cube qualifies, report `DUAL_ENERGY_QA_NOT_AVAILABLE` and do not widen the
selection. This branch informs a future dual-energy preregistration and makes
no depth-attribution claim.

## 7. Budget and fail-closed additions

- Hard self-cap: 2,200 reference chunks total (about 4.6 GB). Exceeding the
  plan aborts with `TECHNICAL_FAIL_BUDGET` before any download.
- Per-cube work is resumable: a cube with a written per-cube result file is
  not re-analyzed; chunk files on disk are reused byte-identically (SHA-256
  journaled on first write, re-hashed on reuse; mismatch aborts).
- Acquisition identity, retry, 404-as-masked, and byte-length rules: 5A
  section 5 unchanged.

## 8. Code commitment

Before the first body acquisition: `phase5b_scaleup_run.py` (importing the
frozen 5A analysis functions unmodified) plus unit tests covering the anchor
rule (synthetic prediction with known ranking and tie cases), the selection
stride, and the Wilson interval, in `tests/test_phase5b_rules.py`. Freeze both
SHA-256 values and the exact command in
`PHASE5B_CODE_FREEZE_ADDENDUM_2026-08-07.md`, then run once. The 5A phantom
suite must still pass unchanged at freeze time.
