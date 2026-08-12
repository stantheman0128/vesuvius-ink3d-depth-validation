# Phase 5C preregistration: scale-up rerun with masked-content prescreen

Date frozen: 2026-08-07, Asia/Taipei.
Status: approved continuation of the option-3 mandate (scale-up, then package).

Phase 5B ended `O4b_INCONCLUSIVE_INSUFFICIENT_EVALUABLE_5B`: 89 evaluable
anchors against the frozen floor of 96. The dominant failure mode is now
diagnosed and journaled: 9 of 24 cubes lie entirely in masked air (75 anchors
with masked-zero fraction near 1.0), because geometric eligibility tested
array bounds, not mask content. Per the 5B preregistration, its confirmatory
gate was not computed and is not claimed. All 5B records remain reported as
descriptive.

Phase 5C repeats the 5B design exactly, with two frozen deltas and nothing
else changed (anchors, analysis, statistics, floors, caps, outcomes all
inherit from 5B, which inherits from 5A):

## Delta 1: masked-content prescreen in the selection walk

During the deterministic stride walk, a cube is selected only if all of:

1. HEAD of its prediction chunk returns 200 (as in 5B);
2. the single reference chunk containing the mapped cube center
   (canonical center = cube * 256 + 128, through the frozen chain) has a
   nonzero-voxel fraction >= 0.5. The probe chunk is acquired through the
   shared journaled cache and every probe result is journaled.

Otherwise the walk advances to the next sorted cube, exactly as the 5B
404-fallback walk does.

## Delta 2: exclusion of seen cubes

The 24 Phase 5B cubes are excluded from selection. Their data has been seen;
reusing any of them would contaminate the confirmatory statistic. The 5C
sample is therefore disjoint from 5B by construction.

## Unchanged and restated

- Target N = 24 cubes, up to 192 anchors; evaluability floor 96.
- Confirmatory statistic: G2s only (exact one-sided binomial, K of D <= 3
  among evaluable, p0 = 2/7, alpha 0.05).
- Classes C1 to C4 with Wilson intervals; per-anchor records in full.
- Chunk self-cap 2,200 for 5C's own acquisitions (shared cache reuse of
  already-journaled chunks does not count against it).
- Outcomes: O1c LOCALIZATION_CONFIRMED / O3c NOT_CONFIRMED /
  O4c INCONCLUSIVE_INSUFFICIENT_EVALUABLE / O5c TECHNICAL_FAIL_*.
- The dual-energy QA branch ran in 5B (NCC 0.924, offset (1,-5,2), off-edge)
  and is not repeated.
- Combined 5B+5C proportions may be reported only as descriptive, clearly
  labeled; the confirmatory claim, if any, rests on 5C alone.

## Code commitment

`phase5c_scaleup_run.py` (importing frozen 5A/5B functions unmodified) plus a
unit test of the prescreen predicate in `tests/test_phase5c_rules.py`. Freeze
SHA-256 values and the exact command in
`PHASE5C_CODE_FREEZE_ADDENDUM_2026-08-07.md`; the 5A and 5B suites must still
pass at freeze; then run once (resumable, as frozen in 5B).
