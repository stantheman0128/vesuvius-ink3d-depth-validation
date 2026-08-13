# Phase 6 preregistration: depth attribution against the energy-independent 137keV reference

Date frozen: 2026-08-13, Asia/Taipei.
Status: approved continuation (Stan's ordered list item 3, 2026-08-13).
Execution requires the code-freeze addendum of section 7 and must not start
while the Phase 7 run is still executing (resource discipline).

Inherits the analysis definitions of Phase 5A section 6 (smoothing, structure
tensor, profile, surface prominence, D, FWHM, evaluability) and the anchor
rule, statistics shape, and acquisition rules of 5B/5C, with the deltas below.

## 1. Question

Phase 5C established surface localization against the 1.129um reference,
which shares the beam energy (78keV) with the canonical volume. The 137keV
volume is an independent scan at a different energy. This phase asks whether
surface localization replicates when the reference is independent in physics,
not only in scan session and processing.

## 2. Registration chain and per-cube correction (frozen)

- Base chain: canonical x_78 -> fixed 20230205180739 via the official 78keV
  transform, then fixed -> 137keV via the inverse of the official 137keV
  transform. Landmark residuals are larger on this chain (mean 1.9 fixed
  voxels, max 6.4), so a per-cube correction is preregistered:
- For each selected cube, acquire its 8 raw-CT octant chunks (canonical
  volume) and the mapped 137keV region, and run the Phase A NCC construct
  (template = raw central 128^3, search pad 16). If the peak NCC >= 0.85 and
  the peak is not on the search edge, the integer offset is applied to the
  cube's chain for all its anchors. Otherwise the whole cube is
  `CUBE_CHAIN_FAIL` and its anchors are non-evaluable with that reason.
  No other correction, and no per-anchor re-registration.

## 3. Sample (deterministic)

1. Eligible cubes: prediction 256^3 chunks whose anchor envelope (local
   [30, 225] box) maps fully inside the 137keV array with the 8-voxel pad
   through the base chain, before any correction.
2. Excluded: every cube selected by Phase 5B or 5C.
3. Selection: lexicographic sort, seed indices floor(i * L / 12) for
   i = 0..11, lazy-HEAD walk exactly as 5C (prediction chunk 200 plus a
   137keV content prescreen: the chunk containing the mapped cube center has
   nonzero fraction >= 0.5). N <= 12 cubes, up to 96 anchors.
4. Anchors per cube: the frozen 5B rule (top prediction, Chebyshev 64,
   eight anchors, lexicographic ties).

## 4. Measurement and statistics

Per anchor: D and FWHM exactly as in 5A section 6, with the reference window
resampled from the 137keV volume through the corrected chain. Note the
reference here has the same 2.4um voxel pitch as the canonical volume, so
surface localization precision is grid-limited rather than resolution-gained.

- **G2e (energy-independent localization):** K = anchors with D <= 3 among
  evaluable anchors; exact one-sided binomial against p0 = 2/7, alpha 0.05.
- Evaluability floor: 48 evaluable anchors; below it the gate is not
  computed (`INCONCLUSIVE_INSUFFICIENT_EVALUABLE_6`).
- Classes C1 to C4 and per-cube medians reported with Wilson intervals,
  descriptive, as in 5C. Per-cube NCC and offsets reported in full.

## 5. Outcomes

`O1e_LOCALIZATION_CONFIRMED_DUAL_ENERGY` / `O3e_NOT_CONFIRMED` /
`O4e_INCONCLUSIVE_INSUFFICIENT_EVALUABLE_6` / `O5e_TECHNICAL_FAIL_*`.
A pass upgrades the 5C claim from scan-independent to physics-independent
geometric validation; it still makes no ink-identity claim.

## 6. Budget

Per cube about 8 raw chunks (16 MB) + a 137keV region (roughly 27 to 64
chunks, 56 to 134 MB) + one prediction chunk. Twelve cubes fit in about
1.5 to 2 GB. Hard cap: 1,100 new chunks for this phase. 137keV chunks are
cached and journaled in their own directory (they are a different volume
from the shared 1.129um cache).

## 7. Code commitment

`phase6_dual_energy_run.py` (importing frozen 5A/5B functions; the 137keV
chain helper and per-cube correction are new code) plus
`tests/test_phase6_rules.py` covering the chain direction on the official
landmark pairs, the NCC gate thresholds, and the correction application.
Freeze SHA-256 values and the exact command in
`PHASE6_CODE_FREEZE_ADDENDUM_2026-08-13.md`; the 5A phantom suite must still
pass; then run once (resumable).
