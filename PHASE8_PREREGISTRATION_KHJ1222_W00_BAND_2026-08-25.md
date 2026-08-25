# Phase 8 preregistration: D/FWHM of khj1222's w00 measured band

Date frozen: 2026-08-25, Asia/Taipei.
Status: continuation of the #192 offer (comment 2026-08-23). Not a confirmatory
hypothesis test with a go/no-go gate. Produces per-anchor records.

Inherits the frozen Phase 5A registration chain (official 1.129um transform
plus the Phase A one-voxel NCC offset). Does not inherit 5B cube selection.

## 1. Purpose and claim boundary

Score the `v4` measured band khj1222 exported in
`submission/depth_anchors/` against the independently scanned 1.129um volume.

This is distance to an independently observed surface. It is not an
ink-identity test. Spatially coherent D versus scatter is the discriminator
between "the per-pixel wander is estimator noise" and "the geometry is right
and still loses in training".

FWHM here is of the independent CT |gradient| peak along the same line, not
of the ink-3d prediction used in Phase 5A. Their declared `half_width_layers`
is reported beside it, not mixed into D.

## 2. Inputs (fail-closed)

1. Sidecar + CSV from
   `https://github.com/khj1222/vesuvius-challenge/tree/main/submission/depth_anchors`
   (downloaded at run time; SHA-256 journaled).
2. Keep sidecar conventions: `layer_step_voxels = 1.0`; do not flip
   `normal_orientation`.
3. Strict anchor set: `sampled_from == center` only (sidecar
   `readback_caveat`).
4. Registration chain: exactly Phase 5A.

## 3. Point construction

Band centre in canonical (x, y, z) voxels:

`P = (base_x, base_y, base_z) + offset_from_plane * (normal_x, normal_y, normal_z)`

Normal is used as exported (unit-normalized, sign not flipped). Profile
samples `s = -24..+24` along that canonical normal, mapped through the frozen
chain, trilinear on the 1.129um volume.

## 4. Eligibility and node cap

A centre row is **eligible** iff at least one profile sample maps inside the
1.129um array (same corner test as 5A, no pad). The independent scan only
covers a sliver of `w00`; a uniform subsample of all 5396 centre rows would
mostly score uncovered cells and is not a measurement of the band.

If more than 400 eligible rows exist, keep a deterministic subsample:
sort by (cell_row, cell_col) and take indices floor(i * L / 400) for
i = 0..399. Chunk self-cap 2200 (same number as Phase 7). Exceeding the cap
is `TECHNICAL_FAIL_BUDGET`, not a scientific statement.

Addendum 2026-08-25: first `--plan` used the uniform 400-cap on all centre
rows (389/400 `outside_reference_footprint`). Eligibility now filters to the
footprint first.

## 5. Per-anchor measurement

- Evaluable: masked-zero fraction of the sampled window <= 0.10 and at least
  one qualifying |gradient| peak (25% prominence), same constants as 5A.
- **D** = |s*| of the qualifying surface nearest to s = 0, in 1.129um voxels.
- **FWHM** = linear-interpolated width of |gradient| about that peak;
  censored at 48 if either side never crosses half-max.

No anisotropy gate: the normal is theirs, not a structure-tensor estimate.
