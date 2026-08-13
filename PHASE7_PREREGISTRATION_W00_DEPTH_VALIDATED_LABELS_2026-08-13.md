# Phase 7 preregistration: depth-validated 3D ink labels for w00_20231016151002

Date frozen: 2026-08-13, Asia/Taipei.
Status: approved continuation (Stan's ordered list item 1, 2026-08-13).
Execution requires the code-freeze addendum of section 9 before any new body
acquisition.

Inherits the frozen registration chain of Phase 5A section 4 (official
1.129um transform composed with the Phase A one-voxel offset) and the
acquisition identity rules of 5A section 5. Statistical machinery from 5B/5C
is not inherited: this phase produces an artifact with QC, not a confirmatory
test.

## 1. Purpose and claim boundary

Issue #192 asks for ready-to-run image/label pairs that are true 3D. This
phase emits, for segment `w00_20231016151002`, a 3D label volume whose depth
placement is measured against the independently scanned 1.129um volume at a
preregistered grid of nodes, together with per-node QC. Claims are geometric
only: the labels' 2D ink footprint is the published human annotation,
unchanged; ink identity is not measured; training benefit is not claimed
(khj1222's 2026-08-08 result on this segment shows a per-pixel measured band
losing to a constant band, and nothing here presumes the opposite).

## 2. Inputs and identity (fail-closed)

1. Local, already acquired (SHA-256 must match the phase 3 inspection record
   `outputs/phase3_2d_maps_header_inspection_2026-07-26.json`):
   - `w00_20231016151002_inklabels.tif`
     `31810cd7fd8948618a29232cd0839cf249f0f3906a8218d3408156f315c2d014`
     (32249 x 51380, uint8)
   - `w00_20231016151002_supervision_mask.tif`
     `70d5f3d7c849fdc24c04e903791adf047ec95f590246b7472e33f411ebc19a1a`
2. To acquire, journaled (bytes + SHA-256): `x.tif`, `y.tif`, `z.tif`
   (float32 coordinate maps at scale about 0.1; expected shape near
   3225 x 5138). Their coordinate frame was verified from the local
   `meta.json`: values fit the canonical 2.4um volume (bbox z up to 73948
   against a canonical extent of 75784) and not the 7.81um volume, so map
   values are treated as canonical (x, y, z) voxel coordinates. The maps mark
   invalid pixels with the value -1 (about 16% of pixels, verified after
   download). A map pixel is **valid** iff all three coordinates are finite
   and >= 0. Fail-closed check at run time: at least 99% of valid coordinates
   must lie inside the canonical volume bounds, else stop with
   `TECHNICAL_FAIL_FRAME`.
3. The reference volume and chain: exactly as frozen in 5A. Reference chunks
   go through the shared journaled cache.

## 3. Node grid (deterministic)

1. Work on the coordinate-map grid (H x W near 3225 x 5138). Node stride: 8
   map pixels in each direction (about 80 full-resolution segment pixels).
2. A node is **eligible** iff:
   - its 8 x 8 map-pixel cell, upscaled to the full-resolution ink label,
     contains at least 1% ink-positive pixels;
   - its map coordinates and those of its 4-neighbours are finite;
   - the canonical box (node coordinate +/- 34 voxels), mapped through the
     frozen chain with the 8-voxel pad, lies inside the reference array
     (same corner test as 5B).
3. If more than 800 nodes are eligible, keep a deterministic subsample:
   sort nodes lexicographically by (row, col) and take indices
   floor(i * L / 800) for i = 0..799. No other selection freedom.

## 4. Per-node measurement

- **Normal:** unit cross product of the coordinate map's central-difference
  tangents at the node (d/dcol x d/drow of the (x,y,z) map). This is the
  segment's own surface normal; the structure-tensor estimator of 5A is not
  used here, so no anisotropy gate applies.
- **Profile:** the reference volume resampled on the canonical grid at
  integer steps s = -24..+24 along that normal through the node's canonical
  coordinate, exactly as in 5A section 6 (Gaussian sigma 2.0 smoothing,
  trilinear sampling, gradient-peak surfaces with 25% prominence).
- **s\*** = signed s of the qualifying surface nearest to s = 0 (the
  independently measured surface position relative to the published one).
  **D** = |s\*|.
- **Evaluable node:** masked-zero fraction of the sampled window <= 10% and
  at least one qualifying surface exists. Non-evaluable nodes are recorded
  with reasons and excluded from correction fitting.

## 5. Label artifact (frozen construction)

1. Read the layer count L_z of the official surface volume
   `w00_20231016151002.zarr` from its `.zarray`/`.zattrs` metadata (metadata
   request only); L_mid = floor(L_z / 2). A layer step is treated as one
   canonical voxel along the surface normal. If metadata is absent or
   inconsistent, stop with `TECHNICAL_FAIL_SURFACE_META`.
2. Correction field: s\*(u, v) interpolated bilinearly from evaluable nodes
   onto the map grid, clamped to [-8, +8] voxels; zero where the nearest
   evaluable node is farther than 24 map pixels (flagged `unvalidated`).
3. Output zarr `w00_20231016151002_inklabels_v5.zarr`, full segment
   resolution, shape (L_z, 32249, 51380), uint8, chunks (L_z, 256, 256):
   for each full-resolution ink-positive pixel, label voxels
   L_mid + round(interp s\*) +/- 4 layers (a 9-layer constant-width band,
   clamped to [0, L_z-1]) are 255; everything else 0. Band width is fixed at
   the khj1222 v3-equivalent half-width 4; width is deliberately not
   per-pixel.
4. QC sidecar `w00_20231016151002_inklabels_v5_qc.zarr` at map resolution,
   float32 channels: s\* (interpolated), D at nodes (NaN elsewhere),
   validation status code (0 unvalidated, 1 interpolated, 2 node).
5. Everything is emitted under `data_phase7_w00_2026-08-13/`; nothing is
   uploaded anywhere under this preregistration.

## 6. Reported statistics (descriptive only)

Node count, evaluable count and reasons; distribution of s\* and D
(quartiles); fraction of evaluable nodes with D <= 3 voxels; fraction of
full-resolution ink pixels covered by status 1 or 2. No confirmatory gate:
the 5C result already established localization, and this phase's deliverable
is the artifact plus honest QC.

## 7. Budget and caps

- New reference chunks for node windows: hard cap 2,200 (about 4.6 GB);
  exceeding the plan aborts with `TECHNICAL_FAIL_BUDGET` before download.
- Coordinate tifs: about 200 MB, journaled.
- Everything runs on CPU; resumability per the 5B frozen design.

## 8. Honest limitations (to ship with the artifact)

Node density is sparse (one node per about 0.19 mm of segment surface at
stride 8, further thinned by the 800-node cap), so the correction field is a
smooth interpolation, not a per-pixel measurement. Coverage is partial: only
the part of the segment inside the 1.129um footprint (canonical z roughly
46,985 to 56,359 for this segment) can be validated; the rest of the segment
keeps the published plane and is flagged. The published 2D footprint is
trusted as-is.

## 9. Code commitment

`phase7_w00_depth_labels_run.py` plus unit tests in
`tests/test_phase7_rules.py` covering: normal computation on a synthetic
planar map, node eligibility and the 800-cap stride, correction clamping and
the unvalidated radius, and band construction at layer edges. Freeze SHA-256
values and the exact command in
`PHASE7_CODE_FREEZE_ADDENDUM_2026-08-13.md`; the 5A phantom suite must still
pass; then run once (resumable).
