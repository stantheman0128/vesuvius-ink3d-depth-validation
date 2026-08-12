# Phase 5A preregistration: depth-attribution pilot against the 1.129um independent 3D reference

Date frozen: 2026-08-07, Asia/Taipei.
Status: **DRAFT_PENDING_APPROVAL**. No data body may be acquired and no analysis
may touch reference data at candidate locations under this document until (a)
Stan approves this design and (b) the code-freeze addendum described in
section 11 exists.

Phase 5A data bodies downloaded at freeze: 0 bytes beyond the already-journaled
Phase A QA download (`data_crossscan_1129_p50_2026-08-07/`, 75.5 MB, aggregate
NCC use only).

## 1. Purpose and hard claim boundary

Issue #192 asks for true-3D ink labels. The maintainer rejection of PR #1295
(erdpx, 2026-08-06) states the operative bar: assigned depths must be
"validated against an independent 3d reference". This pilot asks one narrow
question:

> Do the eight highest-prediction anchors of the P50_INTERIOR cube sit at a
> papyrus surface as measured in an independently scanned, independently
> reconstructed 1.129um volume, and is the prediction depth profile thin,
> rather than smeared through the sheet?

This pilot does **not** establish ink presence, anatomical truth, issue #192
completion, or prize eligibility. Carbon ink identity is not measurable by this
construct. A pass means only: the model's depth attribution is geometrically
consistent with an independent 3D reference at the tested anchors. Surface
consistency and ink identity must never be conflated in any report.

## 2. Relationship to Phase 4B

Phase 4B ended `NO_GO_CURRENT_ROUTE` because its frozen matching rule produced
zero prediction-low controls. This document does not resurrect that design and
does not reuse its gates. The question here is different (depth attribution,
not candidate/control enrichment), the comparator is different (an analytic
null over depth positions, not matched control crops), and the evidence source
is different (an independent scan, not human axial review). The Phase 4B
candidate anchors are reused unchanged as the test set; their selection rule
(top prediction values, Chebyshev-64 separation, eight per cube) was frozen on
2026-07-27 before any Phase 5 concept existed, so reusing them adds no new
selection freedom.

## 3. Test set (frozen)

Cube: P50_INTERIOR, prediction chunk (z,y,x) = (148, 58, 75) of
`PHercParis4/representations/predictions/ink-3d/20260411134726-ink3d-20260428123845-v3-78k-fullsup.zarr`,
canonical-volume origin voxel (37888, 14848, 19200), all in the 78keV 2.4um
volume `20260411134726`. Local anchor centers (z,y,x), from the frozen
custodian bundle (sha-recorded in the 2026-08-06 session):

| pair_id | local center | global center (z,y,x) |
|---|---|---|
| P50:1 | 118,136,173 | 38006, 14984, 19373 |
| P50:2 | 136,186, 82 | 38024, 15034, 19282 |
| P50:3 | 182,156,172 | 38070, 15004, 19372 |
| P50:4 | 129, 64,144 | 38017, 14912, 19344 |
| P50:5 | 102, 74, 80 | 37990, 14922, 19280 |
| P50:6 | 64, 64,144 | 37952, 14912, 19344 |
| P50:7 | 64,138, 64 | 37952, 14986, 19264 |
| P50:8 | 166, 64, 65 | 38054, 14912, 19265 |

P25 and P75 are outside the 1.129um reference and are excluded by geometry,
recorded here so their exclusion cannot be read as case deletion.

## 4. Independent reference and registration chain (frozen)

Reference: `PHercParis4/volumes/20260608103018-1.129um-0.2m-78keV-masked.zarr`,
level 0, shape (59969, 36006, 32354), uint8, chunks 128^3, compressor null,
fill_value 0. Independence: separate scan session (2026-06), separate optics
and magnification, separate reconstruction and masking run. Shared: same
beamline and energy, same physical scroll. This is therefore a
scan-and-processing-independent geometric reference, not an energy-independent
one; the dual-energy 137keV chain is a separate future design.

Mapping: official `transform.json` of the reference
(x_canonical_xyz = A @ x_ref_xyz + t, axis order (x,y,z)), composed with the
Phase A NCC-measured residual offset of (z,y,x) = (0, -1, 0) canonical voxels.
The composed map is frozen; no per-candidate re-registration is allowed.
Phase A measured NCC 0.950 at that offset with a sharp off-edge peak; landmark
residual of the official transform is mean 2.73um. Total expected registration
error is about one 2.4um voxel.

Disclosure: the analyst has seen (a) aggregate NCC values of the central 128^3
of P50 and (b) the candidate coordinates. The analyst has not computed or
viewed any reference intensity profile, surface position, or prediction
profile at any candidate location. All constants in section 6 were fixed
before any such computation. This is an unblinded computational pilot; no
human-reviewer claim is made.

## 5. Acquisition rule (fail-closed)

1. Compute, from the frozen geometry alone, the set of reference 128^3 chunks
   intersecting the union of per-candidate boxes: global center +/- 34 voxels
   (2.4um grid) mapped through the frozen chain, padded 8 reference voxels.
2. Acquire exactly that set (GET, 3 bounded retries), journaling key, bytes,
   SHA-256, HTTP status. 404 means an all-zero masked chunk and is filled with
   fill_value; any other non-200 aborts with `TECHNICAL_FAIL_ACQUISITION`.
3. A chunk with byte length other than 2,097,152 aborts.
4. Self-cap: if the computed set exceeds 500 chunks (about 1.0 GB), stop with
   `TECHNICAL_FAIL_BUDGET` and amend rather than trim silently.
5. Chunks already journaled by Phase A may be reused only if their recorded
   SHA-256 matches a fresh re-read of the local file; mismatch aborts.
6. No other object may be downloaded under this preregistration.

## 6. Frozen analysis definitions

All computation on the 2.4um canonical grid, (z,y,x) voxel units, trilinear
sampling of the reference through the frozen chain. Software: CPython 3.13.3,
NumPy 2.3.1, SciPy 1.16.2, single deterministic run, no seeds, no tuning.

- **Smoothed reference** R_s: Gaussian sigma 2.0 voxels.
- **Normal estimation:** structure tensor of R_s (gradient sigma 2.0,
  integration sigma 4.0) at the anchor; normal = eigenvector of the largest
  eigenvalue. Anisotropy requirement: lambda1 >= 2 * lambda2.
- **Profile:** R_s and the raw prediction cube sampled at integer steps
  s = -24..+24 voxels along the normal through the anchor (57.6um half-length).
- **Surface positions:** local maxima of |d/ds R_s(s)| with prominence >= 25%
  of the profile's maximum |gradient|. Qualifying surfaces are such maxima.
- **D (localization error):** D = min over qualifying surfaces of
  |s_surface|, the distance from the anchor (s = 0) to the nearest qualifying
  surface. Threshold: D <= 3 voxels (7.2um).
- **Prediction peak and FWHM:** peak = argmax of the prediction profile;
  FWHM = linear-interpolated width at (max + min)/2 of the profile; if the
  profile never falls below half height inside +/-24, FWHM is censored at 48
  voxels and counts as a thinness failure.
- **Evaluable anchor:** profile fully inside acquired data; masked-zero
  fraction of the sampled reference window <= 10%; anisotropy holds; at least
  one qualifying surface exists in +/-24.

## 7. Hypotheses and gates

- **G1 (evaluability):** at least 4 of 8 anchors evaluable. Fewer stops the
  pilot as `INCONCLUSIVE_INSUFFICIENT_EVALUABLE`.
- **G2 (surface localization):** K = number of evaluable anchors with
  D <= 3. Null: anchor depth is arbitrary within the +/-24 band with both
  sheet interfaces available, modeled conservatively as per-anchor success
  probability p0 = 2/7. Gate: K meets the exact binomial one-sided
  alpha = 0.05 threshold frozen here:

  | n_eval | 4 | 5 | 6 | 7 | 8 |
  |---|---|---|---|---|---|
  | K required | 4 | 4 | 5 | 5 | 5 |

  (Exact tail probabilities 0.0067, 0.0257, 0.0087, 0.0233, 0.0476.)
  A secondary descriptive report uses p0 = 7/49 (thresholds 3,3,3,4,4) and
  must be labeled descriptive.
- **G3 (depth thinness):** median FWHM over evaluable anchors <= 12 voxels
  (28.8um, the order of a physical ink layer). Censored values enter the
  median at 48.

## 8. Outcomes (exhaustive)

- **O1 GO_DEPTH_ATTRIBUTION_SUPPORTED:** G1, G2, G3 all pass. Next step is a
  scale-up preregistration (more cubes inside the 1.129um footprint, and the
  dual-energy 137keV chain), not a submission claim.
- **O2 DEPTH_SMEARED:** G1 and G2 pass, G3 fails. The prediction localizes at
  a surface but is smeared in depth. Documented as a quantified instance of
  the exact defect #192 names; potentially submittable as a diagnostic.
- **O3 NO_SURFACE_PREFERENCE:** G1 passes, G2 fails. The depth attribution is
  not distinguishable from arbitrary placement at these anchors. Honest
  negative; route review before any further spend.
- **O4 INCONCLUSIVE_INSUFFICIENT_EVALUABLE:** G1 fails.
- **O5 TECHNICAL_FAIL_*:** acquisition, identity, or budget failure before
  any scientific statement; the scientific question stays unanswered.

Whatever the outcome, all eight per-anchor records (D, FWHM, evaluability,
reasons) are reported; no anchor may be dropped post hoc.

## 9. Explicit alternative explanations

A G2 pass can reflect: the model tracking readable papyrus surfaces rather
than ink; surfaces being dense in the sampled region (partially absorbed by
the conservative p0); registration luck (bounded by Phase A residual); or the
anchor selection rule preferring surface-adjacent voxels for reasons unrelated
to ink. A pass therefore supports depth attribution, not ink detection.

## 10. Prohibited moves

No threshold tuning after seeing profiles. No anchor replacement or addition.
No switching cubes. No re-registration per candidate. No renaming of surface
consistency as ink truth. No second statistic swapped in after unblinding. Any
deviation requires a written amendment before the deviating step runs.

## 11. Code commitment

Before the first run on real candidate data: implement
`phase5a_depth_attribution_run.py` plus a phantom test
(`tests/test_phase5a_phantom.py`) that (a) recovers a known synthetic surface
distance within 0.5 voxel, (b) recovers a known synthetic FWHM within 1
voxel, and (c) correctly censors a non-dropping profile. Freeze both files'
SHA-256 in `PHASE5A_CODE_FREEZE_ADDENDUM_2026-08-07.md` together with the
exact command line. Only then acquire and run, once.
