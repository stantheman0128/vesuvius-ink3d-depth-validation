# Validating 3D ink-label depth against an independent scan: method and first measurements on PHercParis4

Status: submission draft, pending Stan's review. Not yet published anywhere.
Date: 2026-08-08. All analysis ran on a desktop CPU; total download 3.9 GB
from the official open-data bucket.

## Summary

Issue #192 asks for true-3D ink labels. When PR #1295 was closed on
2026-08-06, the maintainer stated the missing piece: assigned depths must be
validated against an independent 3D reference. This work supplies that
validation as a runnable, preregistered pipeline and reports first
measurements. Using the official landmark transforms that ship inside the
PHercParis4 volume zarrs, we register the canonical 2.4um 78keV volume to the
independently scanned and reconstructed 1.129um volume, verify the chain by
normalized cross-correlation to about one 2.4um voxel, and then test whether
the depth positions asserted by the official ink-3d prediction coincide with
papyrus surfaces measured only in the independent scan.

Across 24 preregistered cubes disjoint from all earlier samples (192 anchors,
142 evaluable), 99 of 142 anchors sit within 3 voxels (7.2um) of an
independently measured surface. Under a conservative chance model in which a
depth claim could land anywhere in the sampled band with both sheet faces
counted as favorable, the exact one-sided binomial p-value is 4.1e-24. At the
same time, 24 of 142 anchors carry depth profiles wider than 12 voxels
(28.8um), and 43 are not surface-localized at all. The depth quality of the
current labels is therefore mixed in a way this pipeline can now measure
point by point.

## What the maintainers asked for

The #192 issue text asks for labels that are "in true 3d rather than a single
image projected across multiple layers". The closing comment on PR #1295 says
depth windows placed by CT brightness "may produce useful experimental
pseudo-labels", that quality control which only shows assigned depths vary is
insufficient, that 2D cross-scan comparison "cannot validate the 3d depth",
and that the work could be reconsidered "if the full set can be published and
the depths validated against an independent 3d reference". This report is a
direct response to that last clause.

## The unused asset: official cross-scan transforms

Each PHercParis4 volume zarr contains a `transform.json` with a 3x4 affine
and paired landmark sets. The 2.4um 78keV canonical volume and the 2.4um
137keV volume both register to the original 7.81um volume
(PHercParis4-20230205180739); the 1.129um 78keV volume registers directly to
the canonical 2.4um volume. The issue and PR searches performed for this work
did not surface community use of these files for label validation, though we
have not audited every community repository.

Landmark residuals: 2.73um mean for the 1.129um transform, 5.5um mean for the
137keV transform, 14.9um mean (50um worst case) for the 78keV transform.

## Method in brief

Full preregistrations, frozen constants, code hashes, phantom tests, and
acquisition journals are in the repository. The chain of experiments:

1. **Registration QA (Phase A).** Resample the 1.129um volume onto the 2.4um
   grid through the official transform and cross-correlate against an
   already-acquired raw 2.4um cube. Result: NCC peak 0.950 at an offset of
   one voxel, off the search edge. The composed transform plus this fixed
   one-voxel correction is then frozen for all later phases.
2. **Depth-attribution construct.** For each anchor (a local maximum of the
   official ink-3d prediction under the frozen Phase 4B sampling rule), take
   the surface normal from the structure tensor of the resampled independent
   volume, extract a 49-sample profile along that normal, and measure
   D, the distance from the anchor to the nearest qualifying gradient peak
   (an independently observed surface), and the FWHM of the prediction
   profile along the same line. All thresholds (D <= 3 voxels, FWHM <= 12
   voxels, prominence, anisotropy, masked-content caps) were frozen before
   any real profile was computed, and the analysis code passed synthetic
   phantom tests (surface recovery within 0.5 voxel, FWHM within 1 voxel)
   before its hash was frozen.
3. **Pilot (Phase 5A, n = 8).** One cube. Surface localization passed exactly
   at its preregistered threshold (5 of 8 anchors, exact p = 0.048); depth
   thinness failed (median FWHM 17.3 voxels). Preregistered outcome:
   DEPTH_SMEARED.
4. **Scale-up (Phase 5B, 24 cubes).** Ended inconclusive by its own rule: 89
   evaluable anchors against a floor of 96, because 9 of 24 cubes fell in
   masked air. The design defect (bounds checked, mask content not) is
   journaled, and 5B numbers are reported as descriptive only.
5. **Confirmatory rerun (Phase 5C, 24 fresh cubes).** Adds a masked-content
   prescreen and excludes every previously seen cube. 142 of 192 anchors
   evaluable. K = 99 with D <= 3, exact p = 4.1e-24 against the conservative
   chance model. Preregistered outcome: LOCALIZATION_CONFIRMED.

Class proportions among the 142 evaluable 5C anchors, Wilson 95% intervals:

| Class | n | Share | Wilson 95% |
|---|---:|---:|---|
| Localized and thin (D <= 3, FWHM <= 12) | 75 | 52.8% | 44.6% to 60.8% |
| Localized but depth-smeared (D <= 3, FWHM > 12) | 24 | 16.9% | 11.6% to 23.9% |
| Not surface-localized (D > 3) | 43 | 30.3% | 23.3% to 38.3% |

Non-evaluable anchors: 50 of 192 (26.0%), each with journaled reasons. The
5B descriptive sample, disjoint from 5C, shows the same structure: 76.4% of
its evaluable anchors localized, 51.7% localized and thin.

## Dual-energy chain

The 137keV volume registers through the shared 7.81um fixed volume, which
gives a second, energy-independent reference covering a band of the scroll
(canonical z roughly 42,100 to 48,800). A QA cross-correlation on one cube in
that band gives NCC 0.924 at an offset of about 5 voxels, so this chain is
usable after a local NCC correction. Depth validation against a
different-energy scan is the natural next preregistration; this report only
measures the chain quality.

## What this does and does not show

It shows, with a preregistered confirmatory test on data disjoint from every
earlier phase, that the depth positions asserted by the current ink-3d
prediction agree with independently measured papyrus surfaces far more often
than a conservative chance model allows, and that roughly half of the tested
anchors carry depth profiles as thin as a physical ink layer. It also
quantifies the defect #192 names: about one sixth of evaluable anchors are
surface-localized but smeared in depth, and about three in ten are not
surface-localized at all.

It does not show ink identity. Carbon ink contrast is not directly measurable
by this construct, so surface consistency must not be renamed ink truth. It
does not cover the whole scroll: anchors are prediction-selected maxima
inside the 1.129um footprint, so nothing here estimates recall or false
negatives. The 2.4um voxel grid bounds the precision of both D and FWHM, and
the conservative chance model is a modeling choice, with the sensitivity
alternative reported in the frozen preregistrations.

## Why this is useful for #192

A label set that ships with per-voxel D and FWHM measured against an
independent scan answers the reopening condition stated on PR #1295. The
pipeline here runs on a desktop CPU, downloads only the chunks it needs, and
fails closed on identity drift, coverage gaps, and budget caps. Applied at
scale it can (a) filter existing labels to their depth-validated subset, (b)
score any candidate 3D label set produced by others, and (c) extend to the
dual-energy chain for a reference that is independent in physics as well as
in scan session.

## Relation to concurrent work on this issue

After our runs completed, khj1222 reported on this issue (2026-08-08) a
training comparison in which a per-pixel measured depth band loses 0.038 F1 to
a constant band on `w00_20231016151002`, with a one-voxel plane and an
eight-voxel constant band tied. That work measures the downstream training
utility of label shapes; ours measures whether asserted depths are
geometrically real against an independent scan, which is the reopening
condition stated on PR #1295. The two results point the same way from
different directions: their per-pixel depths came from model profiling and
lost to a constant band, and we find that roughly three in ten
prediction-derived anchors are not surface-localized, with another sixth
smeared in depth. Their depth source differs from the prediction volume we
audited, so we treat this as consistent context rather than an explanation.
A natural joint next step is to score any candidate band construction with
this pipeline before spending a training run on it.

## Reproducibility

Every phase has a preregistration file frozen before its data was touched, a
code-freeze addendum with SHA-256 hashes, phantom or rule tests (16 passing),
and a JSON result with full per-anchor records. Phase A and 5A persisted full
per-request HTTP journals; 5B and 5C persisted selection, probe, and QA
journals plus the chunk cache itself, and a post-hoc SHA-256 inventory of
every cached chunk (computed 2026-08-13, labeled as such) is included. Total spend: 3.9 GB of downloads, no
GPU, no model training, and no tuning after unblinding at any phase. One
technical defect (an array-shape error in untested selection glue) is
documented in an amendment; it crashed before any acquisition and its fix
changed no rule.
