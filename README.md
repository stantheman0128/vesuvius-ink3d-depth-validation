# Depth-attribution validation for 3D ink labels (villa #192)

This repository validates the depth component of 3D ink labels for the
Vesuvius Challenge against an independent scan, using the official landmark
transforms that ship inside the PHercParis4 volume zarrs. It responds to the
reopening condition stated when PR ScrollPrize/villa#1295 was closed: depths
must be "validated against an independent 3d reference".

Headline result: across 24 preregistered cubes (192 anchors, 142 evaluable),
99 of 142 anchors asserted by the official ink-3d prediction sit within 3
voxels (7.2um) of a papyrus surface measured only in the independently
scanned 1.129um volume; the exact one-sided binomial p-value against a
conservative chance model is 4.1e-24. About half the anchors carry depth
profiles as thin as a physical ink layer, one sixth are surface-localized but
depth-smeared, and three in ten are not surface-localized. Full method,
numbers, and claim limits are in [REPORT.md](REPORT.md).

This validates depth geometry only. It never measures ink identity, and the
report states that boundary explicitly.

## Layout

- `REPORT.md`: the full write-up.
- `PHASE5*_PREREGISTRATION_*.md`, `PHASE5*_CODE_FREEZE_ADDENDUM_*.md`: frozen
  designs, constants, gates, and code hashes, written before each run.
- `phaseA_crossres_registration_qa_2026-08-07.py`: registration chain QA
  (NCC 0.950 at a one-voxel offset).
- `crossscan_overlap_analysis_2026-08-07.py`: transform audit and footprint
  overlap computation.
- `phase5a_depth_attribution_run.py`, `phase5b_scaleup_run.py`,
  `phase5c_scaleup_run.py`: the pilot, the scale-up, and the confirmatory
  rerun. 5B and 5C fetch the chunks they need from the official open-data
  bucket and cache them locally.
- `tests/`: phantom and rule tests (16 tests; surface recovery within 0.5
  voxel, FWHM within 1 voxel, sampler tie policy, NCC shift recovery).
- `data_crossscan_meta_2026-08-07/`: copies of the three official
  `transform.json` files the chain is built on.
- `results/`: result JSONs with full per-anchor records, acquisition and
  selection journals, and post-hoc SHA-256 inventories of the chunk caches.

## Running

Python 3.13 with numpy, scipy, numcodecs. No GPU. From the repository root:

```bash
python -m pytest tests -q
python phase5b_scaleup_run.py --plan
python phase5c_scaleup_run.py
```

`phase5c_scaleup_run.py` reproduces the confirmatory result end to end,
downloading about 2 GB of chunks on first run. `phase5a_depth_attribution_run.py`
additionally expects the eight raw and one prediction chunk of the P50 cube
(chunk keys and hashes are pinned in the Phase 5A preregistration and
journals); it was run against a locally acquired snapshot that this
repository does not redistribute.

No volume data is redistributed here. Every binary chunk is fetched from the
official `vesuvius-challenge-open-data` bucket at run time, and the journals
pin the identity of everything that was fetched.

## License

MIT for all code and text in this repository. The underlying scan data
belongs to the Vesuvius Challenge and is governed by its own terms.
