# Phase 7 code freeze addendum

Date: 2026-08-13, Asia/Taipei.
Preregistration: `PHASE7_PREREGISTRATION_W00_DEPTH_VALIDATED_LABELS_2026-08-13.md`.

Frozen implementation, SHA-256 computed before any node measurement or label
emission (the only prior Phase 7 data actions were the journaled x/y/z.tif
download, the surface zarr metadata read, and the map frame checks recorded
in PROGRESS.md):

- `phase7_w00_depth_labels_run.py`
  `449fa5d554beb358d2f54706f207a9a7489cfec163232a18ffd2b5d176b7af50`
- `tests/test_phase7_rules.py`
  `74a94bd2f6727e85a47672a1b454d052e19f61ad795124aca96a30d449877062`

Test state at freeze: `12 passed` (Phase 7 rules plus the unchanged 5A
phantom suite). One pre-freeze fix during test development: the correction
interpolator now falls back to nearest-neighbour when fewer than three nodes
exist or triangulation fails; caught by the unit test on a two-node case.

Environment addition, recorded: `imagecodecs 2026.6.26` installed into the
project venv (required to decode the LZW-tiled ink label TIFF). No other
package changes.

Dry-run plan at freeze: 148 eligible nodes (below the 800 cap; the cap rule
is inert for this segment), 853 new reference chunks, about 1.8 GB, within
the 2,200-chunk cap.

Exact single run command:

```text
.venv\Scripts\python.exe phase7_w00_depth_labels_run.py
```

Resumable per the frozen design (node selection, measurements, and chunk
cache are reused on re-invocation).
