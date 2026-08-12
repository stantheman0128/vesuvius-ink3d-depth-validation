# Phase 5A code freeze addendum

Date: 2026-08-07, Asia/Taipei.
Preregistration: `PHASE5A_PREREGISTRATION_DEPTH_ATTRIBUTION_1129_2026-08-07.md`
(approved by Stan in session before this freeze).

Frozen implementation, SHA-256 computed before any acquisition or analysis on
real candidate data:

- `phase5a_depth_attribution_run.py`
  `18af22711e7960bfe0e94714f52a2d503ab1b2238b9b1f90f88ff6bef9327011`
- `tests/test_phase5a_phantom.py`
  `2a8ce5a29942c5f5127c15bf10110d998117999d7c6ba3183da3c80be38005af`

Phantom test result before freeze: `6 passed in 5.58s`
(`.venv\Scripts\python.exe -m pytest tests/test_phase5a_phantom.py -q`),
covering planar-surface recovery within 0.5 voxel (tilted and centered),
FWHM recovery within 1 voxel, censoring behavior, no-surface behavior, and
gate-table equality with the preregistration.

Dry-run acquisition plan: 72 reference chunks, about 151 MB, within the
500-chunk self-cap.

Exact single run command:

```text
.venv\Scripts\python.exe phase5a_depth_attribution_run.py
```

Environment: CPython 3.13.3 (.venv), NumPy 2.3.1, SciPy 1.16.2,
numcodecs 0.16.5. Deterministic, no seeds. Any code change after this freeze
requires a new addendum before re-running.
