# Phase 5B code freeze addendum

Date: 2026-08-07, Asia/Taipei.
Preregistration: `PHASE5B_PREREGISTRATION_SCALEUP_2026-08-07.md` (approved as
the option-3 continuation).

Frozen implementation, SHA-256 computed before any body acquisition under 5B:

- `phase5b_scaleup_run.py`
  `fb88a101eafeec41b6ea6f471f2a10b592de87ff732b9f13bccfcae2bdf43850`
- `tests/test_phase5b_rules.py`
  `f96455230e086fcad8fe0fa75ab0477dd088f2c70b388d72255b6ae1ba48546b`

Test state at freeze: `13 passed` for the combined 5B rule suite and the
unchanged 5A phantom suite
(`.venv\Scripts\python.exe -m pytest tests/test_phase5b_rules.py tests/test_phase5a_phantom.py -q`).
Two earlier failures were test-fixture defects (a second synthetic peak placed
at Chebyshev distance 60, correctly rejected by the frozen rule; a negative
slice index in the NCC fixture) and were fixed in the tests only; the frozen
rule implementations were not changed in response to any real data.

Selection plan output (deterministic, cached in
`data_phase5b_2026-08-07/cube_selection.json`): 422,859 geometrically eligible
cubes; 24 selected by the frozen stride walk, z-chunks 110 to 215.

Imports from the frozen Phase 5A module are byte-identical: 5A file SHA-256
remains `18af22711e7960bfe0e94714f52a2d503ab1b2238b9b1f90f88ff6bef9327011`.

Exact single run command:

```text
.venv\Scripts\python.exe phase5b_scaleup_run.py
```

Environment unchanged from 5A. Resumability (per-cube result files, on-disk
chunk cache) is part of the frozen design; a resumed invocation of the same
command is the same run, not a re-run.
