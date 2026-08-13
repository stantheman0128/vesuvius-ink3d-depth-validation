# Phase 6 code freeze addendum

Date: 2026-08-13, Asia/Taipei.
Preregistration: `PHASE6_PREREGISTRATION_DUAL_ENERGY_2026-08-13.md`.

Frozen implementation, SHA-256 computed before any Phase 6 measurement (the
only prior Phase 6 data actions were the selection walk's journaled HEAD and
content probes):

- `phase6_dual_energy_run.py`
  `f464a48be10e8ee8137d5b6756b287593a81be4aceb11e1303faa525aac263a0`
- `tests/test_phase6_rules.py`
  `037d68f64ec8f39fa1813dd360d359a6ca705cc12885fab410845bac58a2a30b`

Test state at freeze: `9 passed` (Phase 6 rules, including a chain-direction
round trip on the official landmark pairs and the input-shift offset
convention, plus the unchanged 5A phantom suite). One pre-freeze fix during
implementation review: the per-cube NCC correction is applied as an input
shift, fn(v + offset), matching the frozen Phase A convention; the first
draft added the offset in 137keV space, which differs by the chain rotation.

Selection at freeze (cached in `data_phase6_2026-08-13/cube_selection.json`):
12 cubes spanning z-chunks 165 to 187 across the 137keV band, all passing the
prediction HEAD and 137keV content prescreen, disjoint from every 5B/5C cube.

Exact single run command:

```text
.venv\Scripts\python.exe phase6_dual_energy_run.py
```

Resumable per the frozen 5B design.
