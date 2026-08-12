# Phase 5C code freeze addendum

Date: 2026-08-07, Asia/Taipei.
Preregistration: `PHASE5C_PREREGISTRATION_CONTENT_PRESCREEN_2026-08-07.md`.

Frozen implementation, SHA-256 computed before any 5C body acquisition:

- `phase5c_scaleup_run.py`
  `d083ec282652ee7a4fe0a14a442d1f90044360cd25965146f79be40ee6139c37`
- `tests/test_phase5c_rules.py`
  `8f243d9317fefc3987160803e532f69a3a5bdbd228b86db8539e04877aec2635`

Imported frozen modules unchanged at freeze time:

- `phase5a_depth_attribution_run.py`
  `18af22711e7960bfe0e94714f52a2d503ab1b2238b9b1f90f88ff6bef9327011`
- `phase5b_scaleup_run.py`
  `fb88a101eafeec41b6ea6f471f2a10b592de87ff732b9f13bccfcae2bdf43850`

Test state at freeze: `16 passed` across the 5A phantom, 5B rule, and 5C
prescreen suites.

Exact single run command:

```text
.venv\Scripts\python.exe phase5c_scaleup_run.py
```

Environment unchanged. Resumability as frozen in 5B.

## Amendment 1 (2026-08-07, same day)

The first invocation crashed with a `TypeError` in `plan_selection_5c` before
any 5C acquisition or selection output was written: the frozen 5A helper
`canonical_to_ref_zyx` always returns a (N,3) array and the probe-chunk index
computation omitted the `[0]` row selection. This is a technical defect fix in
untested glue, not a rule or threshold change. New SHA-256 of
`phase5c_scaleup_run.py`:
`afcfc87d358f92d1e6a269ff6eaf59b532a7b4a2f169cd241786914881fc3ae0`

Prescreen predicate, selection walk order, exclusion set, and all analysis
rules are unchanged. Rerun uses the same frozen command.
