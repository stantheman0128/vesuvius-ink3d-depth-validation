@echo off
set PY=C:\Users\stans\Projects\math-science-bounty-lab\experiments\vesuvius_192\.venv\Scripts\python.exe
set STG=C:\Users\stans\Projects\villa-pr-staging
set VR=C:\Users\stans\Projects\villa\foundation\volume-registration
set PUB=C:\Users\stans\Projects\vesuvius-ink3d-depth-validation
copy /y %STG%\interface_distance.py %VR%\ || exit /b 1
copy /y %STG%\tests\test_interface_distance.py %VR%\tests\ || exit /b 1
cd /d %VR%
%PY% -m pytest tests/test_interface_distance.py -q 2>&1 || exit /b 1
echo --- schema check of the official 1.129um transform.json
%PY% -c "import json, jsonschema; s=json.load(open('transform_schema.json')); d=json.load(open(r'%PUB%\data_crossscan_meta_2026-08-07\transform_20260608103018.json')); jsonschema.validate(d, s); print('official transform.json validates against villa schema; fixed_volume =', d['fixed_volume'])" || exit /b 1
echo --- local mirror of cached 1.129um chunks
%PY% %STG%\build_local_store.py %PUB%\data_phase5b_2026-08-07\ref_chunks %STG%\local_1129_level0 || exit /b 1
echo --- real-data run: khj1222 w00 band centres vs 1.129um rescan
%PY% interface_distance.py --points %PUB%\data_phase8_khj1222_2026-08-25\w00_20231016151002_depth_anchors.csv --step inv:%PUB%\data_crossscan_meta_2026-08-07\transform_20260608103018.json --moving-offset-zyx 0 -1 0 --reference %STG%\local_1129_level0 --normal points --output %STG%\w00_interface_distance.json --csv %STG%\w00_interface_distance.csv 2>%STG%\run_stderr.txt || exit /b 1
