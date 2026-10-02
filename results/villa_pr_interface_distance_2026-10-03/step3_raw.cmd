@echo off
set PY=C:\Users\stans\Projects\math-science-bounty-lab\experiments\vesuvius_192\.venv\Scripts\python.exe
set STG=C:\Users\stans\Projects\villa-pr-staging
set VR=C:\Users\stans\Projects\villa\foundation\volume-registration
set PUB=C:\Users\stans\Projects\vesuvius-ink3d-depth-validation
cd /d %VR%
%PY% interface_distance.py --points %PUB%\data_phase8_khj1222_2026-08-25\w00_20231016151002_depth_anchors.csv --step inv:%PUB%\data_crossscan_meta_2026-08-07\transform_20260608103018.json --moving-offset-zyx 0 -1 0 --reference %STG%\local_1129_level0 --normal points --smooth-sigma 0 --output %STG%\w00_interface_distance_sigma0.json --csv %STG%\w00_interface_distance_sigma0.csv 2>%STG%\run_stderr_sigma0.txt || exit /b 1
echo --- volume-normal mode, default smoothing
%PY% interface_distance.py --points %PUB%\data_phase8_khj1222_2026-08-25\w00_20231016151002_depth_anchors.csv --step inv:%PUB%\data_crossscan_meta_2026-08-07\transform_20260608103018.json --moving-offset-zyx 0 -1 0 --reference %STG%\local_1129_level0 --normal volume --output %STG%\w00_interface_distance_volnormal.json 2>%STG%\run_stderr_vol.txt || exit /b 1
