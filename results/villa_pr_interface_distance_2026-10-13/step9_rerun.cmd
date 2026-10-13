@echo off
rem Rerun the three PR table rows on the completed mirror. Outputs overwrite the
rem October-2 runs; the old files are kept with a _pre-fetch suffix for diffing.
set PY=C:\Users\stans\Projects\math-science-bounty-lab\experiments\vesuvius_192\.venv\Scripts\python.exe
set STG=C:\Users\stans\Projects\villa-pr-staging
set VR=C:\Users\stans\Projects\villa\foundation\volume-registration
set PUB=C:\Users\stans\Projects\vesuvius-ink3d-depth-validation
set PTS=%PUB%\data_phase8_khj1222_2026-08-25\w00_20231016151002_depth_anchors.csv
set TR=inv:%PUB%\data_crossscan_meta_2026-08-07\transform_20260608103018.json
cd /d %VR%
for %%f in (w00_interface_distance w00_interface_distance_sigma0 w00_interface_distance_volnormal) do (
  if exist %STG%\%%f.json if not exist %STG%\%%f_pre-fetch.json copy /y %STG%\%%f.json %STG%\%%f_pre-fetch.json >nul
)
if exist %STG%\w00_interface_distance.csv if not exist %STG%\w00_interface_distance_pre-fetch.csv copy /y %STG%\w00_interface_distance.csv %STG%\w00_interface_distance_pre-fetch.csv >nul
if exist %STG%\w00_interface_distance_sigma0.csv if not exist %STG%\w00_interface_distance_sigma0_pre-fetch.csv copy /y %STG%\w00_interface_distance_sigma0.csv %STG%\w00_interface_distance_sigma0_pre-fetch.csv >nul
echo --- row 1: points normal, default smoothing (sigma 2)
%PY% interface_distance.py --points %PTS% --step %TR% --moving-offset-zyx 0 -1 0 --reference %STG%\local_1129_level0 --normal points --output %STG%\w00_interface_distance.json --csv %STG%\w00_interface_distance.csv 2>%STG%\run_stderr.txt || exit /b 1
echo --- row 2: points normal, raw profile (sigma 0)
%PY% interface_distance.py --points %PTS% --step %TR% --moving-offset-zyx 0 -1 0 --reference %STG%\local_1129_level0 --normal points --smooth-sigma 0 --output %STG%\w00_interface_distance_sigma0.json --csv %STG%\w00_interface_distance_sigma0.csv 2>%STG%\run_stderr_sigma0.txt || exit /b 1
echo --- row 3: volume normal (structure tensor), default smoothing
%PY% interface_distance.py --points %PTS% --step %TR% --moving-offset-zyx 0 -1 0 --reference %STG%\local_1129_level0 --normal volume --output %STG%\w00_interface_distance_volnormal.json 2>%STG%\run_stderr_vol.txt || exit /b 1
echo --- summary lines (old _pre-fetch beside new)
%PY% %STG%\summarize_runs.py
