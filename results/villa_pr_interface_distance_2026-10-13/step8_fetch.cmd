@echo off
rem Download the 366 level-0 chunks the w00 read boxes need (about 730 MB) and
rem link them into local_1129_level0. Safe to rerun: cached chunks are reused.
set PY=C:\Users\stans\Projects\math-science-bounty-lab\experiments\vesuvius_192\.venv\Scripts\python.exe
set STG=C:\Users\stans\Projects\villa-pr-staging
%PY% %STG%\fetch_missing_chunks.py || exit /b 1
echo --- re-audit: should report 0 distinct_box_chunks_missing unless some chunks are 404 (masked)
%PY% %STG%\audit_missing_chunks.py || exit /b 1
