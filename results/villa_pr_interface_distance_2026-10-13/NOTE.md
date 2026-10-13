# interface_distance runs on the complete local mirror (2026-10-13)

These replace `../villa_pr_interface_distance_2026-10-03/`. A pre-submission
check of villa PR #1945 found that the local zarr mirror used on 2026-10-03
held only the 322 level-0 chunks the August scripts had fetched. The PR
script reads a larger box per point (34 point voxels, about 72 reference
voxels, each side), so for 159 of the 172 in-footprint points at least one
chunk of that box was absent and read as zero, which the script treats as
masked. `audit_missing_chunks.py` lists the chunks each point's box and
profile tube touch; before the fetch it reported 366 distinct box chunks
missing (101 for the profile tubes alone).

`fetch_missing_chunks.py` downloaded those 366 chunks from
`vesuvius-challenge-open-data` (all HTTP 200, none masked; SHA-256 per chunk
in `fetch_missing_chunks_journal.json`) into the August cache and linked them
into the mirror, now 688 chunks. `audit_after_fetch.json` shows 0 missing.
`step9_rerun.cmd` then reran the three table rows of the PR with the script at
villa commit 5628683.

| run | 2026-10-03 (322 chunks) | 2026-10-13 (688 chunks) |
|---|---|---|
| points normal, sigma 0 | 123 / 164 within 3 vox, median 2.0 | 123 / 164, median 2.0 |
| points normal, sigma 2 | 97 / 163, median 3.0 | 98 / 163, median 3.0 |
| volume normal, sigma 2 | 36 / 51, median 2.0 | 92 / 147, median 3.0 |

The points-normal rows barely move because the profile tube of most points
was already inside the 322 chunks. The volume-normal row changes because the
structure tensor uses the whole window; with the incomplete mirror 114 points
were rejected as masked, with the complete one 9 are. The 8 points rejected
as masked in the points rows hit genuinely masked (zero) voxels of the
reference and are unchanged by the fetch.
