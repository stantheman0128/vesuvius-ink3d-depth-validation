# interface_distance band fields, original and negated normals (2026-10-18)

Companion to `../villa_pr_interface_distance_2026-10-13/`. Two things changed
since that directory was written.

1. villa PR #1945 gained a third commit (`21c920a`) that adds the band fields
   `d_minus`, `d_plus`, `band_width`, `band_position` and `inside_brighter` to
   every record and `bracketed`, `median_band_width_voxels` and
   `inside_brighter` to the summary. `w00_points_sigma0.*` and
   `w00_points_sigma2.*` are the two points-normal rows rerun with that commit
   on the same 688-chunk mirror. `D_voxels` is identical id by id to the
   2026-10-13 files.
2. On villa issue #192 (2026-10-14) khj1222 pointed out that the sign of the
   exported normal is the one convention their sidecar marks as unverified,
   and that with `d_minus`/`d_plus` the sign now decides which side is minus
   and also where the point lands, since the point is
   `base + offset_from_plane * normal`. They suggested running the same CSV
   with the three `normal_*` columns negated. `negate_normals.py` does that;
   `w00_points_*_normals_negated.*` are the results.

All runs: script at villa commit `21c920a`, `--step inv:transform_20260608103018.json
--moving-offset-zyx 0 -1 0 --normal points`, reference = local level-0 mirror
of the 1.129 um PHerc. Paris 4 volume (688 chunks, see the 2026-10-13 note),
threshold 3 voxels. Every evaluable point is in region 15, the only region the
reference footprint covers.

| run | evaluable | within 3 vox | median D | bracketed | median band width | inside_brighter |
|---|---|---|---|---|---|---|
| sigma 0, normals as exported | 164 | 123 | 2.0 | 154 | 7.0 | 70 / 131 = 0.534 [0.449, 0.618] |
| sigma 0, normals negated | 163 | 113 | 2.0 | 141 | 7.0 | 60 / 112 = 0.536 [0.444, 0.625] |
| sigma 2, normals as exported | 163 | 98 | 3.0 | 146 | 11.0 | 81 / 133 = 0.609 [0.524, 0.688] |
| sigma 2, normals negated | 162 | 93 | 3.0 | 134 | 10.0 | 70 / 116 = 0.603 [0.513, 0.688] |

`compare_readings.py` compares the two readings point by point
(`compare_readings.txt`):

- Negating the normal moves a point by a median 5.5 reference voxels
  (min 0.3, max 31.3), which matches khj1222's estimate of about 5 voxels
  from the median offset in region 15.
- The within-3 verdict is the same for only 104 of 163 points (sigma 0) and
  100 of 161 (sigma 2). Equal numbers of points move farther from and closer
  to an interface (69 and 69; 67 and 67), which is why the totals barely move.
- `inside_brighter` does not depend on the sign by construction (it compares
  the mean inside the band with the mean outside it), so where it differs
  between the readings the displaced point sits in a different band. It agrees
  on 45 of 91 points (sigma 0) and 54 of 97 (sigma 2), about chance.

Reading: the aggregate within-threshold count is insensitive to the sign of
the normal, so this check cannot settle which sign is right. The two readings
are different point sets that happen to score alike. The sign stays
unverified, and anyone using the band fields on this export should carry both
readings until it is settled by something that sees the surface directly.
