"""Scratch helper (not part of the PR): expose the cached 1.129um chunks as a
zarr v2 array directory so interface_distance.py can read them like the bucket.

Cached files are raw uint8 128^3 blocks named chunk_z_y_x.bin (zero-length for
404/masked chunks). Hard links, no copies.
"""
import json
import os
import sys

src, dst = sys.argv[1], sys.argv[2]
os.makedirs(dst, exist_ok=True)
with open(os.path.join(dst, ".zarray"), "w") as f:
    json.dump({"zarr_format": 2, "shape": [59969, 36006, 32354],
               "chunks": [128, 128, 128], "dtype": "|u1", "compressor": None,
               "fill_value": 0, "order": "C", "filters": None,
               "dimension_separator": "."}, f)
n = 0
for name in os.listdir(src):
    if not name.startswith("chunk_") or not name.endswith(".bin"):
        continue
    p = os.path.join(src, name)
    if os.path.getsize(p) != 128 ** 3:
        continue
    z, y, x = name[6:-4].split("_")
    out = os.path.join(dst, f"{z}.{y}.{x}")
    if not os.path.exists(out):
        os.link(p, out)
    n += 1
print(f"linked {n} chunks into {dst}")
