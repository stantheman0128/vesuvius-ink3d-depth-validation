"""Fetch the 1.129um level-0 chunks the w00 read boxes need but the local cache lacks.

Reads audit_missing_chunks.json (box_missing per point), downloads each chunk from
the public S3 bucket into the August cache directory using the same layout and 404
convention as phase5b_scaleup_run.cached_chunk (raw uint8 128^3 as chunk_z_y_x.bin;
zero-length file = 404/masked), writes a journal, then hard-links the new full
chunks into the zarr v2 mirror. Scratch helper, not part of the PR.
"""
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

STG = os.path.dirname(os.path.abspath(__file__))
CACHE = (r"C:\Users\stans\Projects\vesuvius-ink3d-depth-validation"
         r"\data_phase5b_2026-08-07\ref_chunks")
MIRROR = os.path.join(STG, "local_1129_level0")
REF_BASE = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/"
            "PHercParis4/volumes/20260608103018-1.129um-0.2m-78keV-masked.zarr/0")
CHUNK_BYTES = 128 ** 3

audit = json.load(open(os.path.join(STG, "audit_missing_chunks.json")))
keys = sorted({k for r in audit["rows"] for k in r["box_missing"]})
print(f"{len(keys)} chunks to fetch", flush=True)


def fetch(key):
    z, y, x = key.split(".")
    path = os.path.join(CACHE, f"chunk_{z}_{y}_{x}.bin")
    if os.path.exists(path):
        return {"key": key, "status": "local_reuse", "bytes": os.path.getsize(path)}
    url = f"{REF_BASE}/{z}/{y}/{x}"
    err = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = r.read()
            if len(data) != CHUNK_BYTES:
                raise RuntimeError(f"size {len(data)}")
            with open(path + ".part", "wb") as f:
                f.write(data)
            os.replace(path + ".part", path)
            return {"key": key, "status": 200, "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest()}
        except urllib.error.HTTPError as e:
            if e.code == 404:
                open(path, "wb").close()
                return {"key": key, "status": 404, "bytes": 0}
            err = f"http {e.code}"
        except Exception as e:  # network hiccup: retry with backoff
            err = repr(e)
        time.sleep(2 * (attempt + 1))
    return {"key": key, "status": "error", "error": err}


t0 = time.time()
journal = []
with ThreadPoolExecutor(max_workers=8) as ex:
    for i, rec in enumerate(ex.map(fetch, keys), 1):
        journal.append(rec)
        if i % 25 == 0 or rec["status"] == "error":
            print(f"{i}/{len(keys)} {rec}", flush=True)
counts = {}
for r in journal:
    counts[str(r["status"])] = counts.get(str(r["status"]), 0) + 1
print("status counts:", counts, f"{time.time() - t0:.0f}s", flush=True)
with open(os.path.join(STG, "fetch_missing_chunks_journal.json"), "w") as f:
    json.dump({"ref_base": REF_BASE, "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "counts": counts, "journal": journal}, f, indent=1)

linked = 0
for r in journal:
    if r.get("bytes") != CHUNK_BYTES:
        continue
    z, y, x = r["key"].split(".")
    dst = os.path.join(MIRROR, r["key"])
    if not os.path.exists(dst):
        os.link(os.path.join(CACHE, f"chunk_{z}_{y}_{x}.bin"), dst)
        linked += 1
print(f"linked {linked} new chunks into mirror; mirror now has",
      sum(1 for k in os.listdir(MIRROR) if k.count(".") == 2), "chunk files")
