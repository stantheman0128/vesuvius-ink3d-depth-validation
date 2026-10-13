"""Print the summary block of each w00 run, old (_pre-fetch) beside new when present."""
import json
import os

STG = os.path.dirname(os.path.abspath(__file__))
for name in ("w00_interface_distance", "w00_interface_distance_sigma0",
             "w00_interface_distance_volnormal"):
    for suffix in ("_pre-fetch", ""):
        path = os.path.join(STG, f"{name}{suffix}.json")
        if not os.path.exists(path):
            continue
        s = json.load(open(path))["summary"]
        print(f"{name}{suffix or ' (new)'}: evaluable {s['n_evaluable']}, "
              f"within {s['n_within_threshold']} ({s['fraction_within_threshold']}), "
              f"wilson {s['wilson95_fraction_within']}, median {s['median_D_voxels']}, "
              f"reasons {s['non_evaluable_reasons']}")
