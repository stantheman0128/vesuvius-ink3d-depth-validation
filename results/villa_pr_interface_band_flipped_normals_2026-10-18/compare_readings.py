"""Point-by-point comparison of the original and negated-normal runs.
Usage: python compare_readings.py  (run inside this directory)"""
import csv, math, statistics as st
def load(p): return {r["id"]: r for r in csv.DictReader(open(p, newline="", encoding="utf-8"))}
def fl(x):
    try: return float(x)
    except ValueError: return None
for tag in ("sigma0", "sigma2"):
    O, F = load(f"w00_points_{tag}.csv"), load(f"w00_points_{tag}_normals_negated.csv")
    both = [i for i in O if i in F and O[i]["evaluable"] == "True" and F[i]["evaluable"] == "True"]
    P = lambda r: [fl(r[k]) for k in ("point_z", "point_y", "point_x")]
    disp = [math.dist(P(O[i]), P(F[i])) for i in both]
    dD = [fl(F[i]["D_voxels"]) - fl(O[i]["D_voxels"]) for i in both]
    same = sum((fl(O[i]["D_voxels"]) <= 3) == (fl(F[i]["D_voxels"]) <= 3) for i in both)
    ib = [(O[i]["inside_brighter"], F[i]["inside_brighter"]) for i in both
          if O[i]["inside_brighter"] in ("True", "False") and F[i]["inside_brighter"] in ("True", "False")]
    print(f"[{tag}] evaluable in both runs: {len(both)} (rows with evaluable == True in both files)")
    print(f"  point displacement in reference voxels: median {st.median(disp):.2f}, min {min(disp):.2f}, max {max(disp):.2f}")
    print(f"  within 3 voxels on the common set: original {sum(fl(O[i]['D_voxels'])<=3 for i in both)}, negated {sum(fl(F[i]['D_voxels'])<=3 for i in both)}; same verdict {same}/{len(both)}")
    print(f"  D unchanged {sum(abs(x)<1e-9 for x in dD)}, farther {sum(x>0 for x in dD)}, closer {sum(x<0 for x in dD)}, median |dD| {st.median(abs(x) for x in dD):.2f}")
    print(f"  inside_brighter defined in both: {len(ib)}; same {sum(a==b for a,b in ib)}, opposite {sum(a!=b for a,b in ib)}")
