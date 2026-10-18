"""Write a copy of a depth_anchors CSV with normal_x, normal_y, normal_z negated.
Usage: python negate_normals.py in.csv out.csv"""
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, newline="", encoding="utf-8") as f, open(dst, "w", newline="", encoding="utf-8") as g:
    r = csv.reader(f); w = csv.writer(g); h = next(r); w.writerow(h)
    idx = [h.index(c) for c in ("normal_x", "normal_y", "normal_z")]
    for row in r:
        for i in idx:
            v = row[i].strip()
            row[i] = v[1:] if v.startswith("-") else "-" + v
        w.writerow(row)
