"""Cross-scan overlap analysis for villa #192 dual-energy validation route.

Metadata-only: consumes the three official transform.json files downloaded to
data_crossscan_meta_2026-08-07/ and the known zarr level-0 shapes. No volume
bodies are read. Convention verified by hand on the first 137keV landmark pair:
x_fixed = M @ x_moving + t, landmark/matrix axis order is (x, y, z).
"""
import json
import numpy as np

META = "data_crossscan_meta_2026-08-07"

# zarr level-0 shapes, (z, y, x) as stored
SHAPES_ZYX = {
    "137keV": (6625, 8431, 8431),
    "78keV": (75784, 32693, 32693),
    "1.129um": (59969, 36006, 32354),
}

# phase4b prediction 256^3 chunk coords (z, y, x) from the frozen 2026-08-06 snapshot
PHASE4B_PRED_CHUNKS_ZYX = {
    "P25_INTERIOR": (73, 43, 40),
    "P50_INTERIOR": (148, 58, 75),
    "P75_INTERIOR": (222, 89, 80),
}


def load(tag):
    with open(f"{META}/transform_{tag}.json") as f:
        d = json.load(f)
    M = np.array(d["transformation_matrix"], dtype=float)
    A, t = M[:, :3], M[:, 3]
    fixed = np.array(d["fixed_landmarks"], dtype=float)
    moving = np.array(d["moving_landmarks"], dtype=float)
    return d["fixed_volume"], A, t, fixed, moving


def residuals(A, t, fixed, moving):
    pred = moving @ A.T + t
    r = np.linalg.norm(pred - fixed, axis=1)
    return r


def corners_xyz(shape_zyx):
    z, y, x = shape_zyx
    return np.array([[xi, yi, zi] for zi in (0, z - 1) for yi in (0, y - 1) for xi in (0, x - 1)], dtype=float)


def main():
    out = {}
    fv137, A137, t137, f137, m137 = load("20260323153942")
    fv78, A78, t78, f78, m78 = load("20260411134726")
    fv1129, A1129, t1129, f1129, m1129 = load("20260608103018")

    print(f"137keV  fixed_volume: {fv137}")
    print(f"78keV   fixed_volume: {fv78}")
    print(f"1.129um fixed_volume: {fv1129}")

    for name, (A, t, f, m, unit) in {
        "137keV->fixed(7.91um?)": (A137, t137, f137, m137, None),
        "78keV->fixed(7.91um?)": (A78, t78, f78, m78, None),
        "1.129um->fixed(2.4um)": (A1129, t1129, f1129, m1129, 2.4),
    }.items():
        r = residuals(A, t, f, m)
        scale = np.linalg.norm(A, axis=0).mean()
        print(f"\n{name}: n={len(r)} scale~{scale:.5f} "
              f"residual fixed-voxels mean={r.mean():.2f} max={r.max():.2f}")
        if unit:
            print(f"  in um (fixed voxel={unit}um): mean={r.mean()*unit:.2f}um max={r.max()*unit:.2f}um")

    # implied fixed voxel size if nominal moving voxel is exact
    s137 = np.linalg.norm(A137, axis=0).mean()
    s78 = np.linalg.norm(A78, axis=0).mean()
    print(f"\nimplied fixed voxel um via 137keV: {2.4/s137:.4f}; via 78keV: {2.4/s78:.4f}")

    # --- chain 78keV -> fixed -> 137keV ---
    A137i = np.linalg.inv(A137)

    def to137(p78_xyz):
        p = np.atleast_2d(np.asarray(p78_xyz, dtype=float))
        out = ((p @ A78.T) + t78 - t137) @ A137i.T
        return out[0] if np.asarray(p78_xyz).ndim == 1 else out

    print("\n== phase4b cube centers mapped to 137keV voxel space (x,y,z) ==")
    lim137 = np.array([SHAPES_ZYX["137keV"][2], SHAPES_ZYX["137keV"][1], SHAPES_ZYX["137keV"][0]], float)
    for name, (cz, cy, cx) in PHASE4B_PRED_CHUNKS_ZYX.items():
        c78 = np.array([cx * 256 + 128, cy * 256 + 128, cz * 256 + 128], float)
        p = to137(c78)
        inside = bool(np.all(p >= 0) and np.all(p < lim137))
        print(f"  {name}: 78keV(x,y,z)={c78} -> 137keV={np.round(p,1)} inside={inside}")

    # --- 137keV coverage band expressed in 78keV voxel space ---
    A78i = np.linalg.inv(A78)

    def to78(p137_xyz):
        return (A78i @ ((A137 @ p137_xyz.T).T + t137 - t78).T).T

    c137 = corners_xyz(SHAPES_ZYX["137keV"])
    c78cov = to78(c137)
    lo, hi = c78cov.min(axis=0), c78cov.max(axis=0)
    print("\n== 137keV volume corners mapped into 78keV voxel space (x,y,z) ==")
    print(f"  bounding box lo={np.round(lo,1)} hi={np.round(hi,1)}")
    print(f"  => 78keV z band ~[{lo[2]:.0f}, {hi[2]:.0f}] "
          f"(256-chunk z index [{lo[2]/256:.1f}, {hi[2]/256:.1f}])")
    out["band_z_chunks"] = (lo[2] / 256, hi[2] / 256)

    # --- 1.129um coverage in 2.4um space (direct transform) ---
    c1129 = corners_xyz(SHAPES_ZYX["1.129um"])
    cov = (c1129 @ A1129.T) + t1129
    lo2, hi2 = cov.min(axis=0), cov.max(axis=0)
    print("\n== 1.129um volume corners mapped into 2.4um fixed space (x,y,z) ==")
    print(f"  bounding box lo={np.round(lo2,1)} hi={np.round(hi2,1)}")
    A1129i = np.linalg.inv(A1129)
    lim1129 = np.array([SHAPES_ZYX["1.129um"][2], SHAPES_ZYX["1.129um"][1], SHAPES_ZYX["1.129um"][0]], float)
    for name, (cz, cy, cx) in PHASE4B_PRED_CHUNKS_ZYX.items():
        c78 = np.array([cx * 256 + 128, cy * 256 + 128, cz * 256 + 128], float)
        p = A1129i @ (c78 - t1129)
        inside = bool(np.all(p >= 0) and np.all(p < lim1129))
        print(f"  {name} center {c78} -> 1.129um voxel {np.round(p,1)} inside-volume={inside}")

    # --- count prediction 256^3 chunks fully inside the 137keV volume (margin voxels) ---
    margin = 32.0
    zc, yc, xc = (np.array(SHAPES_ZYX["78keV"]) // 256)
    z_range = range(max(0, int(lo[2] // 256)), min(int(zc), int(hi[2] // 256) + 1))
    y_range = range(max(0, int(lo[1] // 256)), min(int(yc), int(hi[1] // 256) + 1))
    x_range = range(max(0, int(lo[0] // 256)), min(int(xc), int(hi[0] // 256) + 1))
    inside_chunks = []
    for zi in z_range:
        for yi in y_range:
            for xi in x_range:
                # all 8 corners of the 256^3 chunk must map inside 137keV with margin
                cs = np.array([[xi * 256 + dx, yi * 256 + dy, zi * 256 + dz]
                               for dz in (0, 256) for dy in (0, 256) for dx in (0, 256)], float)
                p = to137(cs)
                if np.all(p >= margin) and np.all(p < lim137 - margin):
                    inside_chunks.append((zi, yi, xi))
    print(f"\n== prediction 256^3 chunks fully inside 137keV (margin {margin:g} vox) ==")
    print(f"  count={len(inside_chunks)}")
    if inside_chunks:
        arr = np.array(inside_chunks)
        print(f"  z-chunk range [{arr[:,0].min()},{arr[:,0].max()}], "
              f"y [{arr[:,1].min()},{arr[:,1].max()}], x [{arr[:,2].min()},{arr[:,2].max()}]")


if __name__ == "__main__":
    main()
