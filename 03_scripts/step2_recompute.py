# -*- coding: utf-8 -*-
"""第 2 步：按温度取帧 -> 归档 -> 算细网格 RDF
每帧先用 dr=0.1 与原作者 RDF/<T> 比对，确认取到的是同一帧，再算 dr=0.01 / 0.02
"""
import os, sys, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rdf_engine import read_dump_frame, extract_frame, compute_rdf, write_rdf

BASE = r"E:\论文\GaAs纳米液滴结晶\Semi-GaAs-droplet\cooling\NVT-2000-200-8ns"
BIG  = os.path.join(BASE, "DzhbGaAs-DeepMDwork1Semi-GaAs-dropletcoolingNVT-2000-200-8nsGaAs-200k.dump")
RDFDIR = os.path.join(BASE, "RDF")
OUT = r"E:\论文\GaAs纳米液滴结晶\reCalculate"

# T(t) = 2000 - 225*t[ns],  t[ns] = step*1e-6  (timestep 0.001 ps)
TEMPS = [1800, 1000, 800, 700, 400, 200]
KEYS = ["As-As", "As-Ga", "Ga-Ga"]

def ref_rdf(T):
    p = os.path.join(RDFDIR, "%dK" % T)
    return np.array([[float(v) for v in l.split()] for l in open(p, errors="ignore")
                     if not l.startswith("#") and len(l.split()) == 4])

def est_step(T):
    s = (2000 - T) / 225.0 * 1e6
    return int(round(s / 50.0) * 50)

os.makedirs(os.path.join(OUT, "frames"), exist_ok=True)
report = []
import csv as _csv
_fi = {}
_fp = os.path.join(OUT, "frames", "frame_index.csv")
if os.path.exists(_fp):
    for row in _csv.DictReader(open(_fp)):
        _fi[int(row["T_nominal_K"])] = int(row["step"])

for T in TEMPS:
    ref = ref_rdf(T)
    if T in _fi:
        s, box, typ, xyz = extract_frame(BIG, _fi[T])
        res = compute_rdf(box, typ, xyz, rmax=10.0, dr=0.1)
        err = max(np.abs(res[k][1] - ref[:, j + 1]).max() for j, k in enumerate(KEYS))
    else:
        err, s, box, typ, xyz = 9.99, *((lambda q: (q[1], q[2], q[3], q[4]))(extract_frame(BIG, est_step(T))))
    Tactual = 2000 - 225 * s * 1e-6
    print("T=%4d K  用 step=%d (T_actual=%.1f K)  与原作者 RDF 的 max|Δ| = %.4f %s"
          % (T, s, Tactual, err, "OK" if err < 0.05 else "<<< 存疑"))
    # 归档该帧
    frm = os.path.join(OUT, "frames", "%dK_step%d.dump" % (T, s))
    with open(frm, "w") as g:
        g.write("ITEM: TIMESTEP\n%d\n" % s)
        g.write("ITEM: NUMBER OF ATOMS\n%d\n" % len(typ))
        g.write("ITEM: BOX BOUNDS pp pp pp\n")
        for k in range(3):
            g.write("%.8g %.8g\n" % (0.0, box[k]))
        g.write("ITEM: ATOMS type x y z id\n")
        for k in range(len(typ)):
            g.write("%d %.10g %.10g %.10g %d\n" % (typ[k], xyz[k,0], xyz[k,1], xyz[k,2], k + 1))
    # 细网格
    d = os.path.join(OUT, "data", "%dK" % T)
    os.makedirs(d, exist_ok=True)
    for dr in (0.1, 0.02, 0.01):
        res = compute_rdf(box, typ, xyz, rmax=10.0, dr=dr)
        write_rdf(os.path.join(d, "rdf_dr%.3f.dat" % dr), res,
                  header_extra="T=%dK (target), step=%d, T_actual=%.2f K, single frame, dr=%.3f A"
                               % (T, s, Tactual, dr))
    # 首峰区域单独导一份 CSV（方便直接看/画）
    res01 = compute_rdf(box, typ, xyz, rmax=6.0, dr=0.01)
    r = res01[KEYS[0]][0]
    with open(os.path.join(d, "first_peak_dr0.01.csv"), "w") as g:
        g.write("r_A,As-As,As-Ga,Ga-Ga\n")
        for i in range(len(r)):
            g.write("%.4f,%s\n" % (r[i], ",".join("%.6f" % res01[k][1][i] for k in KEYS)))
    report.append((T, s, Tactual, err))

print("\n完成。汇总：")
for T, s, Ta, e in report:
    print("  %4d K  step=%-9d T_actual=%7.2f K   max|Δ| vs 原作者 = %.4f" % (T, s, Ta, e))
