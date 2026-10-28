# -*- coding: utf-8 -*-
"""第 3 步（修订）：出图 + 统计误差带
  fig1: 三个 pair 的 g(r)（6 个温度，dr=0.01）
  fig2: 首峰放大（200/400/700/800 K）
  fig3: 200 K 粗网格 vs 细网格（细网格带 Poisson 误差带）
  fig4: 200 K 首峰 —— 用 ±3 帧做时间平均，看噪声降下去后还剩什么
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rdf_engine import extract_frame, compute_rdf

BASE = r"E:\论文\GaAs纳米液滴结晶\Semi-GaAs-droplet\cooling\NVT-2000-200-8ns"
BIG  = os.path.join(BASE, "DzhbGaAs-DeepMDwork1Semi-GaAs-dropletcoolingNVT-2000-200-8nsGaAs-200k.dump")
OUT = r"E:\论文\GaAs纳米液滴结晶\reCalculate"
KEYS = ["As-As", "As-Ga", "Ga-Ga"]
TEMPS = [1800, 1000, 800, 700, 400, 200]
COL = {1800: "#DAA520", 1000: "purple", 800: "green", 700: "blue", 400: "red", 200: "black"}

rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False,
                 "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
                 "legend.fontsize": 9, "axes.linewidth": 0.9,
                 "xtick.direction": "in", "ytick.direction": "in",
                 "xtick.top": True, "ytick.right": True, "pdf.fonttype": 42})

def load(T, dr):
    p = os.path.join(OUT, "data", "%dK" % T, "rdf_dr%.3f.dat" % dr)
    a = np.array([[float(v) for v in l.split()] for l in open(p) if not l.startswith("#")])
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]

CFG = [(1, "(a) Ga–Ga", (3.2, 5.2)), (2, "(b) Ga–As", (2.1, 3.1)), (3, "(c) As–As", (3.2, 5.2))]

# ---------------- fig1 ----------------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
for ax, (col, name, _) in zip(axs, CFG):
    for T in TEMPS:
        r, a1, a2, a3 = load(T, 0.01)
        ax.plot(r, [a1, a2, a3][col - 1], lw=1.0, color=COL[T], label="%d K" % T)
    ax.set_xlim(0, 10); ax.set_ylim(0, None)
    ax.set_xlabel("Pair separation distance (Å)"); ax.set_ylabel("g(r)")
    ax.set_title(name, fontsize=11)
axs[0].legend(loc="upper right", frameon=False)
fig.suptitle("重算 g(r)，dr = 0.01 Å（单帧，与原图同一帧）", fontsize=12)
fig.savefig(os.path.join(OUT, "figures", "fig1_rdf_three_pairs.png"), dpi=150)

# ---------------- fig2 ----------------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
for ax, (col, name, xr) in zip(axs, CFG):
    for T in (200, 400, 700, 800):
        r, a1, a2, a3 = load(T, 0.01)
        y = [a1, a2, a3][col - 1]; m = (r >= xr[0]) & (r <= xr[1])
        ax.plot(r[m], y[m], lw=1.1, color=COL[T], label="%d K" % T)
    ax.set_xlim(*xr); ax.set_xlabel("r (Å)"); ax.set_ylabel("g(r)")
    ax.set_title(name, fontsize=11); ax.legend(frameon=False)
fig.suptitle("首峰区域放大（dr = 0.01 Å）", fontsize=12)
fig.savefig(os.path.join(OUT, "figures", "fig2_first_peak_zoom.png"), dpi=150)

# ---------------- fig3: 200K 粗 vs 细 + 误差带 ----------------
s, box, typ, xyz = extract_frame(BIG, 8000000)
fine = compute_rdf(box, typ, xyz, rmax=10.0, dr=0.01, return_counts=True)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
for ax, (col, name, xr) in zip(axs, CFG):
    r1, b1, b2, b3 = load(200, 0.1)
    y1 = [b1, b2, b3][col - 1]; m1 = (r1 >= xr[0]) & (r1 <= xr[1])
    rf, g, h = fine[KEYS[col - 1]]
    err = g / np.sqrt(np.maximum(h, 1))
    m2 = (rf >= xr[0]) & (rf <= xr[1])
    ax.fill_between(rf[m2], (g - err)[m2], (g + err)[m2], color="crimson", alpha=0.22,
                    lw=0, label="±1σ (Poisson)")
    ax.plot(r1[m1], y1[m1], "o-", color="0.45", lw=1.0, ms=6, label="原网格 dr=0.1 Å（100 点）")
    ax.plot(rf[m2], g[m2], "-", color="crimson", lw=1.0, label="细网格 dr=0.01 Å（1000 点）")
    ax.set_xlim(*xr); ax.set_xlabel("r (Å)"); ax.set_ylabel("g(r)")
    ax.set_title(name, fontsize=11); ax.legend(frameon=False, fontsize=8.5)
fig.suptitle("200 K 首峰：粗网格 vs 细网格（细网格的起伏基本落在单帧统计噪声内）", fontsize=12)
fig.savefig(os.path.join(OUT, "figures", "fig3_200K_coarse_vs_fine.png"), dpi=150)

# ---------------- fig4: 200K 用 ±3 帧平均降噪 ----------------
acc = {k: None for k in KEYS}
n = 0
for ds in (-150, -100, -50, 0, 50, 100, 150):
    st = 8000000 + ds
    if st < 0 or st > 8000000:                          # 越过末帧的直接跳过（不做重复计数）
        continue
    _, bx, ty, xz = extract_frame(BIG, st)
    res = compute_rdf(bx, ty, xz, rmax=10.0, dr=0.01, return_counts=True)
    for k in KEYS:
        r_, g_, h_ = res[k]
        acc[k] = h_ if acc[k] is None else acc[k] + h_
    n += 1
def gm_of(k):
    a, b = {"As-As": (1, 1), "As-Ga": (1, 2), "Ga-Ga": (2, 2)}[k]
    return (acc[k] / n) * (V / (idx[a] * idx[b])) / phi

print("200 K 平均用了 %d 帧（step %d..%d，T 跨 %.1f K）"
      % (n, 8000000 - 150, 8000000, 225 * 150 * 1e-6))

fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
for ax, (col, name, xr) in zip(axs, CFG):
    h = acc[KEYS[col - 1]]
    V = float(np.prod(box))
    idx = {1: (typ == 1).sum(), 2: (typ == 2).sum()}
    a, b = {1: (1, 1), 2: (1, 2), 3: (2, 2)}[col]
    r_, g = fine[KEYS[col - 1]][0], None
    phi = 4 * np.pi * r_ ** 2 * 0.01
    gm = gm_of(KEYS[col - 1])          # 已除以帧数 n
    g1 = fine[KEYS[col - 1]][1]
    errm = gm / np.sqrt(np.maximum(h, 1)) / np.sqrt(n)
    m = (r_ >= xr[0]) & (r_ <= xr[1])
    ax.fill_between(r_[m], (gm - errm)[m], (gm + errm)[m], color="crimson", alpha=0.25, lw=0)
    ax.plot(r_[m], g1[m], "-", color="0.65", lw=0.9, label="单帧 (step 8000000)")
    ax.plot(r_[m], gm[m], "-", color="crimson", lw=1.5, label="%d 帧平均" % n)
    ax.set_xlim(*xr); ax.set_xlabel("r (Å)"); ax.set_ylabel("g(r)")
    ax.set_title(name, fontsize=11); ax.legend(frameon=False, fontsize=9)
fig.suptitle("200 K 首峰：单帧 vs 时间平均（平均后噪声降到 ~1/√%d）" % n, fontsize=12)
fig.savefig(os.path.join(OUT, "figures", "fig4_200K_timeaveraged.png"), dpi=150)

# 写出平均后的数据
d = os.path.join(OUT, "data", "200K_avg4frames")
os.makedirs(d, exist_ok=True)
with open(os.path.join(d, "rdf_dr0.01_avg.dat"), "w") as f:
    f.write("# 200 K, average of %d frames (step %d..%d), dr=0.01 A\n" % (n, 8000000 - 150, 8000000))
    f.write('# "Pair separation distance" ' + " ".join(KEYS) + "\n")
    for i in range(len(r_)):
        f.write("%.4f " % r_[i] + " ".join("%.10g" % (acc[k][i] * (V / (idx[{1: (1,1), 2: (1,2), 3: (2,2)}[[("As-As").find("x")][0]] if False else 1, 1][0]) / phi[i])) for k in [] if False) + " \n" if False else "")
# 上面写法太绕，改成清晰版
with open(os.path.join(d, "rdf_dr0.01_avg.dat"), "w") as f:
    f.write("# 200 K, average of %d frames (step %d..%d), dr=0.01 A, single-frame identical geometry\n"
            % (n, 8000000 - 150, 8000000))
    f.write('# "Pair separation distance" ' + " ".join(KEYS) + "\n")
    for i in range(len(r_)):
        f.write("%.4f " % r_[i] + " ".join("%.10g" % gm_of(k)[i] for k in KEYS) + " \n")
print("figures saved; 平均数据 -> data/200K_avg4frames/rdf_dr0.01_avg.dat")

from scipy.signal import find_peaks
print("\n=== 200 K 首峰：局部极大（prominence > 2% 峰高） ===")
for col, name, xr in CFG:
    for tag, yy_ in (("单帧   ", fine[KEYS[col - 1]][1]), ("7帧平均", gm_of(KEYS[col - 1]))):
        rr = fine[KEYS[col - 1]][0]
        m = (rr >= xr[0]) & (rr <= xr[1]); rr2, yy2 = rr[m], yy_ [m]
        pk, prop = find_peaks(yy2, prominence=0.02 * yy2.max())
        print("  %-12s %s 峰: %s" % (name, tag,
              "  ".join("%.2f(%.1f, prom %.1f)" % (rr2[i], yy2[i], p)
                        for i, p in zip(pk, prop["prominences"])) or "无"))
        # 主峰两侧 ±0.25 Å 的极小值（判断是否真的分开）
        if len(pk) >= 1:
            k0 = pk[np.argmax(yy2[pk])]
            lo, hi = max(0, k0 - 25), min(len(yy2) - 1, k0 + 25)
            print("       主峰 r=%.2f g=%.1f ；±0.25 Å 内最小 g=%.1f（相对落差 %.1f%%）"
                  % (rr2[k0], yy2[k0], yy2[lo:hi+1].min(),
                     100 * (yy2[k0] - yy2[lo:hi+1].min()) / yy2[k0]))
