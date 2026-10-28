# -*- coding: utf-8 -*-
"""第 6b 步：Ga-As 首峰左侧「台阶/肩」的全温度对照 + 统计显著性

只用已经算好的细网格数据（data/<T>K/rdf_dr0.01_avg*.dat），不重新取帧。

g(r) 的 Poisson 误差可直接由 g 本身给出，不需要逐 bin 计数：
    sigma_g = sqrt( g * K / phi ) / sqrt(n),   K = V/(N_A N_B), phi = 4*pi*r^2*dr
推导：hbar = g*phi/K，sigma_hbar = sqrt(hbar)/sqrt(n)。
"""
import os, re, glob, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

OUT = r"E:\论文\GaAs纳米液滴结晶\reCalculate"
K = 1000000.0 / (1990 * 2026)          # V/(N_As*N_Ga)，4016 原子，盒 100^3
DR = 0.01
TEMPS = [200, 400, 700, 1000, 1800]
COL = {200: "k", 400: "red", 700: "blue", 1000: "purple", 1800: "#DAA520"}


def load(T):
    d = os.path.join(OUT, "data", "%dK" % T)
    f = sorted(glob.glob(os.path.join(d, "rdf_dr0.01_avg*.dat")))[-1]
    n = int(re.search(r"avg(\d+)", f).group(1))
    a = np.array([[float(v) for v in l.split()] for l in open(f) if not l.startswith("#")])
    return a, n


def sigma(g, r, n, dr=DR):
    phi = 4 * np.pi * r ** 2 * dr
    return np.sqrt(np.abs(g) * K / phi) / np.sqrt(n)


def regrid(r, g, fac):
    """把 dr=0.01 的曲线体积加权平均到 dr=fac*0.01"""
    phi = 4 * np.pi * r ** 2 * DR
    nb = (len(g) // fac) * fac
    G = g[:nb].reshape(-1, fac); P = phi[:nb].reshape(-1, fac); R = r[:nb].reshape(-1, fac)
    gc = (G * P).sum(1) / P.sum(1)
    rc = (R ** 3).sum(1) / (R ** 2).sum(1)
    return rc, gc


def analyze(r, g, sg, lo, hi):
    """返回 (局部极大列表, 台阶位置, 台阶深度, 深度噪声, 显著性)
    台阶深度 = min(左斜率, 右斜率) - 本段斜率，取窗口内最深者"""
    m = (r >= lo) & (r <= hi)
    rr, gg, ss = r[m], g[m], sg[m]
    d = np.diff(gg)
    loc = [(rr[i + 1], gg[i + 1], ss[i + 1]) for i in range(1, len(d)) if d[i - 1] > 0 > d[i]]
    best = (np.nan, -np.inf, np.nan, np.nan)
    for k in range(1, len(d) - 1):
        depth = min(d[k - 1], d[k + 1]) - d[k]
        sd = np.sqrt(ss[k - 1] ** 2 + 2 * ss[k] ** 2 + ss[k + 1] ** 2)
        if depth > best[1]:
            best = (rr[k], depth, sd, depth / sd if sd > 0 else np.nan)
    return (loc,) + best


print("=" * 100)
print("K = V/(N_As*N_Ga) = %.5f    dr = %.2f A    窗口 2.30-2.50 A" % (K, DR))
print("=" * 100)
print("\n[全温度] Ga-As 上升沿（2.30-2.50 A）中的局部极大与斜率台阶")
print("%6s %4s | %-26s | %-42s" % ("T/K", "n", "局部极大（=严格意义的肩峰）", "斜率最小处（=视觉上的台阶）"))
print("-" * 100)

rows, curves = [], {}
for T in TEMPS:
    a, n = load(T)
    r, g = a[:, 0], a[:, 2]
    curves[T] = (r, g, n)
    sg = sigma(g, r, n)
    loc, rk, depth, sdep, nsig = analyze(r, g, sg, 2.30, 2.50)
    mx = ", ".join("%.3f (g=%.1f)" % (x, y) for x, y, _ in loc) or "无"
    if not np.isfinite(depth):
        stp = "窗口内无台阶"
    else:
        m = (r >= 2.30) & (r <= 2.50)
        rr, gg = r[m], g[m]
        k = int(np.argmin(np.abs(rr - rk)))
        stp = ("r=%.3f Å | 左斜率 %.0f，台阶处 %.0f，右斜率 %.0f /Å | 深度 %.1f±%.1f (%.1fσ)"
               % (rk, (gg[k] - gg[k - 1]) / DR, (gg[k + 1] - gg[k]) / DR,
                  (gg[k + 2] - gg[k + 1]) / DR, depth, sdep, nsig))
    rows.append((T, n, len(loc), round(float(rk), 3) if np.isfinite(rk) else "",
                 round(float(depth), 1) if np.isfinite(depth) else "",
                 round(float(nsig), 2) if np.isfinite(nsig) else ""))
    print("%6d %4d | %-26s | %-42s" % (T, n, mx, stp))

with open(os.path.join(OUT, "data", "shoulder_summary.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["T_nominal_K", "n_frames", "n_local_max_in_2.30_2.50A", "step_r_A",
                "step_depth_perA", "step_significance_sigma"])
    w.writerows(rows)
print("\n-> 写出 data/shoulder_summary.csv")

# ---------------------------------------------------------------- 网格粗化
print("\n[网格粗化测试] 200 K Ga-As：台阶随网格变粗如何演化")
r, g, n = curves[200]
for fac in (1, 2, 5, 10):
    rc, gc = regrid(r, g, fac)
    sg = sigma(gc, rc, n, DR * fac)
    loc, rk, depth, sdep, nsig = analyze(rc, gc, sg, 2.30, 2.50)
    mm = (rc >= 2.30) & (rc <= 2.50)
    d = np.diff(gc[mm]) / (DR * fac)
    print("   dr=%5.3f Å : 台阶 %.3f Å | 窗口内斜率 %.0f → %.0f → ... → %.0f /Å | 深度 %.1f /Å (%.1fσ) | 局部极大 %d 个"
          % (DR * fac, rk, d[0], d.min(), d.max(), depth, nsig, len(loc)))

# ---------------------------------------------------------------- 图
rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False,
                 "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9,
                 "legend.fontsize": 8.5, "axes.linewidth": 0.9,
                 "xtick.direction": "in", "ytick.direction": "in",
                 "xtick.top": True, "ytick.right": True, "pdf.fonttype": 42})
fig, axs = plt.subplots(2, 2, figsize=(11.6, 8.5), constrained_layout=True)

# (a)
ax = axs[0, 0]
r, g, n = curves[200]
sg = sigma(g, r, n)
m = (r >= 2.15) & (r <= 2.62)
ax.fill_between(r[m], (g - sg)[m], (g + sg)[m], color="0.8", lw=0,
                label="±1σ（Poisson，%d 帧平均）" % n)
ax.plot(r[m], g[m], "k-", lw=1.9, label="g(r)")
ax.axvspan(2.433, 2.452, color="orange", alpha=.22, lw=0)
ax.annotate("斜率台阶 2.44 Å\n（拐点，g 仍在涨）", xy=(2.442, 96), xytext=(2.50, 34),
            arrowprops=dict(arrowstyle="->", color="darkorange", lw=1.0),
            color="darkorange", fontsize=9)
ax.annotate("真主峰 2.505 Å", xy=(2.505, 145), xytext=(2.55, 130),
            arrowprops=dict(arrowstyle="->", color="crimson", lw=1.0),
            color="crimson", fontsize=9)
ax.set_xlim(2.15, 2.65); ax.set_ylim(0, 165)
ax.set_xlabel("Pair separation distance (Å)"); ax.set_ylabel("g(r)")
ax.set_title("(a) 200 K Ga–As：拐点有，第二极大没有", fontsize=10)
ax.legend(frameon=False, loc="upper left")

# (b)
ax = axs[0, 1]
for T in TEMPS:
    r, g, n = curves[T]
    m = (r >= 2.2) & (r <= 2.62)
    ax.plot(r[m], g[m], lw=1.7 if T in (200, 400, 700) else 1.1, color=COL[T], label="%d K" % T)
ax.axvspan(2.433, 2.452, color="orange", alpha=.15, lw=0)
ax.axvline(2.448, color="0.35", ls=":", lw=1.0)
ax.text(2.455, 118, "理想 ZB 键长 2.448 Å", color="0.3", fontsize=8)
ax.set_xlim(2.2, 2.62); ax.set_ylim(0, 165)
ax.set_xlabel("Pair separation distance (Å)"); ax.set_ylabel("g(r)")
ax.set_title("(b) 固态各温度台阶都在同一位置", fontsize=10)
ax.legend(frameon=False, loc="upper left")

# (c)
ax = axs[1, 0]
r, g, n = curves[200]
for fac, lab, cc in ((1, "dr = 0.01 Å", "k"), (2, "dr = 0.02 Å", "red"),
                     (5, "dr = 0.05 Å", "blue"), (10, "dr = 0.10 Å（原作者）", "#c44e52")):
    rc, gc = regrid(r, g, fac)
    m = (rc >= 2.2) & (rc <= 2.62)
    ax.plot(rc[m], gc[m], "-", color=cc, lw=1.6, marker="o" if fac >= 5 else None, ms=3.2, label=lab)
ax.set_xlim(2.2, 2.62); ax.set_ylim(0, 165)
ax.set_xlabel("Pair separation distance (Å)"); ax.set_ylabel("g(r)")
ax.set_title("(c) 粗网格把台阶抹成一个平台", fontsize=10)
ax.legend(frameon=False, loc="upper left")

# (d)
ax = axs[1, 1]
r, g, n = curves[200]
phi = 4 * np.pi * r ** 2 * DR
h = g * phi / K
m = (r >= 2.2) & (r <= 2.62)
ax.semilogy(r[m], h[m] / n, "k-", lw=1.5, label="单帧每 bin 计数")
ax.semilogy(r[m], h[m], "--", color="crimson", lw=1.3, label="%d 帧合计计数" % n)
ax.set_xlim(2.2, 2.62)
ax.set_xlabel("Pair separation distance (Å)"); ax.set_ylabel("Ga–As 对数（0.01 Å 壳内）")
ax.set_title("(d) 分辨率尺子：薄壳里的原子对够不够多", fontsize=10)
ax.legend(frameon=False, loc="upper left")

fig.suptitle("Ga–As 首峰左侧「肩峰」检验：是斜率台阶，不是第二个极大", fontsize=12.5)
p = os.path.join(OUT, "figures", "fig12_shoulder_all_temps.png")
fig.savefig(p, dpi=170)
print("-> 写出 %s" % p)
