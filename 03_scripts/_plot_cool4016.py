# -*- coding: utf-8 -*-
"""4016 原子 GaAs 纳米液滴：225 vs 900 K/ns 降温结晶对比
版式与 cryst_vs_rate_alt100new_v2 一致（14.65 cm 宽，1:1 插入 SI/正文）
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SRC = r"E:\论文\GaAs纳米液滴结晶\Semi-GaAs-droplet\cooling\rate_compare_4016"

FS_LAB, FS_TICK, FS_LEG, FS_PANEL = 10.0, 9.5, 9.0, 10.5

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": FS_LAB, "axes.labelsize": FS_LAB, "axes.labelpad": 3.2,
    "axes.linewidth": 0.9,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "xtick.major.size": 3.6, "ytick.major.size": 3.6,
    "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.minor.width": 0.7, "ytick.minor.width": 0.7,
    "xtick.major.pad": 3.0, "ytick.major.pad": 3.0,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "legend.frameon": False, "legend.fontsize": FS_LEG,
    "legend.handlelength": 2.6, "legend.handletextpad": 0.7,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": False,
    "savefig.facecolor": "white"})


def smooth(y, w=5):
    if len(y) < w:
        return y
    yp = np.pad(y, (w // 2, w // 2), mode="edge")
    return np.convolve(yp, np.ones(w) / w, mode="valid")


SERIES = [("run8ns", 225, "#1F4E79", "-",            "225 K/ns"),
          ("run2ns", 900, "#A50F15", (0, (4.2, 1.5)), "900 K/ns")]

W_CM, H_CM = 14.65, 6.80
fig = plt.figure(figsize=(W_CM / 2.54, H_CM / 2.54))
gs = fig.add_gridspec(1, 2, width_ratios=[2.05, 1.0], wspace=0.32,
                      left=0.085, right=0.978, top=0.955, bottom=0.185)
axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

fin = {}
for tag, rate, col, ls, lab in SERIES:
    d = np.load(os.path.join(SRC, f"cool4016_{tag}.npz"))
    axA.plot(d["T"], smooth(d["zb"] + d["wz"]), color=col, lw=1.5, ls=ls,
             solid_capstyle="round", dash_capstyle="round", label=lab)
    fin[rate] = (float(d["zb"][-1]), float(d["wz"][-1]),
                 float(d["zb"][-1] + d["wz"][-1]))

axA.set_xlim(2030, 170)
axA.set_ylim(-0.02, 1.02)
axA.set_xticks(np.arange(2000, 199, -300))
axA.set_yticks(np.arange(0, 1.01, 0.2))
axA.set_xlabel("Temperature (K)")
axA.set_ylabel("Crystallinity (ZB + WZ)")
axA.legend(loc="upper left", bbox_to_anchor=(0.015, 0.995), title="Cooling rate",
           title_fontsize=FS_LEG, labelspacing=0.32, borderaxespad=0.25,
           borderpad=0.2)
axA.text(0.975, 0.035, "(a)", transform=axA.transAxes, ha="right", va="bottom",
         fontsize=FS_PANEL)

rates = [225, 900]
xs = np.array(rates, float)
zbv = np.array([fin[r][0] for r in rates])
wzv = np.array([fin[r][1] for r in rates])
cv = np.array([fin[r][2] for r in rates])

axB.plot(xs, cv, "-s", color="#333333", ms=4.6, lw=1.1, mfc="none", mew=1.1,
         label="ZB + WZ", zorder=3)
axB.plot(xs, zbv, "-o", color="#C0392B", ms=4.6, lw=1.1, label="ZB", zorder=4)
axB.plot(xs, wzv, "--o", color="#2471A3", ms=4.6, lw=1.1, mfc="white", mew=1.1,
         label="WZ", zorder=4)
axB.set_xscale("log")
axB.xaxis.set_minor_formatter(mticker.NullFormatter())
axB.set_xticks(xs)
axB.set_xticklabels([str(r) for r in rates])
axB.set_xlim(150, 1400)
axB.set_ylim(-0.03, 1.06)
axB.set_yticks(np.arange(0, 1.01, 0.2))
axB.set_xlabel("Cooling rate (K/ns)")
axB.set_ylabel("Fraction at 200 K")
axB.legend(loc="center right", bbox_to_anchor=(1.0, 0.80), labelspacing=0.26,
           handlelength=1.75, handletextpad=0.6, borderaxespad=0.2, borderpad=0.15)
axB.text(0.04, 0.035, "(b)", transform=axB.transAxes, ha="left", va="bottom",
         fontsize=FS_PANEL)

for ext, kw in [("png", dict(dpi=600)), ("pdf", {}), ("eps", {})]:
    fig.savefig(os.path.join(SRC, f"cryst_4016_225_vs_900.{ext}"), **kw)

print("saved:", os.path.join(SRC, "cryst_4016_225_vs_900.png"))
for r in rates:
    z, w, c = fin[r]
    print(f"{r:>4} K/ns  ZB {z*100:5.1f}  WZ {w*100:5.1f}  other {(1-z-w)*100:5.1f}  "
          f"cryst {c:.3f}  ZB:WZ {z/max(w, 1e-9):.2f}")
