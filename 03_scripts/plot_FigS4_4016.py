# -*- coding: utf-8 -*-
"""SI Figure S4：4016 原子富镓 GaAs 纳米液滴在 225 K/ns 冷却过程中
ZB / WZ / 总结晶度 随【温度】的变化（单面板，无插图）。

版式与 SI Figure S2（cryst_vs_rate_alt100new_v2）完全对齐：
  - 画布 14.65 x 6.80 cm（= SI 正文宽度，插入 1:1 不缩放）
  - 轴标题 10 pt / 刻度 9.5 pt / 图例 9 pt（与 S2 同一组常量）
  - 线宽、刻度长度、四边朝内刻度等 rcParams 逐项相同
  - 温度轴降序（冷却方向从左到右）
  - 配色：ZB 红 #C0392B、WZ 蓝 #2471A3、ZB+WZ 灰 #666666
数据：rate_compare_4016/cool4016_run8ns.npz（源自 NVT-2000-200-8ns 的 OVITO 逐帧统计）
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = r"E:/论文/GaAs纳米液滴结晶/Semi-GaAs-droplet/cooling/rate_compare_4016"
OUT = SRC
NPZ = os.path.join(SRC, "cool4016_run8ns.npz")

# ---- 字号：与 _plot_sci_v3.py（Figure S2）完全相同 ----
FS_LAB, FS_TICK, FS_LEG = 10.0, 9.5, 9.0

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

d = np.load(NPZ)
T = d["T"]
# 不做平滑：末帧数值须与正文/SI 引用的 39.9 / 48.7 / 88.6 严格一致
zb = np.asarray(d["zb"], float)
wz = np.asarray(d["wz"], float)
tot = np.asarray(d["zb"] + d["wz"], float)

W_CM, H_CM = 14.65, 6.80
fig = plt.figure(figsize=(W_CM / 2.54, H_CM / 2.54))
ax = fig.add_axes([0.090, 0.185, 0.881, 0.770])

ax.plot(T, tot * 100, color="#666666", lw=1.5, zorder=2, label="ZB + WZ")
ax.plot(T, zb * 100, color="#C0392B", lw=1.5, zorder=3, label="ZB")
ax.plot(T, wz * 100, color="#2471A3", lw=1.5, zorder=3, label="WZ")

ax.set_xlim(2000, 200)
ax.set_ylim(-1.5, 108)
ax.set_xticks(np.arange(2000, 199, -200))
ax.set_xticks(np.arange(1900, 199, -100), minor=True)
ax.set_yticks(np.arange(0, 101, 20))
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Crystallinity (%)")
ax.legend(loc="upper left", bbox_to_anchor=(0.015, 0.995), labelspacing=0.32,
          borderaxespad=0.25, borderpad=0.2, ncol=1)

for ext, kw in [("png", dict(dpi=600)), ("pdf", {}), ("eps", {})]:
    fig.savefig(os.path.join(OUT, f"FigS4_4016_ZB_WZ_vs_T.{ext}"), **kw)

print("canvas %.2f x %.2f cm  |  fonts %.1f / %.1f / %.1f pt" % (W_CM, H_CM, FS_LAB, FS_TICK, FS_LEG))
print("末帧 T=%.0f K  ZB %.1f%%  WZ %.1f%%  总 %.1f%%"
      % (T[-1], zb[-1] * 100, wz[-1] * 100, tot[-1] * 100))
print("起晶(>2%%) T=%.0f K；>50%% T=%.0f K"
      % (T[np.argmax(tot > 0.02)], T[np.argmax(tot > 0.5)]))
print("saved:", os.path.join(OUT, "FigS4_4016_ZB_WZ_vs_T.png"))
