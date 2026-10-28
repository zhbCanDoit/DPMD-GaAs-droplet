# -*- coding: utf-8 -*-
"""
绘制 GaAs 纳米液滴升温过程中的 As/Ga 优先蒸发曲线。
数据源: evaporation_curve.csv (extract_evap.py 生成, 35 A 半径阈值)
蒸发判据: 距质心 > 35 A 的原子计为"已蒸发"(飞离液滴)。
低温段(300~1000 K)的非零值来自初始立方体的角原子(表面失稳/球化前),
真正的蒸发从 ~1287 K 开始, 因此主图取 1100~2000 K 窗口。
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 读取 CSV
data = np.genfromtxt("evaporation_curve.csv", delimiter=",", skip_header=1,
                     names=["step", "T", "nGa", "nAs", "nGa_evap", "nAs_evap"])
T = data["T"]
ga = data["nGa_evap"].astype(float)
a_s = data["nAs_evap"].astype(float)

# 蒸发起始温度 (step 1741000 = 1287 K, 论文正文 "beyond 1741 ps")
onset_T = 1287.0

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1.15]})

# --- 左图: 全程 300~2000 K ---
ax = axes[0]
ax.plot(T, ga, color="#c0392b", lw=1.8, label="Ga")
ax.plot(T, a_s, color="#2471a3", lw=1.8, label="As")
ax.axvline(onset_T, color="gray", ls="--", lw=1.0)
ax.annotate("onset $\\approx$ 1287 K", xy=(onset_T, ax.get_ylim()[1] * 0.92),
            xytext=(1300, ax.get_ylim()[1] * 0.95), fontsize=8, color="gray")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Evaporated atoms")
ax.set_title("(a) full heating range")
ax.legend(frameon=False, fontsize=9)
ax.annotate("cube-corner artifact\n(pre-melting)", xy=(450, 25), fontsize=7.5,
            color="gray", ha="center")

# --- 右图: 1100~2000 K 放大(蒸发窗口) ---
ax = axes[1]
m = T >= 1100
ax.plot(T[m], ga[m], color="#c0392b", lw=2.0, label="Ga")
ax.plot(T[m], a_s[m], color="#2471a3", lw=2.0, label="As")
ax.axvline(onset_T, color="gray", ls="--", lw=1.0)
ax.annotate("onset $\\approx$ 1287 K", xy=(onset_T, 2), xytext=(1350, 12),
            fontsize=8.5, color="gray",
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))
# 末点标注
ax.annotate(f"As: {a_s[-1]:.0f}", xy=(T[-1], a_s[-1]), xytext=(-30, 6),
            textcoords="offset points", fontsize=8.5, color="#2471a3")
ax.annotate(f"Ga: {ga[-1]:.0f}", xy=(T[-1], ga[-1]), xytext=(-30, -12),
            textcoords="offset points", fontsize=8.5, color="#c0392b")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Evaporated atoms")
ax.set_title("(b) evaporation window")
ax.legend(frameon=False, fontsize=9, loc="upper left")

fig.suptitle("Preferential As evaporation during GaAs nanodroplet melting (NVT 300$\\rightarrow$2000 K)",
             fontsize=10.5, y=1.02)
fig.tight_layout()

for ext in ("png", "pdf", "svg"):
    fig.savefig(f"evaporation_curve.{ext}", bbox_inches="tight", dpi=300)

print("saved evaporation_curve.png / .pdf / .svg")
print(f"final @ {T[-1]:.0f} K: Ga={ga[-1]:.0f}, As={a_s[-1]:.0f}, As/Ga={a_s[-1]/max(ga[-1],1):.2f}")
