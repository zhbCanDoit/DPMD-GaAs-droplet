# -*- coding: utf-8 -*-
"""4016 原子 GaAs 纳米液滴降温轨迹的结晶度扫描（两档速率：225 / 900 K/ns）

用法: python _cool4016_scan.py <tag>
      tag = run8ns (225 K/ns) | run2ns (900 K/ns)

输出: cool4016_<tag>.npz  (step, T, zb, wz, other) 与 csv 摘要
"""
import os
import sys
import time
import numpy as np
from ovito.io import import_file
from ovito.modifiers import IdentifyDiamondModifier

RUNS = {
    "run8ns": dict(
        path=r"E:\论文\GaAs纳米液滴结晶\Semi-GaAs-droplet\cooling\NVT-2000-200-8ns"
             r"\DzhbGaAs-DeepMDwork1Semi-GaAs-dropletcoolingNVT-2000-200-8nsGaAs-200k.dump",
        steps=8_000_000, rate=225.0, T0=2000.0, T1=200.0),
    "run2ns": dict(
        path=r"E:\论文\GaAs纳米液滴结晶\Semi-GaAs-droplet\cooling\nvt-2000-200k-2ns"
             r"\GaAs-200k.xyz",
        steps=2_000_000, rate=900.0, T0=2000.0, T1=200.0),
}

OUT = r"E:\论文\GaAs纳米液滴结晶\Semi-GaAs-droplet\cooling\rate_compare_4016"
os.makedirs(OUT, exist_ok=True)

tag = sys.argv[1]
cfg = RUNS[tag]
DUMPS = 50.0          # 每 50 步写一帧

t0 = time.time()
pipe = import_file(cfg["path"])
pipe.modifiers.append(IdentifyDiamondModifier())
nf = pipe.source.num_frames
print(f"[{tag}] frames={nf}  scan={time.time()-t0:.0f}s", flush=True)


def t_of(idx):
    """帧序号 -> 程序设定温度（线性降温）"""
    return cfg["T0"] - (cfg["T0"] - cfg["T1"]) * (idx / (nf - 1))


def measure(idx):
    data = pipe.compute(int(idx))
    st = np.asarray(data.particles["Structure Type"])
    cnt = np.bincount(st, minlength=7)
    n = st.size
    zb = (cnt[1] + cnt[2] + cnt[3]) / n
    wz = (cnt[4] + cnt[5] + cnt[6]) / n
    return zb, wz, 1.0 - zb - wz


# ---------- 第一轮：粗扫 120 帧定位转变区 ----------
coarse = np.unique(np.linspace(0, nf - 1, 120).astype(int))
res = {}
for i in coarse:
    res[int(i)] = measure(i)
tot = np.array([res[i][0] + res[i][1] for i in coarse])
peak = float(tot.max())
on = coarse[np.argmax(tot > 0.02)]
off = coarse[np.argmax(tot > 0.85 * peak)]
print(f"[{tag}] 粗扫: 起晶≈frame {on} (T≈{t_of(on):.0f} K), "
      f"到 85%峰值≈frame {off} (T≈{t_of(off):.0f} K), 峰值 {peak*100:.1f}%  "
      f"{time.time()-t0:.0f}s", flush=True)

# ---------- 第二轮：转变区加密 150 帧 ----------
span = max(off - on, 8)
lo = max(0, on - int(0.60 * span))
hi = min(nf - 1, off + int(0.60 * span))
dense = np.unique(np.linspace(lo, hi, 150).astype(int))
for i in dense:
    if int(i) not in res:
        res[int(i)] = measure(i)
print(f"[{tag}] 加密区 frame {lo}..{hi}  共 {len(res)} 帧  {time.time()-t0:.0f}s", flush=True)

idx = np.array(sorted(res))
zb = np.array([res[i][0] for i in idx])
wz = np.array([res[i][1] for i in idx])
ot = np.array([res[i][2] for i in idx])
Ts = np.array([t_of(i) for i in idx])
steps = idx * DUMPS

np.savez(os.path.join(OUT, f"cool4016_{tag}.npz"),
         step=steps, T=Ts, zb=zb, wz=wz, other=ot, rate=cfg["rate"])
print(f"[{tag}] 末帧: ZB {zb[-1]*100:.1f}%  WZ {wz[-1]*100:.1f}%  "
      f"cryst {(zb[-1]+wz[-1])*100:.1f}%  T_end={Ts[-1]:.0f} K")
print(f"[{tag}] total {time.time()-t0:.0f}s")
