# -*- coding: utf-8 -*-
"""提取弛豫段(2000 K 恒温 2 ns)的 As/Ga 蒸发曲线, 判据与熔化段一致: 距质心 > 35 A。"""
import numpy as np

PATH = r"E:/论文/GaAs纳米液滴结晶/Semi-GaAs-droplet/relax/GaAs-2000k.xyz"
OUT = r"E:/论文/GaAs纳米液滴结晶/Semi-GaAs-droplet/relax/relaxation_evap.csv"
R = 35.0
BOX = 100.0

def process(coords, types, step):
    c = np.array(coords)
    com = c.mean(0)
    d = c - com
    d = d - np.round(d / BOX) * BOX
    r = np.sqrt((d ** 2).sum(1))
    t = np.array(types)
    return int((t[r > R] == 1).sum()), int((t[r > R] == 2).sum())  # Ga, As

with open(OUT, "w") as fo:
    fo.write("step,time_ps,nGa_evap,nAs_evap\n")
    step = None
    coords, types = [], []
    with open(PATH) as f:
        for ln in f:
            s = ln.strip()
            if s.startswith("ITEM: TIMESTEP"):
                if step is not None and coords:
                    g, a = process(coords, types, step)
                    fo.write(f"{step},{step*0.001:.2f},{g},{a}\n")
                coords, types = [], []
                step = None
            elif s.startswith("ITEM:"):
                pass
            elif step is None:
                try:
                    step = int(s)
                except ValueError:
                    pass
            else:
                p = s.split()
                if len(p) >= 5:
                    types.append(int(p[1]))
                    coords.append([float(p[2]), float(p[3]), float(p[4])])
        if step is not None and coords:
            g, a = process(coords, types, step)
            fo.write(f"{step},{step*0.001:.2f},{g},{a}\n")

print("done ->", OUT)
