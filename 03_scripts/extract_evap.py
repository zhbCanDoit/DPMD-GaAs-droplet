#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract_evap.py —— 从熔化轨迹 GaAs-2000k.xyz（LAMMPS dump 格式）统计蒸发曲线。

type 约定（in.lammps）: type 1 = Ga, type 2 = As
升温: 300 -> 2000 K, 0.567 K/ps, timestep 0.001 ps, dump 每 50 步一帧

判据: 原子距质心 > threshold (35 A) 记为"蒸发/游离"。
自洽验证: 论文最终移除 22 Ga + 58 As，脚本末帧应复现该值。
"""
import numpy as np

PATH = r'E:/论文/GaAs纳米液滴结晶/Semi-GaAs-droplet/melt_nvt-300-2000k-3ns/GaAs-2000k.xyz'
OUT  = r'E:/论文/GaAs纳米液滴结晶/Semi-GaAs-droplet/melt_nvt-300-2000k-3ns/evaporation_curve.csv'

THRESHOLD = 35.0          # Angstrom，液滴半径 ~28 A，留 7 A 缓冲
T0, RAMP, DT = 300.0, 0.567, 0.001

DIAG_STEPS = {0, 300000, 600000, 1000000, 1500000, 1741000, 2000000, 2500000, 3000000}

rows = []
with open(PATH, 'r', errors='ignore') as f:
    while True:
        line = f.readline()
        if not line:
            break
        if not line.startswith('ITEM: TIMESTEP'):
            continue
        step = int(f.readline())
        f.readline()                       # ITEM: NUMBER OF ATOMS
        natom = int(f.readline())
        f.readline()                       # ITEM: BOX BOUNDS
        for _ in range(3):
            f.readline()
        f.readline()                       # ITEM: ATOMS id type x y z

        typ = np.empty(natom, dtype=np.int32)
        xyz = np.empty((natom, 3), dtype=np.float64)
        for i in range(natom):
            p = f.readline().split()
            typ[i] = int(p[1])
            xyz[i, 0] = float(p[2]); xyz[i, 1] = float(p[3]); xyz[i, 2] = float(p[4])

        com = xyz.mean(axis=0)
        d = np.sqrt(((xyz - com) ** 2).sum(axis=1))
        ev = d > THRESHOLD

        nGa = int((typ == 1).sum()); nAs = int((typ == 2).sum())
        nGe = int(((typ == 1) & ev).sum()); nAe = int(((typ == 2) & ev).sum())
        temp = T0 + step * DT * RAMP
        rows.append((step, round(temp, 1), nGa, nAs, nGe, nAe))

        if step in DIAG_STEPS:
            q = np.percentile(d, [0, 5, 25, 50, 75, 90, 99, 100])
            print(f'step={step:8d} T={temp:6.0f}K  Ga={nGa} As={nAs}  '
                  f'游离 Ga={nGe} As={nAe}  (Ga+As={nGe+nAe})', flush=True)
            print('    dist[A] p0/5/25/50/75/90/99/100 = '
                  + ' '.join(f'{v:.1f}' for v in q), flush=True)
        if step % 250000 == 0:
            print(f'[progress] step={step} T={temp:.0f}K frames={len(rows)}', flush=True)

with open(OUT, 'w') as g:
    g.write('step,temp_K,nGa_total,nAs_total,nGa_evap,nAs_evap\n')
    for r in rows:
        g.write(','.join(map(str, r)) + '\n')

# 末帧汇总
last = rows[-1]
print('\n=== 末帧 step=%d T=%.0fK: 游离 Ga=%d As=%d (合计 %d) ===' %
      (last[0], last[1], last[4], last[5], last[4] + last[5]))
print('论文期望: 移除 22 Ga + 58 As = 80 原子')
print('DONE  frames=%d  csv=%s' % (len(rows), OUT))
