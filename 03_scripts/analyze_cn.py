#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
analyze_cn.py —— 冷却速率扫描的配位数(CN)统计

读 LAMMPS custom dump（多帧），取最后 N 帧，计算每个原子的异类配位数
（最小镜像 + 截断 4.0 Å：As 只数 Ga 近邻、Ga 只数 As 近邻），
统计 CN=3(表面)/4(体相)/5(缺陷) 的原子数（分 As / Ga）。

对应正文 Figure 5 的口径（2026-10-03 实测破解）：
  - 异类配位（不是总配位），截断 4.0 Å（= ZB 次近邻距离 a/√2，a=5.65 Å）
  - 该口径对 200K_step8000000.dump 实测复现正文：Ga CN3=218/CN5=244，As CN3=201/CN5=311
  - OVITO 等价设置：Coordination analysis, cutoff=4.0 Å，Ga-Ga、As-As 同类对截断设 0
type 1 = As，type 2 = Ga（与 in.cool_*.lammps 的 group 定义一致）。

用法:
  python analyze_cn.py --dump rate50.dump rate100.dump rate225.dump rate500.dump --last 10
"""
import argparse
import csv
import os

import numpy as np

CUT = 4.0   # 异类配位截断（= ZB 次近邻 a/√2≈3.996 Å；只覆盖第一壳 Ga-As 键 2.45 Å）


def read_frames(path):
    """读多帧 custom dump，按时间顺序返回 [(step, box, typ, xyz)]。"""
    with open(path, errors="ignore") as f:
        lines = f.read().splitlines()
    frames, i, n = [], 0, len(lines)
    while i < n:
        if lines[i].strip() == "ITEM: TIMESTEP":
            step = int(lines[i + 1]); i += 2
            natom = int(lines[i + 1]); i += 2
            i += 1                              # 跳过 "ITEM: BOX BOUNDS ..." 行
            box = []
            for _ in range(3):
                p = lines[i].split()
                box.append(float(p[1]) - float(p[0])); i += 1
            assert lines[i].startswith("ITEM: ATOMS"), lines[i]
            cols = lines[i].split()[2:]; i += 1
            a = np.array([[float(v) for v in l.split()]
                          for l in lines[i:i + natom] if l.split()])
            ix = {c: k for k, c in enumerate(cols)}
            typ = a[:, ix["type"]].astype(int)
            xyz = np.stack([a[:, ix["x"]], a[:, ix["y"]], a[:, ix["z"]]], axis=1)
            frames.append((step, np.array(box, float), typ, xyz))
            i += natom
        else:
            i += 1
    return frames


def cn_counts(box, typ, xyz, cut=CUT):
    """异类配位：As 只数 Ga 近邻、Ga 只数 As 近邻（正文 Figure 5 口径）。
    最小镜像，分块避免 O(N^2) 内存。"""
    natom = len(typ)
    boxa = box.astype(float)
    cn = np.zeros(natom, int)
    iAs = np.where(typ == 1)[0]   # As 原子
    iGa = np.where(typ == 2)[0]   # Ga 原子
    CH = 512
    # As 原子 -> 数 Ga 近邻
    for s in range(0, len(iAs), CH):
        e = min(s + CH, len(iAs))
        d = xyz[iAs[s:e]][:, None, :] - xyz[iGa][None, :, :]
        d -= boxa * np.round(d / boxa)
        dist = np.sqrt((d ** 2).sum(-1))
        cn[iAs[s:e]] = (dist < cut).sum(axis=1)
    # Ga 原子 -> 数 As 近邻
    for s in range(0, len(iGa), CH):
        e = min(s + CH, len(iGa))
        d = xyz[iGa[s:e]][:, None, :] - xyz[iAs][None, :, :]
        d -= boxa * np.round(d / boxa)
        dist = np.sqrt((d ** 2).sum(-1))
        cn[iGa[s:e]] = (dist < cut).sum(axis=1)
    return cn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", nargs="+", required=True)
    ap.add_argument("--last", type=int, default=10, help="取最后 N 帧做平均")
    ap.add_argument("--out", default="cn_summary.csv")
    args = ap.parse_args()

    header = ("dump", "frames", "As_CN3", "As_CN4", "As_CN5",
              "Ga_CN3", "Ga_CN4", "Ga_CN5")
    print("%-14s %-6s %-7s %-7s %-7s %-7s %-7s %-7s" % header)
    rows = []
    for p in args.dump:
        frames = read_frames(p)
        last = frames[-args.last:]
        acc = {1: np.zeros(9, int), 2: np.zeros(9, int)}
        for step, box, typ, xyz in last:
            cn = cn_counts(box, typ, xyz)
            for t in (1, 2):
                m = typ == t
                c = np.bincount(cn[m], minlength=9)
                acc[t][:len(c)] += c
        avg = {t: np.round(acc[t] / len(last)).astype(int) for t in (1, 2)}
        a3, a4, a5 = avg[1][3], avg[1][4], avg[1][5]
        g3, g4, g5 = avg[2][3], avg[2][4], avg[2][5]
        print("%-14s %-6d %-7d %-7d %-7d %-7d %-7d %-7d"
              % (os.path.basename(p), len(last), a3, a4, a5, g3, g4, g5))
        rows.append([os.path.basename(p), len(last), a3, a4, a5, g3, g4, g5])

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(header))
        w.writerows(rows)
    print("-> 写出 %s" % args.out)


if __name__ == "__main__":
    main()
