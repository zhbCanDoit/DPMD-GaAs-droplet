# -*- coding: utf-8 -*-
"""对比熔化轨迹末帧 与 弛豫轨迹首帧 是否同一构型(逐原子 id+type+最小镜像坐标)。"""
import numpy as np, os

MELT = r"E:/论文/GaAs纳米液滴结晶/Semi-GaAs-droplet/melt_nvt-300-2000k-3ns/GaAs-2000k.xyz"
RELAX = r"E:/论文/GaAs纳米液滴结晶/Semi-GaAs-droplet/relax/GaAs-2000k.xyz"
N = 4096
BOX = 100.0

def read_first_frame(path):
    atoms = {}
    with open(path) as f:
        in_atoms = False
        for ln in f:
            if "ITEM: ATOMS" in ln:
                in_atoms = True
                continue
            if in_atoms:
                p = ln.split()
                if len(p) >= 5:
                    atoms[int(p[0])] = (int(p[1]), float(p[2]), float(p[3]), float(p[4]))
                    if len(atoms) >= N:
                        break
    return atoms

def read_last_frame(path):
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        f.seek(max(0, size - 800 * 1024))
        data = f.read().decode("utf-8", "replace")
    lines = data.split("\n")
    last_atoms_idx = None
    for i, ln in enumerate(lines):
        if "ITEM: ATOMS" in ln:
            last_atoms_idx = i
    atoms = {}
    for ln in lines[last_atoms_idx + 1:]:
        p = ln.split()
        if len(p) >= 5:
            atoms[int(p[0])] = (int(p[1]), float(p[2]), float(p[3]), float(p[4]))
    return atoms

melt = read_last_frame(MELT)
relax = read_first_frame(RELAX)
print("melt 末帧原子数:", len(melt), "  relax 首帧原子数:", len(relax))

common = sorted(set(melt) & set(relax))
print("共同 id 数:", len(common))

type_mismatch = [i for i in common if melt[i][0] != relax[i][0]]
print("type 不一致的原子数:", len(type_mismatch))

m = np.array([melt[i][1:] for i in common])
r = np.array([relax[i][1:] for i in common])
d = r - m
d = d - np.round(d / BOX) * BOX   # 最小镜像
maxd = np.abs(d).max()
meand = np.abs(d).mean()
print(f"最小镜像坐标差: max = {maxd:.6f} A, mean = {meand:.6f} A")

if maxd < 0.01 and len(type_mismatch) == 0:
    print(">>> 结论: 熔化末帧 == 弛豫初始构型 (完全一致)")
else:
    print(">>> 结论: 不一致 (存在偏移或重编号)")
