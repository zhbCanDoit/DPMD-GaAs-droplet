# -*- coding: utf-8 -*-
"""
rdf_engine.py —— 从 LAMMPS dump 重算径向分布函数 g(r)

归一化方式（已用 RDF/200K 反推并完整复现，max|Δ| < 0.015）：
    g_AB(r) = V / (N_A * N_B) * dn_AB(r) / (4*pi*r^2*dr)
  - V   = 盒体积（dump 的 BOX BOUNDS 三边乘积）
  - dn  = 距离落在 [r-dr/2, r+dr/2) 内的有序原子对数，用最小镜像
  - A==B 时排除自配对
  - bin 中心取 dr/2, 3dr/2, ...（与原文件一致：0.05, 0.15, ...）

2026-10-02
"""
import os
import numpy as np

CHUNK = 512          # 分块大小，控制内存


# ---------------------------------------------------------------- 读单帧 dump
def _parse_frame_lines(lines):
    step = int(lines[1])
    lo = np.array([float(lines[5].split()[0]), float(lines[6].split()[0]), float(lines[7].split()[0])])
    hi = np.array([float(lines[5].split()[1]), float(lines[6].split()[1]), float(lines[7].split()[1])])
    i0 = next(k for k, l in enumerate(lines) if l.startswith("ITEM: ATOMS"))
    cols = lines[i0].split()[2:]
    a = np.array([[float(v) for v in l.split()] for l in lines[i0 + 1:] if l.split()])
    ix = {c: k for k, c in enumerate(cols)}
    typ = a[:, ix["type"]].astype(int)
    xyz = np.stack([a[:, ix["x"]], a[:, ix["y"]], a[:, ix["z"]]], axis=1)
    return step, lo, hi, typ, xyz


def read_dump_frame(path):
    with open(path, errors="ignore") as f:
        lines = f.read().splitlines()
    step, lo, hi, typ, xyz = _parse_frame_lines(lines)
    return step, hi - lo, typ, xyz


# ---------------------------------------------------------------- 从大 dump 取帧
def extract_frame(big_dump, target_step, total_steps=8000000, frame_step=50,
                  win=8_000_000, out=None):
    """按估计偏移 + 局部搜索取指定 step 的帧；out 给出则写出单帧 dump。"""
    size = os.path.getsize(big_dump)
    est = int(target_step / total_steps * size)
    with open(big_dump, "rb") as f:
        for w in (win, win * 4, win * 16):
            start = max(0, est - w // 2)
            f.seek(start)
            buf = f.read(w)
            pos = []
            p = buf.find(b"ITEM: TIMESTEP")
            while p != -1:
                pos.append(p)
                p = buf.find(b"ITEM: TIMESTEP", p + 1)
            if not pos:
                continue
            cand = []
            for i, p in enumerate(pos):
                nl = buf.find(b"\n", p) + 1
                nl2 = buf.find(b"\n", nl)
                try:
                    s = int(buf[nl:nl2].strip())
                except ValueError:
                    continue
                end = pos[i + 1] if i + 1 < len(pos) else len(buf)
                cand.append((s, p, end))
            if not cand:
                continue
            s, p, end = min(cand, key=lambda q: abs(q[0] - target_step))
            if abs(s - target_step) > frame_step * 0.8:
                continue
            lines = buf[p:end].decode("utf-8", "ignore").splitlines()
            st, lo, hi, typ, xyz = _parse_frame_lines(lines)
            if out:
                with open(out, "w") as g:
                    g.write("ITEM: TIMESTEP\n%d\n" % st)
                    g.write("ITEM: NUMBER OF ATOMS\n%d\n" % len(typ))
                    g.write("ITEM: BOX BOUNDS pp pp pp\n")
                    for k in range(3):
                        g.write("%.8g %.8g\n" % (lo[k], hi[k]))
                    g.write("ITEM: ATOMS type x y z id\n")
                    for k in range(len(typ)):
                        g.write("%d %.10g %.10g %.10g %d\n"
                                % (typ[k], xyz[k, 0], xyz[k, 1], xyz[k, 2], k + 1))
            return st, hi - lo, typ, xyz
    raise RuntimeError("frame step=%d not found" % target_step)


# ---------------------------------------------------------------- 算 g(r)
PAIRS = ((1, 1, "As-As"), (1, 2, "As-Ga"), (2, 2, "Ga-Ga"))


def compute_rdf(box, typ, xyz, rmax=10.0, dr=0.01, pairs=PAIRS, return_counts=False):
    """返回 {label: (r_center, g)}"""
    V = float(np.prod(box))
    nb = int(round(rmax / dr))
    r_center = (np.arange(nb) + 0.5) * dr
    phi = 4.0 * np.pi * r_center ** 2 * dr
    idx = {1: np.where(typ == 1)[0], 2: np.where(typ == 2)[0]}
    res = {}
    for a, b, lab in pairs:
        ia, ib = idx[a], idx[b]
        h = np.zeros(nb)
        for s in range(0, len(ia), CHUNK):
            sel = ia[s:s + CHUNK]
            d = xyz[sel][:, None, :] - xyz[ib][None, :, :]
            d -= box * np.round(d / box)                     # 最小镜像
            dist = np.sqrt((d ** 2).sum(-1))
            if a == b:
                for k in range(dist.shape[0]):
                    dist[k, s + k] = np.inf                  # 排除自身
            h += np.histogram(dist, bins=nb, range=(0.0, rmax))[0]
        g = h * (V / (len(ia) * len(ib))) / phi
        res[lab] = (r_center, g, h) if return_counts else (r_center, g)
    return res


def write_rdf(path, res, header_extra=""):
    r = next(iter(res.values()))[0]
    keys = [k for _, _, k in PAIRS if k in res]
    with open(path, "w") as f:
        f.write("# Radial distribution function (%d data points), dr = %.4f A\n" % (len(r), r[1] - r[0]))
        f.write("# %s\n" % header_extra)
        f.write('# "Pair separation distance" ' + " ".join(keys) + "\n")
        for i in range(len(r)):
            f.write("%.4f " % r[i] + " ".join("%.10g" % res[k][1][i] for k in keys) + " \n")
