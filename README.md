# Data and code for "Revealing Phase Transition Mechanism of GaAs Nanoparticle via Machine Learning Molecular Dynamics Simulation"

Hongbin Zhang, Mengru Huang, Haichao Li, Zijun Meng, Senze Gu, Lina Xu*, Hongping Xiao, Guoyong Fang*

College of Chemistry and Materials Engineering, Wenzhou University, Wenzhou 325035, China

Correspondence: xulina@wzu.edu.cn, fanggy@wzu.edu.cn

---

## 1. What this archive contains

Every file in this archive is a verbatim copy of the file used to produce the results reported in the
paper and in the Supplementary Information. No file has been edited, renumbered, reformatted or
regenerated for the archive; `_filelist.tsv` gives the path, size and MD5 sum of each one, and the MD5
sums are those of the original files. Files have only been renamed when a directory contained several
files of the same name (see Section 3).

Included: the interatomic potential, the input structures, the LAMMPS input files, the post-processing
scripts, and the numerical data underlying the figures.

Not included: the full MD trajectories (Section 6) and all figure files.

Total size: 40.1 MB in 112 files, plus this README and `_filelist.tsv`.

Only the three systems actually simulated in the paper and the Supplementary Information are
represented here: the 4096-atom crystal used for the melting study (which yields the 4016-atom
droplet), the 4016-atom droplet cooled at 225 K/ns, and the 976-atom droplet used for the
cooling-rate scan of Section S2. Configurations and input files from exploratory runs that do not
underlie any reported result are not part of the archive.

## 2. Atom-type convention

The type indices are **not** the same in every file, and this is a property of the original working
files rather than of the archive. The melting stage was run with type 1 = Ga; the configurations
written out after the melting run were re-exported through OVITO with the elements in the opposite
order, so everything downstream of that export — the 2000 K relaxation, the crystallization run, the
976-atom scan and the frame files — uses type 1 = As.

| File | type 1 | type 2 | How this was established |
|---|---|---|---|
| `GaAs4096atoms_initial_ZB.data` | Ga | As | `in.melt_300-2000K.lammps` maps `pair_coeff * * Ga As` |
| `melt_lastFrame_2000K.data` | As | Ga | `in.relax_2000K.lammps` maps `pair_coeff * * As Ga` with mass 1 = 74.92 |
| `droplet_4016_initial_2000K.data` | As | Ga | 1990 As + 2026 Ga; `in.cooling_...lammps` maps `As Ga` |
| `droplet_4016_final_200K.data` | As | Ga | 1990 As + 2026 Ga; coordinates matched atom by atom against `frames_CN/200K.xyz` |
| `droplet_976atoms_initial_1800K.data` | As | Ga | stated in the header of `in.scan976_100Kns.lammps` |
| `frames_CN/*.xyz` | — | — | element symbols are given explicitly, no type index is used |
| `frames_RDF/*.dump`, `frames_RDF_recomputed/*.dump` | As | Ga | coordinates matched atom by atom against `frames_CN/200K.xyz` |
| `04_data/energy_projection/dumpdropletLast.pe` | As | Ga | written by the crystallization run |

The correct mapping is always given by the `pair_coeff` and `mass` commands of the input file that
produced a given file, and those commands are left exactly as they were run. **Using a structure with
the wrong mapping silently exchanges Ga and As**, which changes the masses, invalidates the evaporation
counts and the coordination-number statistics, and is not detectable from the energies alone.

## 3. How to run: file names used inside the inputs

The LAMMPS inputs refer to local file names, as they did on the machine where the runs were performed.
The correspondence with this archive is

| Input file | name used in `read_data` | file in this archive |
|---|---|---|
| `in.melt_300-2000K.lammps` | `./GaAs4096atoms` | `01_structures/GaAs4096atoms_initial_ZB.data` |
| `in.relax_2000K.lammps` | `./lastFrame` | `01_structures/melt_lastFrame_2000K.data` |
| `in.cooling_2000-200K_225Kns.lammps` | `./droplet2` | `01_structures/droplet_4016_initial_2000K.data` |
| `in.scan976_{50,100,200,400}Kns.lammps` | `./GaAs-last-noC` | `01_structures/droplet_976atoms_initial_1800K.data` |

The potential is referenced under two names, `./GaAs-semi.pb` in all inputs except
`in.scan976_400Kns.lammps`, which uses
`../GaAs_0-3050K_0-50000bar_7045frame_train800w_AddVirial.pth`. Both denote the single file

```
05_model/GaAs_0-3050K_0-50000bar_7045frame_train800w_AddVirial.pb
```

which is the GaAs model of the universal IIB-VIA semiconductor potential set of Liu et al. (ref. 24),
taken unmodified from the `5-GaAs` directory of that release (MD5
`38a694b599c859dc3bde20b841a72768`, 6,833,628 bytes). Each of the four copies present in the original
working directories has this same MD5. No training or fine-tuning was performed on it for this work,
and the file is reproduced here only so that the runs can be repeated. The `.pth` extension in that
one input is a historical naming artefact of the cluster on which the run was made; the model itself
is the same.

To run a case, place the structures and the potential where the input expects them, or edit only the
`read_data` and `pair_style` lines to point at `01_structures/` and `05_model/`.

## 4. Contents

```
01_structures/      input and selected configurations (LAMMPS data / xyz / dump)
02_inputs/          LAMMPS input files
03_scripts/         post-processing scripts (Python 3)
04_data/            numerical data underlying the figures
05_model/           interatomic potential
_filelist.tsv       path, size and MD5 of every file
```

### 01_structures

| File | Description |
|---|---|
| `GaAs4096atoms_initial_ZB.data` | Initial zincblende GaAs crystal, 4096 atoms, in a 100 x 100 x 100 A box. Starting point of the melting run. |
| `melt_lastFrame_2000K.data` | Last configuration of the melting run (2000 K). Input of the 200 ps relaxation. |
| `relax2000K_final.xyz` | Last frame of the 2000 K relaxation, before the overflowing atoms were removed. |
| `droplet_4016_initial_2000K.data` | Ga-rich droplet after removal of the overflowing atoms, 4016 atoms (2026 Ga + 1990 As). Starting point of the 225 K/ns crystallization run. |
| `droplet_4016_final_200K.data` | Final crystallized configuration at 200 K; the structure shown in Fig. 3(d) and Fig. 4(c-e). |
| `droplet_976atoms_initial_1800K.data` | 976-atom Ga-rich droplet used for the cooling-rate scan of Section S2. |
| `frames_CN/*.xyz` | Frames of the 225 K/ns run at 200, 400, 600, 700, 800, 900, 1000 and 1800 K, used for the coordination-number analysis. |
| `frames_RDF/200K.dump` | 200 K frame used for the pair distribution functions of Fig. 3. |
| `frames_RDF_recomputed/*.dump` | Frames at 200, 400, 700, 800, 1000 and 1800 K used for the recomputed, finer-grid PDFs. |

### 02_inputs

| File | Description |
|---|---|
| `in.melt_300-2000K.lammps` | Heating 300 to 2000 K over 3 ns, NVT, Nose-Hoover, 0.1 ps damping. |
| `in.relax_2000K.lammps` | Relaxation at 2000 K for 200 ps. |
| `in.cooling_2000-200K_225Kns.lammps` | Cooling 2000 to 200 K at 225 K/ns over 8 ns; production run of the main text. |
| `in.scan976_{50,100,200,400}Kns.lammps` | Cooling of the 976-atom droplet at 50, 100, 200 and 400 K/ns (Section S2). Same initial configuration, same potential, same seed; only the number of steps differs. |

All runs use `units metal`, `timestep 0.001` ps and `velocity all create <T> 23456789`. The melting,
relaxation and crystallization inputs use `fix nvt` with a 0.1 ps damping constant; the four 976-atom
inputs use a damping constant of 1 ps.

### 03_scripts

| File | Purpose |
|---|---|
| `extract_evap.py` | Counts detached atoms during melting (detached when further than 35 A from the instantaneous centre of mass); produces `04_data/melting/evaporation_curve.csv` and Fig. S3. |
| `plot_evap.py` | Plots the evaporation curves. |
| `compare_frames.py` | Compares the last melting frame with the first relaxation frame atom by atom (Fig. 2 snapshots). |
| `extract_evap_relax.py` | Same evaporation analysis for the relaxation stage. |
| `analyze_cn.py` | Coordination-number statistics with a 4.0 A hetero-coordination cut-off (Fig. 4). |
| `rdf_engine.py` | Pair distribution function engine with a selectable grid spacing. |
| `step2_recompute.py`, `step3_figures.py` | Recompute the PDFs on a fine grid and build the corresponding figures. |
| `step6b_shoulder_all_temps.py` | Resolves the two maxima of the first PDF peak at every temperature (Fig. 3 insets). |
| `_cool4016_scan.py`, `_plot_cool4016.py`, `plot_FigS4_4016.py` | Crystallinity and ZB/WZ analysis of the 4016-atom droplet (Fig. S4). `_cool4016_scan.py` retains the two rate branches used during development; the archive contains the trajectory statistics of the 225 K/ns branch only (`cool4016_run8ns.npz`), which is the one reported in Section S4. |

The scripts contain the absolute paths of the machine on which they were written; they are archived as
run, not cleaned up. ZB and WZ labelling was done in OVITO with *Identify Diamond Structure*
(CUBIC_DIAMOND = ZB, HEX_DIAMOND = WZ) following the scheme of ref. 29; crystallinity is the sum of the
ZB and WZ fractions.

### 04_data

| Path | Content | Figure |
|---|---|---|
| `melting/log.lammps` | Thermo output of the heating run (step, pe, ke, etotal, temp, press, vol). | Fig. 1(a) |
| `melting/{GaAs,Ga,As}_msd_data.txt` | Mean square displacement against time, total and per species. | Fig. 1(b) |
| `melting/evaporation_curve.csv` | step, temp_K, nGa_total, nAs_total, nGa_evap, nAs_evap. | Fig. S3 |
| `relaxation/log.lammps`, `relaxation/relaxation_evap.csv` | Thermo output and evaporation analysis of the 2000 K relaxation. | — |
| `crystallization/log.lammps` | Thermo output of the 225 K/ns run. | Fig. 1(a) |
| `crystallization/GaAs2000-200K-Processed_data.data` | Crystallinity against temperature (OVITO structure counts). | Fig. S4 |
| `crystallization/Ga_msd_data.txt` | Ga mean square displacement during cooling. | — |
| `energy_projection/dumpdropletLast.pe` | Per-atom potential energy of the final configuration. | Fig. 5 |
| `RDF_225Kns/<T>` | PDFs of the 225 K/ns run on the coarse grid, 18 temperatures. Columns: r, As-As, As-Ga, Ga-Ga. | Fig. 3 |
| `RDF_recomputed/<T>/` | PDFs on finer grids (0.01, 0.02, 0.1 A), first-peak fits and the two-maxima decomposition. | Fig. 3 insets |
| `crystallization_4016_rate/cool4016_run8ns.npz` | Per-frame crystallinity and ZB/WZ fractions of the 4016-atom droplet at 225 K/ns. | Fig. S4 |
| `scan_976atom/` | Ga and As mean square displacement of the 976-atom runs. | Fig. S2 |
| `crystal_reference/` | Neighbour-shell statistics of ideal ZB and WZ GaAs used as the labelling reference. | Fig. S4 |

## 5. Reproduction sequence

1. Heat `GaAs4096atoms_initial_ZB.data` with `in.melt_300-2000K.lammps` (300 to 2000 K over 3 ns).
2. Relax the last frame at 2000 K for 200 ps with `in.relax_2000K.lammps`.
3. Remove the atoms lying outside the spherical cut-off to obtain the 4016-atom Ga-rich droplet
   (`droplet_4016_initial_2000K.data`).
4. Cool it with `in.cooling_2000-200K_225Kns.lammps` (2000 to 200 K at 225 K/ns).
5. Post-process with the scripts in `03_scripts`.

The cooling-rate scan of Section S2 runs independently from `droplet_976atoms_initial_1800K.data`
with `in.scan976_{50,100,200,400}Kns.lammps`, i.e. the same initial configuration and seed cooled
from 1800 to 200 K in 32, 16, 8 and 4 ns respectively.

## 6. What is not included

The complete trajectories are too large to distribute with the paper and are available from the authors
on reasonable request.

| Run | File written | Size |
|---|---|---|
| Melting, 300 to 2000 K | `GaAs-2000k.xyz`, every 50 steps | 7.5 GB |
| Relaxation, 2000 K, 200 ps | `GaAs-2000k.xyz` | 250 MB |
| Crystallization, 225 K/ns | `GaAs-200k.dump` | 19.5 GB |
| Rate scan, 976-atom droplet | `GaAs.xyz`, one per rate | 23 to 93 MB each |

The selected frames in `01_structures/` are sufficient to reproduce every figure.

## 7. Notes

- The 4016-atom droplet was obtained by deleting the atoms lying outside a spherical cut-off in the
  last frame of the 2000 K relaxation. The initial crystal is stoichiometric (2048 Ga + 2048 As); the
  last relaxation frame still contains 4096 atoms, of which 80 lie outside the cut-off and are
  removed, leaving 2026 Ga + 1990 As. The count of atoms recorded as detached (further than 35 A from
  the instantaneous centre of mass) is a separate diagnostic and is tabulated in
  `04_data/melting/evaporation_curve.csv` and `04_data/relaxation/relaxation_evap.csv`.
- The energy and force comparison of Fig. S1 uses the training and test sets of ref. 24 and is therefore
  not redistributed here; it can be reproduced from the data set of that work together with the
  potential in `05_model/`.
- The absolute melting temperature obtained with this potential is tied to its PBE reference and would
  shift with a different exchange-correlation functional; the sequence of mechanistic events reported in
  the paper would not.
