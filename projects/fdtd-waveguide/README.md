# FDTD Simulation of a 1-D Cavity and a 2-D Rectangular Waveguide

*Virginia Tech ECE 5106 (EM Fields II) course project, spring 2016. Written in MATLAB.*

A from-scratch finite-difference time-domain (FDTD) solver using the Yee scheme. It finds the resonant modes of a 1-D cavity and the modes of a PEC rectangular waveguide, and the report checks both against the analytical solutions.

![2-D FDTD Hy field snapshot](../../assets/fdtd-2d-waveguide.png)

*H<sub>y</sub> field snapshot from the 2-D simulation (figure from the project report).*

## What it does

**1-D cavity:** `FDTD_1D_Take2_C.m`
- PEC-bounded cavity, L = 30, N = 31 nodes, normalized units (ε = μ = 1), Δt = 0.5 Δz, 10,000 time steps.
- Takes the FFT of E<sub>x</sub> and H<sub>y</sub> to find the resonant frequencies, then compares them with the analytical modes f<sub>m</sub> = m / 2L.
- As committed, it starts from a single m = 5 mode. The report also uses multi-mode and random starting fields and a spatial FFT filter to isolate individual modes; those options are present but commented out in the script.

**2-D waveguide cross-section:** `ECE_5106_FDTD_Project_2D_A6.m`
- 0.1 m × 0.05 m PEC rectangular waveguide on a 51 × 51 grid, with Δt set to 0.4 × the Courant stability limit and 10,000 time steps.
- The update equations include the propagation constant β. As committed, the script runs one case: β = 0, starting from a TE<sub>10</sub> field distribution (H<sub>z</sub> ∝ cos(πx/a)).
- Records the fields over time and takes FFTs to find the mode frequencies, then animates the H<sub>z</sub> field.

**In the report only:** sweeping β to compare the numerical β–f dispersion curve with the analytical one (close agreement), and varying Δt, Δx and the number of time steps to show how each trades accuracy against run time. The committed 2-D script has these parameters fixed.

## Files

| File | Purpose |
|------|---------|
| `FDTD_1D_Take2_C.m` | 1-D cavity FDTD simulation and modal FFT analysis |
| `ECE_5106_FDTD_Project_2D_A6.m` | 2-D waveguide FDTD simulation and dispersion analysis |

## Running

Both scripts are standalone. Open MATLAB in this folder and run either one by name, e.g. `FDTD_1D_Take2_C`. Each script starts with `clear all; close all`, so save your workspace first. Neither has been tested in GNU Octave.

**The 2-D script is heavy.** It stores every time step of all six field components (about 1.25 GB), prints the step number on each of the 10,000 iterations, and ends by animating all 10,000 frames. To skip the animation, set `animate_fields = 0` at line 444 before running.

## Known issues and a modern version

A 2026 code review found these problems in the original scripts. The code here is left as it was written in 2016.

- `FDTD_1D_Take2_C.m`, line 39: the "single-mode" starting field uses `(1:N)` where `(0:N-1)` is needed, so it is shifted by one node and is nonzero on the PEC walls. Zeroing the walls then excites other cavity modes.
- `ECE_5106_FDTD_Project_2D_A6.m`: the FFT probe sits at the guide's center, which is a node of H<sub>z</sub> for the TE<sub>10</sub> mode, so the dominant mode is weak in the computed spectrum.
- `ECE_5106_FDTD_Project_2D_A6.m`: the Gaussian pulse is computed but never used as the source, and β is fixed at 0.

**Modern version:** [FDTD from scratch: a 1-D cavity and a 2-D waveguide](https://github.com/jman4162/RF-and-Antenna-Tutorials/blob/main/FDTD_From_Scratch_Cavity_and_Waveguide_Tutorial.ipynb) (Python notebook). It fixes these issues, sweeps β to build the full dispersion diagram, and checks every resonance against the scheme's exact discrete dispersion relation.

## Report

[Full project report (PDF, 55 pages)](../../Hodge_John_ECE_5106_FDTD_Project_Doc_Draft_C.pdf): derivations, update equations, boundary conditions, results and the parameter studies.
