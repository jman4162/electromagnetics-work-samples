# EM Boundaries: Multilayer Dielectric Reflection and a Slab Waveguide

*Virginia Tech ECE 5105 final exam, December 2015. Written in MATLAB.*

Transfer-matrix analysis of plane waves crossing a three-layer dielectric stack, plus the guided TM mode of a dielectric slab waveguide.

![Reflected and transmitted power vs. middle-layer thickness](../../assets/dielectric-slab-reflection.png)

*Reflected and transmitted power vs. middle-layer thickness. In Case A the middle layer is optically thinner, so the wave tunnels through and transmission decays as the layer gets thicker (frustrated total internal reflection). In Case B it is denser, giving periodic Fabry–Pérot-style resonances. Figure from the exam report.*

## What it does

**Problem 1: TE plane wave through a three-layer stack** (`Hodge_John_ECE_5105_Final_Exam_P1.m`)
- TE-polarized wave incident at 60° with ω = 3 × 10<sup>11</sup> rad/s (λ ≈ 6.3 mm); outer layers ε<sub>1</sub> = ε<sub>3</sub> = 4ε<sub>0</sub>.
- Sweeps the middle-layer thickness from 0.01 to 10 mm for two cases: ε<sub>2</sub> = ε<sub>0</sub> (Case A) and ε<sub>2</sub> = 8ε<sub>0</sub> (Case B).
- Computes the reflected and transmitted power fractions from the boundary-condition transfer matrices.

**Problem 2: TM dielectric slab waveguide** (`..._P2a.m`, `..._P2b.m`)
- 0.5 µm-thick core at λ = 1 µm.
- **2a:** finds the propagation constant β of the fundamental TM mode graphically, from the zero of the transfer-matrix element (T<sub>TM</sub>)<sub>11</sub>.
- **2b:** plots the H<sub>y</sub> field profile and the z-component of the time-averaged Poynting vector across both cladding layers and the core.

## Files

| File | Purpose |
|------|---------|
| `Hodge_John_ECE_5105_Final_Exam_P1.m` | Reflection and transmission vs. layer thickness, Cases A and B |
| `Hodge_John_ECE_5105_Final_Exam_P2a.m` | Graphical solution for the fundamental TM mode's β |
| `Hodge_John_ECE_5105_Final_Exam_P2b.m` | Field profile and Poynting vector of the guided mode |

## Running

Each script is standalone. Open MATLAB in this folder and run it by name, e.g. `Hodge_John_ECE_5105_Final_Exam_P1`. Each script starts with `clear all; close all`, so save your workspace first.

## Known issues and a modern version

A 2026 code review found one error in the original scripts, which are left as written in 2015.

- `Hodge_John_ECE_5105_Final_Exam_P2b.m`, lines 79 and 81: the core's Poynting-vector term divides by the cladding-1 permittivity instead of the core permittivity (making the core power 4.5 times too large), and the cladding-2 term decays with the cladding-1 constant κ<sub>x1</sub> and uses `real(A_cl_2).^2` where `abs(A_cl_2).^2` is needed.

**Modern version:** [Plane waves in multilayer media: the transfer-matrix method](https://github.com/jman4162/RF-and-Antenna-Tutorials/blob/main/Multilayer_Media_Transfer_Matrix_Method_Tutorial.ipynb) (Python notebook). It reproduces both exam problems, including the report's wavenumber table and β = 1.7971 × 10<sup>7</sup> m<sup>−1</sup>, with the Poynting vector corrected and checked against the interface conditions.

## Report

[Full exam report (PDF, 24 pages)](../../EM%20Boundaries%20-%20Hodge_John_ECE_5105_Final_Exam_Draft_A%20%281%29.pdf): derivations, k-vector tables, a total-internal-reflection analysis and the results.
