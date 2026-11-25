# Whispering Gallery Modes in a Spherical Dielectric Resonator

*Virginia Tech ECE 5106 (EM Fields II) final exam, May 2016. Written in MATLAB.*

Finds the resonant whispering gallery modes (WGMs) of a 15 µm-radius dielectric microsphere in the near-infrared (λ ≈ 0.93–0.98 µm) using a spherical-harmonic field expansion.

![WGM field magnitude vs. radius](../../assets/wgm-field-profile.png)

*Normalized |E| for the L = M = 133 resonance at θ = π/2: the field is concentrated just inside the sphere's surface, which is the defining feature of a whispering gallery mode. Figure from the exam report.*

## What it does

The resonator is a sphere of radius a = 15 µm with core index n = 1.4 in a cladding of index n = 1.1, and the search considers TE modes with L = M. Inside the sphere the field is expanded in spherical Bessel functions j<sub>L</sub>(kr) and outside in spherical Hankel functions h<sub>L</sub><sup>(1)</sup>(kr), with vector spherical harmonics **X**<sub>LM</sub> carrying the angular dependence.

In Part A, matching the tangential fields at r = a gives a characteristic equation, and the scripts scan ω across the allowed band for each mode order L. A first scan finds a resonance near L = 133, ω ≈ 2.006 × 10<sup>15</sup> rad/s. The report then refines it by solving for the zero of a boundary-condition matrix determinant, which gives the final answer of L = 133, ω ≈ 2.0195 × 10<sup>15</sup> rad/s (λ ≈ 0.933 µm).

Part B computes the time-averaged electric and magnetic energy stored in the resonant mode. That step exists only in the report; its code isn't in this folder.

## Files

| File | Purpose |
|------|---------|
| `ECE_5106_FDTD_Final_Exam_A_L_search.m` | Part A: scans ω for a given L and plots both sides of the characteristic equation |
| `ECE_5106_FDTD_Final_Exam_Part_A_plot_fields.m` | Part A: plots the mode's field magnitude vs. r and φ |
| `ECE_5106_FDTD_Final_Exam_B_omega_search.m` | Finer frequency scan of the characteristic equation |
| `ECE_5106_FDTD_Final_Exam_B.m` | Evaluates the characteristic equation at a single chosen frequency |
| `plot_YLM.m` | Plots the components of the vector spherical harmonic **X**<sub>LM</sub> |
| `plot_sph_besselj.m` | Plots the spherical Bessel function j<sub>L</sub>(kr) inside the sphere |
| `Y_lm_func.m`, `Y_lm_func2.m`, `func_calc_Y_LM.m` | Scalar spherical harmonics Y<sub>LM</sub> (different normalizations and angle units) |
| `dY_lm_theta_func.m`, `func_calc_dY_LM_dtheta.m` | Numerical ∂Y<sub>LM</sub>/∂θ |
| `function_XLM_theta.m`, `function_XLM_phi.m` (+ `_mod` variants) | θ and φ components of **X**<sub>LM</sub> |
| `besselj_sph.m`, `function_spherical_bessel.m` | Spherical Bessel function j<sub>L</sub>; note the argument orders differ: `besselj_sph(L, x)` vs. `function_spherical_bessel(x, L)` |
| `ece_5106_final_exam_plot_trial_B.jpg` | Example output plot |

Despite "FDTD" in their names, these scripts solve the problem analytically; they are not FDTD simulations.

## Running

Open MATLAB in this folder and run a script by name. **Known issue:** the four `ECE_5106_FDTD_Final_Exam_*` scripts call a spherical Hankel function, `function_hankel_first`, that wasn't in the original upload. The version in [`../multipole-scattering`](../multipole-scattering/function_hankel_first.m) takes its arguments as `(x, L)`, but the characteristic-equation lines here call it as `(L, x)`. These scripts need that reconciled before they will run correctly. `plot_YLM.m` and `plot_sph_besselj.m` don't use that function.

## Known issues

Besides the missing `function_hankel_first` described above, a 2026 code review found two errors. The code is left as written in 2016; the report itself uses the correct material values.

- All five scripts that define the materials set `eps_r = sqrt(n)` where it should be `n^2`. The wavenumbers use n directly and are right, but the wave impedances Z<sub>co</sub> and Z<sub>cl</sub> come out as Z<sub>0</sub>/n<sup>1/4</sup> instead of Z<sub>0</sub>/n, which shifts the characteristic-equation roots. The report uses the correct ε<sub>r</sub> = 1.96 and Z<sub>co</sub> = 269.09 Ω.
- `ECE_5106_FDTD_Final_Exam_B.m` and `_B_omega_search.m` (line 50) compute `XLM_phi` with `function_XLM_theta`, and `_Part_A_plot_fields.m` (line 69) uses the θ function for the φ component, so E<sub>φ</sub> duplicates E<sub>θ</sub> in the field plots.

## Report

[Full exam report (PDF, 23 pages)](../../Whispering%20Gallery%20Mode%20-%20Hodge_John_ECE_5106_Final_Exam_Doc_Draft_E.pdf): derivation of the characteristic equation, the mode search, field plots and the stored-energy calculation.
