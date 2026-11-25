# EM Scattering from a Sphere Using a Multipole (Vector Spherical Harmonic) Expansion

*Virginia Tech ECE 5106 (EM Fields II) homework set 5, spring 2016. Written in MATLAB.*

Expands an incident plane wave in vector spherical harmonics, solves for the multipole scattering coefficients of a perfectly conducting sphere, and maps the total field around it.

![Total field magnitude around the sphere](plots/ECE_5106_HW5_P5c_Take1_lookup_pi_by_pi_fig.jpg)

*Total electric field magnitude |E<sub>tot</sub>| (dB) around a sphere of radius a = 0.6 m with k = 1 m<sup>−1</sup> (ka = 0.6). The field is set to zero inside the sphere.*

## What it does

- **Problem 4: plane-wave expansion.** Writes the incident plane wave as a sum over L of spherical Bessel functions j<sub>L</sub>(kr) times vector spherical harmonics **X**<sub>LM</sub>. It checks convergence by rebuilding Re[**e**<sub>x</sub> · **E**<sub>inc</sub>] over a grid as more L terms are added.
- **Problem 5a: scattering coefficients.** Applies the conducting-sphere boundary conditions at r = a to get the multipole coefficients α<sub>L</sub> = −2 j<sub>L</sub>(ka) / h<sub>L</sub><sup>(1)</sup>(ka) and β<sub>L</sub> = −2 [x j<sub>L</sub>(x)]′ / [x h<sub>L</sub><sup>(1)</sup>(x)]′, evaluated at x = ka.
- **Problem 5b: scattering cross-section.** σ<sub>s</sub> = π/(2k²) Σ (2L+1)(|α<sub>L</sub>|² + |β<sub>L</sub>|²), with a plot of each L-term's contribution.
- **Problem 5c: total field map.** Sums the incident and scattered multipole series to map |E<sub>tot</sub>| on a grid around the sphere (figure above).

Radial derivatives of the Bessel and Hankel functions, and angular derivatives of Y<sub>LM</sub>, are computed with central finite differences (step h = 10<sup>−5</sup>, or 10<sup>−4</sup> in `ECE_5106_HW5_P4a.m`).

## Files

| File | Purpose |
|------|---------|
| `ECE_5106_HW5_P4a.m` | Problem 4a, first attempt: computes Y<sub>LM</sub> and ∂Y<sub>LM</sub>/∂θ at one point, then stops (incomplete) |
| `ECE_5106_HW5_P4a_take2.m` | Problem 4a: plane-wave expansion on an x–z grid |
| `ECE_5106_HW5_P4a_take2b.m`, `ECE_5106_HW5_P4b.m` | Problem 4: Re[**e**<sub>x</sub> · **E**<sub>inc</sub>] maps on the x–z and x–y planes |
| `ECE_5106_HW5_P5a.m`, `ECE_5106_HW5_P5a_Take2.m` | Problem 5a: α<sub>L</sub> and β<sub>L</sub> coefficients (two versions) |
| `ECE_5106_HW5_P5b.m`, `ECE_5106_HW5_P5b_Take2.m` | Problem 5b: scattering cross-section (two versions; `_Take2` loads `L_alpha_beta_coeff_arr.mat`) |
| `ECE_5106_HW5_P5c_Take1.m` | Problem 5c: total field map |
| `L_alpha_beta_coeff_arr.mat` | Saved α<sub>L</sub>, β<sub>L</sub> table used by `ECE_5106_HW5_P5b_Take2.m` |
| `ECE_5106_HW5_P4b.fig` | MATLAB figure file for Problem 4b |
| `function_spherical_bessel.m`, `besselj_sph.m` | Spherical Bessel function j<sub>L</sub>; note the argument orders differ: `function_spherical_bessel(x, L)` vs. `besselj_sph(L, x)` |
| `function_hankel_first.m` | Spherical Hankel function h<sub>L</sub><sup>(1)</sup>, called as `(x, L)` |
| `function_calc_f_L_kr.m`, `function_calc_f_L_kr_B.m`, `function_calc_g_L_kr.m` | Radial functions f<sub>L</sub>(kr) and g<sub>L</sub>(kr) for the total field |
| `func_calc_Y_LM.m`, `Y_lm_func2.m` | Scalar spherical harmonics Y<sub>LM</sub> |
| `func_calc_dY_LM_dtheta.m`, `func_calc_dY_LM_dphi.m` | Numerical ∂Y<sub>LM</sub>/∂θ and ∂Y<sub>LM</sub>/∂φ |
| `function_XLM_theta.m`, `function_XLM_phi.m` | θ and φ components of **X**<sub>LM</sub> |
| `plots/` | Output figures for Problems 4b and 5c |

## Running

Open MATLAB in this folder and run a script by name, e.g. `ECE_5106_HW5_P5c_Take1`. Every helper function the scripts call is in this folder. Each script starts with `clear all`, so save your workspace first.

## Known issues and a modern version

A 2026 code review found two errors in `ECE_5106_HW5_P5c_Take1.m`, the script behind the total-field map. The code is left as written in 2016.

- The grid coordinate `x` (line 30) overwrites `x = k*r` (line 16), so on the first pass of the L loop the L = 1 terms (lines 70 and 75) are evaluated at the Cartesian coordinate instead of kr.
- `function_calc_f_L_kr.m` recomputes β<sub>L</sub> at the observation radius instead of at the sphere radius r = a, so the radial-derivative terms use the wrong coefficient everywhere outside the sphere.

The coefficient and cross-section formulas used in `P5a` and `P5b` are equivalent to the standard Mie result; the notebook below shows the match.

**Modern version:** [Scattering from a conducting sphere: the Mie series](https://github.com/jman4162/RF-and-Antenna-Tutorials/blob/main/Scattering_From_a_Conducting_Sphere_Mie_Series_Tutorial.ipynb) (Python notebook). It computes the field map for the same ka = 0.6 case and checks the plane-wave expansion, the boundary condition on the sphere, the optical theorem and the Rayleigh and optical limits.

## Report

There is no standalone report for this homework set; the figures in [`plots/`](plots/) are its outputs.
