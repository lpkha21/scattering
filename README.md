# scattering

MATLAB reproduction of

> H. Toyama, K. Yasumoto and T. Iwasaki, "Electromagnetic scattering from a dielectric cylinder with multiple eccentric cylindrical inclusions," *Progress In Electromagnetics Research (PIER)*, vol. 40, pp. 113–129, 2003.

A plane wave hits an infinitely long dielectric cylinder (the *host*) that contains smaller off-centre cylinders (the *inclusions*). The inclusions can be dielectric, air holes or perfect conductors. The code computes the aggregate T-matrix of the whole structure for TM (E<sub>z</sub>) and TE (H<sub>z</sub>) polarization. From it, it gets the backscattering, forward-scattering and differential scattering cross-sections. All results figures of the paper (Figs. 2–12) are reproduced.

<p align="center">
  <img src="results/8th%20figure.png" width="48%">
  <img src="results/10th%20figure.png" width="48%">
</p>

## Quick start

Requirements: MATLAB (R2020a or newer for `exportgraphics`). No toolboxes are needed. If the Parallel Computing Toolbox is installed, `parfor` loops in Figs. 2 and 3 run in parallel; otherwise they run as normal loops.

```matlab
cd('path/to/scattering')   % the folder containing the .m files
runAllFigures              % runs Figs. 2-12 and saves results/*.png and *.fig
```

Or run a single figure, e.g. `eighthFigure`.

## Using the solver directly

```matlab
M   = 12;                              % truncation: orders -M..M
L   = 2;                               % number of inclusions
ka  = [0.5, 0.5, 2.5];                 % k0*a  for inclusions 1..L, then the host
n   = [4, Inf, 2];                     % refractive indices (Inf = perfect conductor), host last
kd  = [1.5, 1.0];                      % k0*d  of each inclusion centre
th  = [0, pi];                         % angle of each inclusion centre (rad)

T   = T_matrix_L_mag(L, M, ka, n, kd, th);   % TM;  use T_matrix_L_El for TE
fin = 3*pi/2;                                 % direction the incident wave travels
fi0 = linspace(0, 2*pi, 361);                 % observation angles
plot(rad2deg(fi0), k0sigma(T, fin, fi0))      % k0*sigma_d
```

- Backscattering is `k0sigma(T, fin, fin + pi)`.
- Forward scattering is `k0sigma(T, fin, fin)`.
- With `L = 0` (pass `[]` for `kd` and `th`) you get a homogeneous cylinder.

All lengths are given as k<sub>0</sub>·length (k<sub>0</sub> = 1). All materials have μ = μ<sub>0</sub>. Angles are measured from the x-axis. `fin = 0` means the wave travels in +x.

## Files

| File | Purpose | Paper eq. |
|---|---|---|
| `T_matrix_L_mag.m`, `T_matrix_L_El.m` | Aggregate T-matrix T<sub>L+1</sub> of host + inclusions (TM / TE) | 7, 24–27 |
| `TmLMag.m`, `TmLEl.m` | Recursive T-matrix of each inclusion in the presence of the others | 9–10 |
| `T_matrix_base_case_mag.m`, `T_matrix_base_case_el.m` | T-matrix of one inclusion alone in the host medium | 12–13 |
| `alpha.m`, `beta.m` | Translation matrices between cylinder centres (Graf's addition theorem) | 7–11 |
| `besJ.m`, `besH.m`, `dbesJ.m`, `dbesH.m` | Diagonal matrices of Bessel/Hankel functions and their derivatives | 14–16 |
| `dij.m`, `fij.m` | Distance and angle between two cylinder centres | |
| `k0sigma.m` | Normalised cross-section k<sub>0</sub>σ<sub>d</sub> from a T-matrix | 28–31 |
| `secondFigure.m` … `twelfthFigure.m` | One script per paper figure | |
| `runAllFigures.m` | Runs every figure script and saves the outputs to `results/` | |
| `NoInclusion*.m`, `OneInclusion*.m` | Test scripts (optical theorem, reciprocity, simple cases) | 32–33 |

Naming: **Mag** = TM (E<sub>z</sub>), **El** = TE (H<sub>z</sub>).

## Method in brief

1. **Each inclusion alone.** Each inclusion gets a diagonal T-matrix from field matching on its surface. It is embedded in the host material, so the host wavenumber k<sub>h</sub> = n<sub>h</sub>k<sub>0</sub> is used.
2. **Multiple scattering.** The inclusions scatter waves onto each other. Graf's addition theorem moves waves between centres. Chew's recursive algorithm adds the inclusions one at a time. The result is one matrix T<sub>L</sub> for the whole cluster, seen from the host centre.
3. **Outer wall.** Matching fields on the host's outer boundary gives the aggregate T-matrix T<sub>L+1</sub>.
4. **Far field.** Using the large-argument Hankel asymptotics gives k<sub>0</sub>σ<sub>d</sub>(φ<sub>0</sub>) = 4 |Σ<sub>n</sub> (−i)<sup>n</sup> e<sup>inφ<sub>0</sub></sup> a<sub>n</sub>|², where a = T<sub>L+1</sub>·[i<sup>n</sup> e<sup>−inφ<sub>in</sub></sup>].

The only approximation is truncating the series at |n| ≤ M (M = 12, as in the paper).

## Figures

| Fig. | Pol. | Host (k<sub>0</sub>a, n) | Inclusions | Plotted |
|---|---|---|---|---|
| 2 | TE, TM | 3.0, n<sub>3</sub> = 2 − Δn | two, k<sub>0</sub>a = 0.25, n = 2 | k<sub>0</sub>σ<sub>b</sub> vs Δn |
| 3 | TM | 3.0, 1.5 | three, k<sub>0</sub>a = 0.25, n = 2; cylinder 1 moved radially | k<sub>0</sub>σ<sub>b</sub> vs k<sub>0</sub>d<sub>1</sub> |
| 4, 5 | TM, TE | 2.5, 2 | one or two, k<sub>0</sub>a = 0.5, n = 4 | k<sub>0</sub>σ<sub>b</sub>, k<sub>0</sub>σ<sub>f</sub> vs φ<sub>in</sub> |
| 6, 7 | TM, TE | 2.5, 2 | one or two perfect conductors | k<sub>0</sub>σ<sub>b</sub>, k<sub>0</sub>σ<sub>f</sub> vs φ<sub>in</sub> |
| 8 | both | 2.5, 2 | none | k<sub>0</sub>σ<sub>d</sub> vs φ<sub>0</sub>, φ<sub>in</sub> = 270° |
| 9–12 | both | 2.5, 2 | three, k<sub>0</sub>a = 0.3: air / n = 4 / conductor / mixed | k<sub>0</sub>σ<sub>d</sub> vs φ<sub>0</sub>, φ<sub>in</sub> = 270° |

## Verification

- All of Figs. 2–12 agree with the paper: peak values, lobe positions and curve shapes.
- The results also agree with an independent Python implementation of the same method.
- Reciprocity holds to about 10<sup>−15</sup>; see the check printed by `thirdFigure`.
- The optical theorem holds to about 10<sup>−14</sup>; see `OneInclusion6thfig`.

Benchmark values:

- **Fig. 8:** forward peak k<sub>0</sub>σ<sub>d</sub> = 114.7 (TM) and 61.2 (TE).
- **Fig. 2:** TM peak near Δn ≈ 0.43 is 10.4 (k<sub>0</sub>d<sub>1</sub> = 1.5) and 12.6 (k<sub>0</sub>d<sub>1</sub> = 2.25).

## Notes

- The paper's Eq. (33) prints the optical theorem as −(4/k<sub>0</sub>) Re{f}. Because f carries the factor (1−i)/√(πk<sub>0</sub>), the correct form is k<sub>0</sub>σ<sub>t</sub> = −4 Re Σ<sub>n</sub> (−i)<sup>n</sup> e<sup>inφ<sub>in</sub></sup> a<sub>n</sub>.
- Translation matrices must use the **host** wavenumber: argument n<sub>host</sub>·k<sub>0</sub>·d, not k<sub>0</sub>·d. Earlier versions of this code used k<sub>0</sub>·d, which made every figure with inclusions wrong except at Δn = 0 and Δn = 1. This is fixed in `T_matrix_L_mag.m` / `T_matrix_L_El.m`.
- Earlier versions also used the wrong incidence angle in Figs. 8–12 (120° instead of 270°). They also computed forward scattering as σ<sub>t</sub> − σ<sub>b</sub> in Figs. 4–7. Both are fixed.
