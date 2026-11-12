# Edge-Emitting Laser Diode to Single-Mode Fiber Coupling Analysis

Optical modeling and engineering design of a coupling system between an 800 nm AlGaAs edge-emitting semiconductor laser diode and a 780HP single-mode optical fiber, developed for the *Passive Optical Components* course at Politecnico di Torino.

## Repository Structure

```text
laser-diode-smf-coupling/
├── README.md                  # Project overview, engineering rationale, and key results
├── src/
│   └── sm_wave.m              # Main MATLAB script (modal analysis, free-space & ABCD propagation)
├── docs/
│   └── report.pdf             # Technical report detailing derivations and design choices
└── results/                   # Exported MATLAB figures (.png)
    ├── beta0_search.png
    ├── normalized_mode.png
    ├── coupling_loss_vs_z.png
    ├── beam_vs_fiber_freespace.png
    └── beam_vs_fiber_lenses.png
```

## Project Overview

Coupling an edge-emitting laser diode into a single-mode fiber (SMF) presents a major optical challenge due to the strong **intrinsic astigmatism** of the laser source and the mode-field mismatch with the circular fiber core. 

This repository contains the MATLAB simulation pipeline (`src/sm_wave.m`) used to:
1. **Waveguide Modal Analysis:** Model the vertical laser structure as a 3-layer symmetric slab waveguide ($Al_{0.5}Ga_{0.5}As$ core with $n_{co} = 3.3067$, $Al_{0.7}Ga_{0.3}As$ cladding with $n_{clad} = 3.1710$) to determine the single-mode cutoff thickness ($d = 400\text{ nm} < 426.6\text{ nm}$), solve the eigenvalue equation ($\beta_0 = 25.5729\ \mu\text{m}^{-1}$), and extract the power-normalized fundamental TE mode profile $E_0(y)$ accounting for free-space impedance $Z_0 = 377\ \Omega$.
2. **Astigmatic Gaussian Beam Characterization:** Map the vertical waveguide mode to an equivalent Gaussian waist ($w_{0y} \approx 0.2905\ \mu\text{m}$) at the $1/e$ field amplitude point, while modeling the horizontal slow axis directly from its Mode Field Diameter ($w_{0x} = 25\ \mu\text{m}$).
3. **Direct Coupling (Butt-Coupling) Evaluation:** Compute analytical 2D overlap integrals accounting for both beam spot expansion $w(z)$ and spherical wavefront curvature $R(z)$ at practical packaging distances of $z = 100\ \mu\text{m}$ and $z = 1\text{ mm}$.
4. **Bulk Optics Coupling Design (ABCD Formalism):** Design a two-lens collimating/focusing system using commercial plano-convex N-BK7 lenses, incorporating paraxial ABCD matrix transformations and finite-aperture **beam clipping** losses.

---

## Engineering & Design Decisions

### 1. Fast-Axis vs. Slow-Axis Optimization Trade-off
An edge-emitting laser exhibits extreme diffraction along the vertical **fast axis** ($z_{Ry} \approx 0.3314\ \mu\text{m}$) compared to the well-collimated horizontal **slow axis** ($z_{Rx} \approx 2454.4\ \mu\text{m}$). 
* **Constraint:** Standard spherical plano-convex lenses are rotationally symmetric and apply the exact same linear magnification $M = f_2 / f_1$ to both transverse axes simultaneously.
* **Decision:** Rather than relying on expensive anamorphic/cylindrical optics, we prioritized preserving optical power along the rapidly diverging fast axis ($y$) at the expense of over-expanding the slow axis ($x$). Matching $w_{0y} \approx 0.2905\ \mu\text{m}$ to the 780HP fiber mode radius ($w_f = 2.5\ \mu\text{m}$) requires $M \approx 8.6$. We selected **$M = 10$** using commercial AR-coated N-BK7 plano-convex lenses:
  * **Collimating Lens ($L_1$):** $f_1 = 15\text{ mm}$, $D_1 = 12.7\text{ mm}$ (1/2-inch)
  * **Focusing Lens ($L_2$):** $f_2 = 150\text{ mm}$, $D_2 = 25.4\text{ mm}$ (1-inch)
  * **Intermediate Free-Space Path:** $d_2 = 50\text{ mm}$ (collimated flight region)

### 2. Physical Aperture Limits & Beam Clipping
Ideal ABCD matrix calculations assume infinite lens apertures. However, after propagating $d_1 = 15\text{ mm}$, the fast-axis beam radius expands to $w_{1y} \approx 13.15\text{ mm}$. 
* Standard optical design rules ($D \ge 3w$ for $>99\%$ power transmission) would require a $\sim 40\text{ mm}$ lens diameter with a $15\text{ mm}$ focal length—which is physically impossible due to hemispherical curvature limits ($D \le f_1$).
* By selecting a realistic $1/2\text{-inch}$ lens ($a_1 = 6.35\text{ mm}$), $L_1$ acts as the system's defining aperture stop, introducing a calculated **clipping loss of $4.29\text{ dB}$** ($37.3\%$ transmission). Because the beam is truncated to a maximum radius of $6.35\text{ mm}$, the larger $1\text{-inch}$ focusing lens $L_2$ ($a_2 = 12.7\text{ mm}$) introduces zero additional clipping.

### 3. Packaging & Alignment Tolerance Analysis
While spherical bulk lenses do not reduce total insertion loss compared to direct coupling (since the beam's elliptical aspect ratio is conserved), they are widely adopted in industrial packaging because they trade **longitudinal sensitivity** for **transverse sensitivity**:
* **Longitudinal ($z$-axis):** The collimated beam between $L_1$ and $L_2$ achieves a massive Rayleigh range ($z_{R1y} \approx 679\text{ m}$), relaxing axial placement tolerances from nanometers to millimeters.
* **Transverse ($x,y$-axes):** Because $M = 10$, any lateral misalignment at the laser diode is magnified tenfold at the fiber facet (e.g., a $1\ \mu\text{m}$ shift at the source becomes $10\ \mu\text{m}$ at the fiber, exceeding the $2.5\ \mu\text{m}$ core radius and dropping coupling efficiency to zero).

---

## Key Simulation Results

| Configuration | Horizontal Efficiency ($\eta_x$) | Vertical Efficiency ($\eta_y$) | Clipping Loss | Total Coupling Loss |
| :--- | :---: | :---: | :---: | :---: |
| **Direct Coupling ($z = 100\ \mu\text{m}$)** | — | — | $0\text{ dB}$ | **$39.21\text{ dB}$** ($\eta = 0.0120\%$) |
| **Direct Coupling ($z = 1\text{ mm}$)** | — | — | $0\text{ dB}$ | **$59.60\text{ dB}$** ($\eta = 0.0001\%$) |
| **Two-Lens System ($M = 10$)** | $0.04\%$ | $97.78\%$ | $4.29\text{ dB}$ | **$38.36\text{ dB}$** ($34.08\text{ dB}$ geom. + $4.29\text{ dB}$ clip.) |

---

## How to Run

1. Open MATLAB (tested on R2023b+).
2. Navigate to the `src/` directory and run:
   ```matlab
   sm_wave
   ```
3. The script prints the fitted vertical waist, coupling efficiencies, and system losses to the Command Window, and generates 5 diagnostic plots covering modal profiles, free-space loss evolution, and spatial beam footprints before and after the lens system.

## References
* H. Zhou, H. Xu, J.-a. Duan, *"Review of the technology of a single mode fiber coupling to a laser diode"*, Optical Fiber Technology 55 (2020) 102097.
* J. Alda, *"Laser and Gaussian Beam Propagation and Transformation"*, Encyclopedia of Optical Engineering (2003).
* A. M. Kowalevicz Jr., F. Bucholtz, *"Beam Divergence from an SMF-28 Optical Fiber"*, Naval Research Laboratory Report NRL/MR/5650--06-8996 (2006).
