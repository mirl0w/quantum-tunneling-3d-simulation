# Quantum Tunneling & Atomic Orbitals — A 3D Julia/MATLAB Simulation

A physics computing project combining:
- **Julia**: numerically solves the 3D time-dependent Schrödinger equation
  (split-operator/FFT method) for a wave packet tunneling through a
  potential barrier, and analytically computes the exact hydrogen-atom
  orbitals (1s, 2s, 2p, 3d, 4f).
- **MATLAB**: renders both as publication-quality 3D isosurfaces, saves
  screenshots, and exports an MP4 animation of the tunneling event.

## Why this project
Quantum tunneling and atomic orbitals are two of the most iconic results
in physics — instantly recognizable, conceptually deep, and visually
striking. Solving the TDSE from scratch (not just plotting a textbook
formula) demonstrates real numerical methods (FFT-based PDE solving),
while the hydrogen orbitals are derived analytically from the actual
Laguerre/Legendre special functions rather than looked up — good signal
for a research-oriented application.

