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

## Folder structure
```
quantum3d_project/
├── julia/
│   ├── Project.toml
│   ├── tunneling_solver.jl      -> produces tunneling_frames.mat
│   └── hydrogen_orbitals.jl     -> produces orbitals.mat
├── matlab/
│   ├── visualize_tunneling.m    -> reads tunneling_frames.mat, makes PNGs + MP4
│   └── visualize_orbitals.m     -> reads orbitals.mat, makes PNGs + gallery
└── outputs/                     -> created automatically, holds all screenshots/video
```

## Step 1 — Install Julia
Download from https://julialang.org/downloads (free). Then, in a terminal:
```bash
cd quantum3d_project/julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
This installs FFTW.jl and MAT.jl (only needed once).

## Step 2 — Run the physics simulations
```bash
julia --project=. tunneling_solver.jl
julia --project=. hydrogen_orbitals.jl
```
- `tunneling_solver.jl` takes ~1-3 minutes on a laptop (64³ grid, 800 steps).
  It prints the wave packet's mean energy, the barrier height, and the
  final tunneling transmission probability — a real number you can quote
  in your write-up (e.g. "12% transmission despite E < V0").
- `hydrogen_orbitals.jl` takes well under a minute.

Both scripts write `.mat` files directly into the `julia/` folder.

## Step 3 — Run the MATLAB visualizations
Open MATLAB, `cd` into `matlab/`, then:
```matlab
visualize_tunneling.m
visualize_orbitals.m
```
This populates `outputs/` with:
- `tunneling_before_barrier.png`, `tunneling_at_barrier.png`,
  `tunneling_tunneling.png`, `tunneling_after_barrier.png`
- `tunneling_animation.mp4`
- `orbital_1s.png`, `orbital_2s.png`, `orbital_2pz.png`,
  `orbital_3dz2.png`, `orbital_4fz3.png`
- `orbital_gallery.png` (a single poster-style composite — great as a
  standalone image for a portfolio page or application)

No manual screenshotting needed — `exportgraphics` renders these at
300 DPI automatically, so they're print-quality out of the box.

## Tuning knobs worth knowing (for your write-up / if asked about it)
- `tunneling_solver.jl`: `V0` (barrier height), `kx0` (packet momentum,
  sets energy E = kx0²/2), `barrier_half_width` — raising `V0` above `E`
  further shows the "classically forbidden" regime; try a few values and
  plot transmission probability vs. barrier height/width for an extra
  figure (this reproduces the textbook tunneling-coefficient curve).
- `hydrogen_orbitals.jl`: swap in other `(n,l,m)` triples for other
  orbitals (e.g. `(3,1,0)` for 3p, `(3,2,2)` for a cloverleaf 3d_{x²-y²}
  analog) — the code is fully general, not hardcoded to these five.

## Suggested framing for an application
Lead with the physics question ("does the electron ever go where
classical mechanics forbids?"), show the barrier + wave packet figure,
state the numerical transmission probability from your run, then show
the orbital gallery as the "structure of the atom" payoff. One paragraph
of physics + one number + one striking image tends to land better than
a long code walkthrough.
