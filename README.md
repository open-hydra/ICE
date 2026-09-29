<p align="center">
  <h1 align="center">ICE</h1>
  <p align="center"><b>Integration of a Condensed phase via an Eulerian method</b></p>
</p>

<p align="center">
  <a href="https://open-hydra.github.io/ICE/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran-734f96.svg" alt="Language: Fortran">
  <a href="https://github.com/open-hydra/ICE/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
</p>

---

ICE is an open-source solver for a dispersed condensed phase — droplets or solid particles — on multi-block structured grids. Written in modern Fortran, it carries the particle cloud as a continuum rather than as individual parcels, and runs either on its own or one-way coupled to a frozen gas field produced by another solver. It is the condensed-phase component of the [Hydra](https://github.com/open-hydra) CFD suite.

## Features

- **Eulerian dispersed phase** — the cloud is carried as a bulk density and a number density, from which the particle radius follows; several particle families, each with its own closure, can be advanced in one run.
- **Three kinetic closures** — monokinetic (MK), isotropic Gaussian (IG) and anisotropic Gaussian (AG), trading cost against the ability to represent crossing trajectories and shear.
- **Interphase exchange** — drag, convective heat, grey-body radiation and evaporation, including a non-equilibrium interface and a Stefan-blowing correction.
- **High-order numerics** — MUSCL reconstruction with flux limiters and a Jameson shock sensor, Saurel/Rusanov/HLLE Riemann solvers, SSP Runge–Kutta time integration, implicit residual smoothing and grid sequencing.
- **Multi-block grids** — 1-D, 2-D or 3-D, joined by block connections or by chimera overset interpolation between overlapping non-matching blocks.
- **Parallel execution** — shared-memory parallelism via OpenMP; MPI over whole blocks for distributed-memory runs. The rank and thread counts do not change the answer.
- **Flexible I/O** — solution output in Tecplot (ASCII and binary) and VTK formats. Point probes for time-history recording. Restart capability.

## Quick Start

### Prerequisites

| Requirement | Details |
|---|---|
| **CMake** | ≥ 3.23 |
| **Fortran compiler** | GNU (`gfortran`) or Intel/oneAPI (`ifx`) |
| **C/C++ compiler** | Required by ORION, and by the optional TecIO component |

### Build

```bash
git clone --recurse-submodules https://github.com/open-hydra/ICE.git
cd ICE

# Build with GNU compilers and OpenMP
./install.sh build --compilers=gnu --use-openmp

# — or with Intel compilers and full feature set —
./install.sh build --compilers=intel --use-openmp --use-mpi --use-tecio
```

The executable is placed in `bin/ICE`.

See the [Installation Guide](https://open-hydra.github.io/ICE/getting-started/installation/) for all build options, CMake presets, and troubleshooting.

### Run the Crossing Jets

```bash
cd test/Doisneau/IG
./ICE.sh solve
```

See the [Quick Start](https://open-hydra.github.io/ICE/getting-started/quick-start/) for a full walkthrough.

## Dependencies

ICE is built on top of several companion libraries:

| Library | Role |
|---|---|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot, VTK, Plot3D) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |

Optional external libraries: **OpenMP**, **MPI**, **TecIO**.

The mesh, boundary-condition table and material property data are normally prepared by [ATLAS](https://github.com/open-hydra/ATLAS), the Hydra pre-processor. All of them are plain text and can be written by hand; the formats are documented in the user guide.

## Project Structure

```
ICE/
├── src/
│   ├── app/           # Main application and the registry doc generator
│   └── lib/           # Solver library sources
│       ├── base/        # Derived types and global state
│       ├── config/      # Input registry and INI parsing
│       ├── numerics/    # Fluxes, reconstruction, time integration, multigrid
│       ├── physics/     # Closures and interphase source terms
│       ├── io/          # Solution, probe and boundary-condition I/O
│       ├── parallel/    # MPI decomposition and ghost exchange
│       ├── driver/      # Setup, solve and post-processing stages
│       └── diagnostic/  # Residuals and run-time reporting
├── lib/               # Git submodule dependencies
├── test/              # Verification & validation cases
│   ├── Doisneau/      #   Crossing jets, one case per closure
│   ├── Berthon/       #   Riemann problems for the Gaussian closure
│   ├── verification/  #   Exact-solution cases, inputs generated on the fly
│   └── fast/          #   Short invariant checks
├── docs/              # MkDocs documentation source
├── cmake/             # CMake modules
├── install.sh         # Build helper script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/ICE](https://open-hydra.github.io/ICE/)**, covering:

- Installation & quick start
- User guide & input file reference
- Verification & validation results
- Theory guide (governing equations, numerical methods, particle physics models)

## License

ICE is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
