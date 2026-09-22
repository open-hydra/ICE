---
title: Overview
---

# Overview

ICE (*Integration of a Condensed phase via an Eulerian method*) is an
open-source solver for a dispersed condensed phase — droplets or solid particles —
on multi-block structured grids, written in modern Fortran. It carries the particle
cloud as a continuum rather than as individual parcels, and can run on its own or
against a frozen gas field produced by another solver.

---

## Hydra CFD Suite

ICE is the **condensed-phase solver** of the **Hydra** CFD ecosystem — an integrated
suite of tools for multi-physics simulation of complex systems.

| Component | Role | Status |
|-----------|------|--------|
| [**ATLAS**](https://github.com/open-hydra/ATLAS) | Pre-processor: mesh prep, initial & boundary conditions, material property data | Separate package |
| [**MOSE**](https://github.com/open-hydra/MOSE) | Solver: compressible Euler/Navier–Stokes with finite-rate chemistry | Sister solver |
| [**FUSS**](https://github.com/open-hydra/FUSS) | Solver: transient and steady-state heat conduction in solids | Sister solver |
| **ICE** | Solver: Eulerian condensed-phase equations for dispersed particles | This package |

!!! info "Using ICE without ATLAS"
    The boundary-condition table ICE reads is normally written by ATLAS BCB, and the
    property table by ATLAS. Both are plain text and can be written by hand; the
    formats are given under [Boundary Conditions](user/boundary-conditions.md) and
    [Initial Conditions](user/initial-conditions.md).

---

## What ICE solves

The cloud is described by a bulk density $\rho_p$ (mass of condensed material per unit
volume of *mixture*, kg/m³) and a number density $n$ (particles per m³). The particle
radius is not a state variable: it follows from the two,

$$
R_p = \left(\frac{3}{4\pi}\,\frac{\rho_p}{n\,\rho_{al}}\right)^{1/3},
$$

with $\rho_{al}$ the density of the condensed material itself. One consequence is that
a cell can hold a different particle size from its neighbour without any extra
bookkeeping.

ICE offers three closures for the velocity distribution inside a cell, chosen per
family:

| Closure | Variables | State |
|---|---|---|
| **MK** — monokinetic | 6 | $\rho_p,\ u,\ v,\ w,\ T_p,\ n$ |
| **IG** — isotropic Gaussian | 7 | $\rho_p,\ u,\ v,\ w,\ P,\ T_p,\ n$ |
| **AG** — anisotropic Gaussian | 12 | $\rho_p,\ u,\ v,\ w,\ P_{11},\ P_{12},\ P_{13},\ P_{22},\ P_{23},\ P_{33},\ T_p,\ n$ |

MK assumes every particle in a cell moves at the same velocity, which makes the system
pressureless: two clouds meeting at an angle pass through each other only if the
scheme lets them, and their trajectories cross in a singularity. IG adds one scalar
velocity dispersion, AG the full symmetric dispersion tensor, which is what lets a
crossing or a shear be represented. See
[Governing Equations](theory/governing-equations.md).

Several families can be carried at once, each with its own closure, set by one
`[ICE-FamilyN]` section per family. They share the mesh and the gas field but not
their state, and do not exchange mass or momentum with one another.

---

## Physical models

Source terms are active only when a gas field is present — that is, when
`INPUT/gas.tec` exists. Without it ICE runs 0-way coupled and transports the cloud
with no exchange at all.

| Term | What it does |
|---|---|
| Drag | Momentum and the matching kinetic-energy exchange with the gas, through a relaxation time built from $C_d$. Twelve correlations. |
| Convective heat | Energy exchange through a Nusselt number. Seven correlations. |
| Radiation | Grey-body exchange with the local gas temperature, at emissivity `emiss`. |

The correlations and their expressions are listed under
[Particle Physics](theory/physics.md). Mass transfer between the phases
(vaporisation, combustion) is present in the equations as a term but is not evaluated:
the closures carry it, the source routine sets it to zero.

The material density and specific heat can be constants from `[ICE-Physics]` or
tabulated against temperature in `INPUT/part-properties.dat`; the table wins when it is
present, and ICE says at startup which of the two it is using.

---

## Numerical methods

| | Details |
|-|---------|
| Framework | Cell-centred finite volume on multi-block structured hexahedral grids, two ghost layers per face |
| Reconstruction | First order, or MUSCL with a choice of eleven limiters; MUSCL-SD blends back towards first order near a shock through a Jameson-type density sensor |
| Riemann solver | Chosen from the closure: an upwind flux for MK (pressureless), Rusanov for IG and AG |
| Dimensionality | 1-D, 2-D or 3-D, inferred from the block dimensions |
| Time integration | Explicit forward Euler, SSP-RK2 or SSP-RK3 |
| Time step | Global minimum (time accurate) or per cell (steady state), from a CFL condition capped by `dt-max` |
| Convergence aids | Implicit residual smoothing; grid sequencing from a coarse level up |
| Source terms | Evaluated explicitly, inside the Runge-Kutta stages |

Details in [Spatial Discretization](theory/numerics.md) and
[Time Integration](theory/time-integration.md).

### Parallel computing

| Mode | Details |
|------|---------|
| Shared memory | OpenMP over the cells of each block |
| Distributed memory | MPI over whole blocks: each rank updates the blocks it owns and exchanges the interior cells its neighbours' ghosts read |
| Hybrid | Both together, threads inside each rank |

The number of ranks does not change the answer: a run on four ranks is bit-for-bit
identical to a serial one. See [Multi-block and MPI](vv/multiblock-mpi.md).

---

## Dependencies

### Required

| Library | Role | Source |
|---------|------|--------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Mesh and solution I/O — Tecplot, VTK | Bundled submodule |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser | Bundled submodule |

### Optional

| Library | Role |
|---------|------|
| OpenMP | Shared-memory thread parallelism |
| MPI | Distributed-memory parallelism |
| TecIO | Binary Tecplot output, through ORION |

### Build toolchain

| Tool | Requirement |
|------|-------------|
| CMake | 3.23 or newer |
| Fortran compiler | GNU `gfortran` or Intel `ifx` |
| C / C++ compiler | Needed by ORION and, if enabled, TecIO |

---

## Documentation guide

| Section | What you'll find |
|---------|-----------------|
| [**Getting Started**](getting-started/index.md) | Installation, prerequisites, and first run |
| [**User Guide**](user/index.md) | Case layout, `input.ini`, initial and boundary conditions, output |
| [**Theory Guide**](theory/index.md) | Governing equations, numerical methods, particle physics models |
| [**Verification & Validation**](vv/index.md) | Test cases, exact solutions and analytical benchmarks |
| [**Developer Guide**](development/index.md) | Repository architecture, testing framework, contribution guidelines |
| [**About**](about/index.md) | License, acknowledgements, and contributors |

---

## License

ICE is free and open-source software released under the **[GNU General Public License
v3.0](about/license.md)** (GPL-3.0).

| Permission | |
|------------|-|
| :white_check_mark: Use freely | For any purpose, including commercial |
| :white_check_mark: Modify | Change the source code as needed |
| :white_check_mark: Distribute | Share original or modified versions |
| :white_check_mark: Patent grant | Contributors grant patent rights |
| :warning: Share-alike | Derivative works must use GPL-3.0 |
| :warning: Disclose source | Source code must be provided when distributing |

Full license text: [`LICENSE`](about/license.md)
