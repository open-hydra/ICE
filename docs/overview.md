---
title: Overview
---

# Overview

ICE (**I**nternal **C**ondensed-phase **E**quations) is an open-source Eulerian solver for dispersed condensed phases on multi-block structured grids, written in modern Fortran. It targets a wide range of particulate-phase problems — from mono-disperse aluminium particle suspensions to multi-group polydisperse flows with coupled heat transfer, drag, and phase change.

---

## Hydra CFD Suite

ICE is the **condensed-phase solver** of the **Hydra** CFD ecosystem — an integrated suite of tools for multi-physics simulation of complex systems.

| Component | Role | Status |
|-----------|------|--------|
| [**ATLAS**](https://github.com/open-hydra/ATLAS) | Pre-processor: mesh prep, initial & boundary conditions, material property data | Separate package |
| [**MOSE**](https://github.com/open-hydra/MOSE) | Solver: compressible Euler/Navier–Stokes with finite-rate chemistry | Sister solver |
| [**FUSS**](https://github.com/open-hydra/FUSS) | Solver: transient and steady-state heat conduction in solids | Sister solver |
| **ICE** | Solver: Eulerian condensed-phase equations for dispersed particles | This package |

!!! info "Using ICE without ATLAS"
    The input files required by ICE (initial conditions, boundary condition table, material property data) are **typically produced by ATLAS**. If ATLAS is not available, all input files can be prepared manually; see the [User Guide](user/using.md) for the expected formats.

---

## ICE Capabilities

ICE solves the Eulerian equations for the condensed (particulate) phase using the **monokinetic (MK) closure**, in which all particles within a group share the same velocity. The conservative-variable vector for each group is:

$$
\mathbf{U} = \begin{pmatrix} \alpha_p \rho_p \\ \alpha_p \rho_p u \\ \alpha_p \rho_p v \\ \alpha_p \rho_p w \\ \alpha_p \rho_p e_p \\ \alpha_p \end{pmatrix}
$$

where $\alpha_p$ is the particle volume fraction, $\rho_p$ the material density, $(u,v,w)$ the particle velocity, and $e_p$ the specific total energy. The governing equations include source terms for drag, heat exchange, and phase change.

---

## Physical Models

### Particle properties

ICE models each particle group as a condensed material with constant thermophysical properties.

| Property | Symbol | Description |
|----------|--------|-------------|
| Material density | $\rho_p$ | Mass per unit volume of particle material |
| Specific heat | $c_{s,p}$ | Constant-pressure heat capacity |
| Latent heat | $L_v$ | Latent heat of vaporisation |
| Heat of combustion | $q_p$ | Energy released per unit mass of reacted particle |

### Multi-group polydisperse model

A population of $N_\text{grp}$ independent groups is supported, each representing particles of a distinct size or material. Within each group, the monokinetic closure assumes all particles share the same local velocity. The group populations are specified at initialisation via `npop` and `ncond` arrays.

### Heat transfer

Gas–particle heat exchange is modelled via a Nusselt correlation:

$$
\dot{q}_{gp} = \frac{6 \alpha_p \lambda_g}{d_p^2}\, Nu(Re_p, Pr, Ma) \,(T_g - T_p)
$$

Available Nusselt models:

| Model | Conditions |
|-------|-----------|
| Stokes | $Re_p \ll 1$ |
| Ranz–Marshall | Moderate Reynolds number |

### Drag

Drag source terms are computed from the particle slip velocity and a drag coefficient $C_D(Re_p, Ma)$ selected from several built-in correlations.

### Radiation

Thermal radiation is modelled via the Stefan–Boltzmann law with configurable emissivity $\varepsilon$:

$$
\dot{q}_\text{rad} = \varepsilon \sigma (T_\text{ref}^4 - T_p^4)
$$

---

## Numerical Methods

### Spatial discretisation

| | Details |
|-|---------|
| Framework | Cell-centred finite volume on structured multi-block hexahedral grids |
| Convective flux | High-order reconstruction with flux limiters |
| Source terms | Point-implicit treatment for stiff drag and heat-transfer terms |
| Metric tensor | Computational-to-physical mapping evaluated per face |

### Time integration

| | Details |
|-|---------|
| Explicit scheme | Multi-stage Runge–Kutta (RK) |
| Stability | CFL condition on particle convection |
| Steady-state mode | Local time stepping |
| Time-accurate mode | Global minimum $\Delta t$ |

### Parallel computing

| Mode | Details |
|------|---------|
| Shared memory | OpenMP thread-level parallelism |
| Distributed memory | MPI domain decomposition across blocks (ghost-cell exchange) |
| Hybrid | OpenMP + MPI combined runs on HPC clusters |

---

## Code Dependencies

### Required libraries

| Library | Role | Source |
|---------|------|--------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O — Tecplot, VTK | Bundled submodule |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser | Bundled submodule |

### Optional libraries

| Library | Role |
|---------|------|
| OpenMP | Shared-memory thread parallelism |
| MPI | Distributed-memory parallelism |
| TecIO | Binary Tecplot output |

### Build toolchain

| Tool | Minimum version |
|------|----------------|
| CMake | 3.23 |
| Fortran compiler | GNU gfortran 11+ or Intel ifx/ifort |
| C / C++ compiler | GCC or ICC (for ORION and the optional TecIO) |

---

## Documentation Guide

| Section | What you'll find |
|---------|-----------------|
| [**Getting Started**](getting-started/index.md) | Installation, prerequisites, and first run |
| [**User Guide**](user/index.md) | Running simulations, configuring input files, boundary conditions, output |
| [**Theory Guide**](theory/index.md) | Governing equations, numerical methods, particle physics models |
| [**Verification & Validation**](vv/index.md) | Test cases and analytical benchmarks |
| [**Developer Guide**](development/index.md) | Repository architecture, testing framework, contribution guidelines |
| [**About**](about/index.md) | License, acknowledgements, and contributors |

---

## License

ICE is free and open-source software released under the **[GNU General Public License v3.0](about/license.md)** (GPL-3.0).

| Permission | |
|------------|-|
| :white_check_mark: Use freely | For any purpose, including commercial |
| :white_check_mark: Modify | Change the source code as needed |
| :white_check_mark: Distribute | Share original or modified versions |
| :white_check_mark: Patent grant | Contributors grant patent rights |
| :warning: Share-alike | Derivative works must use GPL-3.0 |
| :warning: Disclose source | Source code must be provided when distributing |

Full license text: [`LICENSE`](about/license.md)
