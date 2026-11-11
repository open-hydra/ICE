# Getting Started

ICE (*Integration of a Condensed phase via an Eulerian method*) solves a dispersed
condensed phase — droplets or solid particles — on multi-block structured grids. It
carries the cloud as a continuum: each cell holds fields that describe the particles
inside it, rather than tracking parcels one by one. It runs on its own, or one-way
coupled to a frozen gas field computed by another solver.

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Prerequisites, the build script and the CMake options.

    [:octicons-arrow-right-24: Install ICE](installation.md)

-   :material-rocket-launch:{ .lg .middle } __Quick Start__

    ---

    Run a shipped case, a few minutes on one thread, and check its result.

    [:octicons-arrow-right-24: Run the Quick Start](quick-start.md)

</div>

## What ICE solves

The cloud is described by two fields: the **bulk density** $\rho_p$, the mass of
particles per unit volume of mixture (kg/m³), and the **number density** $n$, the
number of particles per m³. The particle radius follows from the two,

$$
R_p = \left(\frac{3}{4\pi}\,\frac{\rho_p}{n\,\rho_{al}}\right)^{1/3},
$$

where $\rho_{al}$ is the density of the material itself. The size can therefore
change from one cell to the next with no extra bookkeeping.

Inside a cell, the particle velocities follow one of three **closures**:

| Closure | Variables | State |
|---|---|---|
| **MK** — monokinetic | 6 | $\rho_p,\ u,\ v,\ w,\ T_p,\ n$ |
| **IG** — isotropic Gaussian | 7 | $\rho_p,\ u,\ v,\ w,\ P,\ T_p,\ n$ |
| **AG** — anisotropic Gaussian | 12 | $\rho_p,\ u,\ v,\ w,\ P_{11},\ P_{12},\ P_{13},\ P_{22},\ P_{23},\ P_{33},\ T_p,\ n$ |

MK gives every particle of a cell the same velocity. It is the cheapest, but two
jets that meet cannot cross: they merge into a singular sheet. IG adds a scalar
velocity dispersion and AG a full dispersion tensor, which is what lets crossing
jets and sheared clouds be represented. See [Governing Equations](../theory/governing-equations.md).

A run can carry several **families** at once, each with its own closure. They share
the mesh and the gas field, and do not exchange mass or momentum with each other.

## Exchange with the gas

The source terms are active only when a gas field is given (`INPUT/gas.tec`, or its
binary and VTK forms). Without one, ICE transports the cloud with no exchange at all.

| Term | What it does |
|---|---|
| Drag | Momentum exchange with the gas. Twelve correlations. |
| Convective heat | Energy exchange through a Nusselt number. Seven correlations. |
| Radiation | Grey-body exchange with the local gas temperature. |
| Evaporation | Mass transfer to the gas. Five models, with an optional non-equilibrium interface and a Stefan-blowing correction. |
| Inside the drop | A uniform temperature or a parabolic profile, and an optional boiling branch above the boiling point. |
| Solidification | Freezing and melting of the material, with its latent heat. |

The material properties are constants, or tables against temperature in
`INPUT/part-properties.dat`. Several materials can be carried, each with its own
models. The correlations are given in [Particle Physics](../theory/physics.md);
combustion is not implemented.

## Numerical methods

| | |
|-|-|
| Discretisation | Cell-centred finite volumes on multi-block structured hexahedral grids |
| Reconstruction | First order, or MUSCL with ten limiters; an optional Jameson shock sensor |
| Riemann solver | Saurel, Rusanov or HLLE |
| Dimensions | 1-D, 2-D, axisymmetric or 3-D, from the block dimensions |
| Time integration | Explicit forward Euler, SSP-RK2 or SSP-RK3 |
| Time step | Global (time accurate) or per cell (steady state), from a CFL number |
| Steady-state aids | Implicit residual smoothing and grid sequencing |
| Parallel runs | OpenMP over the cells of a block, MPI over whole blocks, or both |

The number of MPI ranks does not change the answer: a run on four ranks is
bit-for-bit identical to a serial one. See [Spatial Discretization](../theory/numerics.md),
[Time Integration](../theory/time-integration.md) and [Parallel Execution](../user/parallel.md).

## The inputs of a case

ICE reads a mesh with its initial state, a table of boundary conditions and
`input.ini`. [ATLAS](https://github.com/open-hydra/ATLAS), the Hydra pre-processor,
writes the boundary table, the initial field and the material tables. All of them
are plain text, and can also be written by hand: see
[Initial Conditions](../user/initial-conditions.md) and
[Boundary Conditions](../user/boundary-conditions.md).
