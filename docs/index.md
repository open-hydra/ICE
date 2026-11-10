---
title: ICE
hide:
  - navigation
  - toc
---

<div class="ice-hero" markdown>

<div class="ice-hero__text" markdown>

<p class="ice-hero__eyebrow">HYDRA CFD SUITE</p>

<h1 class="ice-hero__title">ICE</h1>

<p class="ice-hero__lead">Open-source solver for droplets and solid particles on multi-block structured grids. The
cloud is a continuum, transported on its own or through a frozen gas field.</p>

<p class="ice-hero__cmd"><span>$</span> ./ICE.sh -m 4 -p 20 solve</p>

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[Verification](vv/index.md){ .md-button .ice-hero__ghost }

<p class="ice-hero__name"><strong>I</strong>ntegration of a <strong>C</strong>ondensed phase via an
<strong>E</strong>ulerian method</p>

</div>

<!-- The hero figure is docs/images/home-hero.png, drawn by extra.css as a background: without the file the
     hero shows only its own background. -->
<div class="ice-hero__shot" aria-hidden="true"><div class="ice-hero__figure"></div></div>

</div>

<p class="ice-kicker">The solver</p>

## A particle cloud as a continuum { .ice-heading }

<p class="ice-sub" markdown>ICE carries the cloud as fields on the mesh, not as individual parcels. Start with the
[Quick Start](getting-started/quick-start.md): a shipped case runs in a few minutes.</p>

<div class="ice-tools" markdown>

[:material-scatter-plot-outline: <b>Three closures</b> <i>Monokinetic, isotropic and anisotropic Gaussian: the Gaussian ones let crossing jets pass through each other.</i>](theory/governing-equations.md){ .ice-tool }
[:material-thermometer-water: <b>Exchange with the gas</b> <i>Drag, convective heat, radiation and evaporation, with several correlations for each.</i>](theory/physics.md){ .ice-tool }
[:material-grid: <b>Multi-block grids</b> <i>1-D, 2-D, axisymmetric or 3-D, joined by connections, periodic faces or overlapping chimera blocks.</i>](user/boundary-conditions.md){ .ice-tool }
[:material-chart-bell-curve: <b>Robust numerics</b> <i>Finite volumes with MUSCL reconstruction, three Riemann solvers and Runge-Kutta time stepping.</i>](theory/numerics.md){ .ice-tool }
[:material-chip: <b>Parallel runs</b> <i>OpenMP and MPI, alone or together; the answer does not depend on the number of ranks.</i>](user/parallel.md){ .ice-tool }
[:material-file-chart-outline: <b>Inputs and outputs</b> <i>Tecplot and VTK files, temperature-dependent material tables, probes and restarts.</i>](user/using.md){ .ice-tool }

</div>

<p class="ice-kicker"></p>

## Hydra CFD Suite { .ice-heading }

<p class="ice-sub">Hydra couples solvers of gases, real fluids, solids, particles and radiation in one
simulation. ICE is its Eulerian particle solver.</p>

<div class="hydra" markdown>

<div class="hydra__case" markdown>

<p class="hydra__label">Your case</p>

<p class="hydra__file" markdown>:material-grid: structured mesh</p>

<p class="hydra__file" markdown>:material-database: [databases](https://open-hydra.github.io/ATLAS/databases/)</p>

</div>

<div class="hydra__atlas" markdown>

<p class="hydra__name">ATLAS</p>

<p class="hydra__tools" markdown>[GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/) [BCB](https://open-hydra.github.io/ATLAS/user-guide/bcb/) [ICB](https://open-hydra.github.io/ATLAS/user-guide/icb/) [STB](https://open-hydra.github.io/ATLAS/user-guide/stb/) [MDB](https://open-hydra.github.io/ATLAS/user-guide/mdb/) [KAnT](https://open-hydra.github.io/ATLAS/user-guide/kant/)</p>

</div>

<div class="hydra__link">
<span class="hydra__link-line"></span>
<span class="hydra__link-bottom">input files for each run</span>
</div>

<div class="hydra__suite" markdown>

<p class="hydra__label" markdown>[Hydra](https://github.com/open-hydra) · the coupled simulation</p>

<div class="hydra__solvers" markdown>

[<b>MOSE</b> compressible flow <i>ideal-gas</i>](https://github.com/open-hydra/MOSE){ .hydra__solver .hydra__solver--gas }
[<b>GROOT</b> radiative transfer <i>no ATLAS phase</i>](https://github.com/open-hydra/GROOT){ .hydra__solver .hydra__solver--radiation }
[<b>ARES</b> real-fluid flow <i>real-fluid</i>](https://github.com/open-hydra/ARES){ .hydra__solver .hydra__solver--fluid }
[<b>FUSS</b> heat conduction in solids <i>solid</i>](https://github.com/open-hydra/FUSS){ .hydra__solver .hydra__solver--solid }
[<b>IGLOO</b> Lagrangian particles <i>dispersed</i>](https://github.com/open-hydra/IGLOO){ .hydra__solver .hydra__solver--particle }
[<b>ICE</b> Eulerian particles <i>dispersed</i>](https://github.com/open-hydra/ICE){ .hydra__solver .hydra__solver--particle .is-here }

</div>

</div>

<div class="hydra__libs" markdown>

<p class="hydra__label">Built on</p>

[<b>ORION</b> Tecplot, VTK and PLOT3D I/O](https://MarcoGrossi92.github.io/ORION/){ .hydra__lib }
[<b>FLINT</b> thermochemistry](https://MarcoGrossi92.github.io/FLINT/){ .hydra__lib }
[<b>OSLO</b> ODE integrators](https://github.com/MarcoGrossi92/OSLO){ .hydra__lib }
[<b>FiNeR</b> INI file parser](https://github.com/szaghi/FiNeR){ .hydra__lib }

</div>

</div>

<div class="ice-links" markdown>

[Quick Start](getting-started/quick-start.md)
[User Guide](user/index.md)
[Input Parameters](user/registry.md)
[Theory](theory/index.md)
[Verification](vv/index.md)
[Developer Guide](development/index.md)

</div>
