# Input Parameters

Generated from the input registry by `bin/DocGen`; regenerate it after
changing any `reg%add` call. Every parameter is optional unless the
Required column says otherwise, and omitting one selects the default.

## [ICE-Parameters]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `phase` |  |  |  no | ATLAS name of the condensed phase to read: INPUT/<phase>-{bc.txt,ic,properties.dat}. Absent = keep the prefix in force (standalone default part-) |
| `newrun` | true | true, false |  no | Start a new simulation (false = restart) |
| `res-threshold` | 1e-10 | >= 0 |  no | Residual convergence threshold (0 = never stop on it) |
| `time-threshold` | 1e30 | > 0 |  no | Maximum simulation time |
| `iter-threshold` | 1000000000 | > 0 |  no | Maximum number of iterations |

## [ICE-IO]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `ic-format` | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk binary, vtk raw |  no | Initial condition (INPUT/part-ic.*) format |
| `sol-format` | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk binary, vtk raw |  no | Solution (OUTPUT/part-field.*) format, also read back on restart |
| `sol-diter` | 1000000000 | > 0 |  no | Solution output iter frequency |
| `sol-dtime` | 1e30 | > 0 |  no | Solution output time frequency |
| `sol-overwrite` | true | true, false |  no | Overwrite solution files |
| `shell-diter` | 10 | > 0 |  no | Shell update iter frequency |
| `res-diter` | 10 | > 0 |  no | Residual history iter frequency |
| `ini-diter` | 1000000000 | > 0 |  no | input.ini update iter frequency |
| `gas-path` | INPUT/ |  |  no | Directory holding the coupling gas file gas.tec |

## [ICE-Probes]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `probe1` |  |  |  no | Name of the section configuring this probe; it also names its output file OUTPUT/<name>.txt |

## [probe-section]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `variables` | none |  |  no | Probe variables to write |
| `dtime` | 1e30 | > 0 |  no | Probe output time frequency |
| `diter` | 1000000000 | > 0 |  no | Probe output iter frequency |
| `index-position` | 0 0 0 0 | >= 0 |  no | Probe location by index (b i j k) |
| `position` | 0.0 0.0 0.0 |  |  no | Probe location by coordinates |

## [ICE-Numerics]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `time-scheme` | RK2 | euler, RK2, RK3 | yes | Time integration solver |
| `cfl` | 0.5 | > 0 | yes | CFL number |
| `dt-max` | 1e-4 | > 0 |  no | Ceiling on the time step [s], applied after the CFL factor |
| `cfl-rise-threshold` | 0 | >= 0 |  no | CFL rise threshold |
| `time-accurate` | .true. | logical | yes | Time accurate switch |
| `irs` | .false. | logical |  no | Implicit Residual Smoothing |
| `irs-beta` | 0.5 | >= 0 |  no | IRS beta parameter |
| `space-reconstruction` | first-order | MUSCL, first-order | yes | Space reconstruction method |
| `flux-limiter` | none | minmod, vanalbada, vanleer, ospre, umist, osher, sweby, mc, koren, superbee, none |  no | Flux limiter for space reconstruction |
| `shock-detector` | none | Jameson, none |  no | Shock detector method |
| `riemann-solver` |  | Saurel, Rusanov, HLLE |  no | Riemann solver (empty: Saurel for MK, Rusanov for IG and AG) |

## [ICE-Multigrid]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `levels` | 1 | > 0 |  no | Number of grid levels; every block dimension must be divisible by 2^(levels-1) |
| `level1-iter` | 1000000000 | > 0 |  no | Iterations for multigrid level 1 |
| `level2-iter` | 1000000000 | > 0 |  no | Iterations for multigrid level 2 |

## [ICE-Physics]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `drag` | none | Newton, Stokes, Schlichting, Schiller-Naumann, Wen-Yu, Putnam, Clift-Gauvin, Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, Hermsen, NoDrag, none |  no | Drag model, global for all families; required for a coupled run (none = not set, NoDrag = no momentum exchange) |
| `heat-transfer` | none | Stokes, JAXA1, JAXA2, JAXA3, JAXA4, Chang, Ranz-Marshall, Kavanau-Drake, NoHeat, none |  no | Convective heat transfer model, global for all families; required for a coupled run (none = not set, NoHeat = no convective exchange; Chang stops with a pointer to JAXA3, which is its formula) |
| `evaporation` | none | d2-law, CEM, CEM-B, ASM, TC, none |  no | Evaporation model, global for all families |
| `evaporation-interface` | VLE | VLE, LK |  no | Vapour-liquid interface: VLE equilibrium, or LK Langmuir-Knudsen non-equilibrium; ignored when evaporation is none |
| `evaporation-blowing` | none | LK, none |  no | Stefan-blowing reduction of the convective heat; LK applies Miller-Harstad-Bellan f2. Ignored under ASM and TC, which carry their own gas-side heat |
| `density` | 2700.0 | > 0 |  no | Condensed-material density [kg/m^3], used when no property table is given |
| `specific-heat` | 1598.0 | > 0 |  no | Condensed-material specific heat [J/(kg K)], used when no property table is given |
| `latent-heat` | 1.08e7 | > 0 |  no | Latent heat of vaporisation [J/kg]; the evaporation models use it both as the energy sink and as the anchor of the saturation curve |
| `emissivity` | 1.0 | >= 0 |  no | Particle surface emissivity; 0 switches radiative exchange off |
| `vapour-molar-mass` | 26.98 | > 0 |  no | Molar mass of the vapour [kg/kmol] |
| `boiling-temperature` | 2792.0 | > 0 |  no | Boiling temperature at 1 atm [K], the anchor of the Clausius-Clapeyron saturation pressure |
| `vapour-specific-heat` | 0.0 | >= 0 |  no | Specific heat of the vapour [J/(kg K)]; 0 falls back to the gas cp |
| `lewis-number` | 1.0 | > 0 |  no | Lewis number of the vapour in the gas, Le = k/(rho cp D) |
| `vapour-mass-fraction` | 0.0 | >= 0 |  no | Vapour mass fraction in the far-field gas; evaporation stops once the surface value falls to it |
| `evaporation-coefficient` | 1.0 | > 0 |  no | Evaporation (accommodation) coefficient of the Langmuir-Knudsen interface; unused under VLE |

## [ICE-Family1]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `closure` |  | MK, IG, AG | yes | Kinetic closure for this family |
