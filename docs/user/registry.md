# Input Parameters

Generated from the input registry by `bin/DocGen`; regenerate it after
changing any `reg%add` call. Every parameter is optional unless the
Required column says otherwise, and omitting one selects the default.

## [ICE-Parameters]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `newrun` | true | true, false |  no | Start a new simulation (false = restart) |
| `res-threshold` | 1e-10 | >= 0 |  no | Residual convergence threshold (0 = never stop on it) |
| `time-threshold` | 1e30 | > 0 |  no | Maximum simulation time |
| `iter-threshold` | 1000000000 | > 0 |  no | Maximum number of iterations |
| `cfl` | 0.5 | > 0 |  no | CFL stability parameter |
| `dt-max` | 1e-4 | > 0 |  no | Ceiling on the local time step [s], applied before the CFL factor: the step never exceeds cfl * dt-max |
| `cfl-rise-threshold` | 0 | >= 0 |  no | Ramp the CFL number linearly over this many iterations (0 = no ramp) |
| `time-accurate` | .true. |  |  no | Advance every cell with the global minimum step (true) or with its own local step, for steady state (false) |
| `irs` | .false. |  |  no | Enable implicit residual smoothing |
| `irs-beta` | 0.5 |  |  no | IRS Jacobi smoothing coefficient |

## [ICE-IO]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `sol-format` | tecplot ascii |  |  no | Solution output format: "tecplot ascii", "tecplot binary" (needs TecIO) or "vtk ascii" / "vtk binary" |
| `bck-format` | native binary |  |  no | Restart-file format, which also selects the reader for the initial condition: same choices as sol-format |
| `sol-diter` | 1000000000 | > 0 |  no | Solution output iteration frequency |
| `sol-dtime` | 1e30 | > 0 |  no | Solution output time frequency |
| `sol-overwrite` | true | true, false |  no | Overwrite solution files |
| `bck-diter` | 1000000000 | > 0 |  no | Backup output iteration frequency |
| `bck-dtime` | 1e30 | > 0 |  no | Backup output time frequency |
| `bck-overwrite` | true | true, false |  no | Overwrite backup files |
| `shell-diter` | 10 | > 0 |  no | Shell update iteration frequency |
| `res-diter` | 10 | > 0 |  no | Residual history write frequency |
| `ini-diter` | 1000000000 | > 0 |  no | Re-read input.ini every n iterations, so a running simulation can be re-steered |
| `gas-path` | INPUT/ |  |  no | Directory holding the one-way-coupling gas file gas.tec |

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

## [ICE-Scheme]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `space-reconstruction` |  |  |  no | Space reconstruction: MUSCL, MUSCL-SD (MUSCL with the density shock detector), or empty for first order |
| `flux-limiter` | none |  |  no | Flux limiter, used only with MUSCL: IORD, MINMOD, VANALBADA, VANLEER, OSPRE, UMIST, OSHER, SWEBY, MC, KOREN, SUPERBEE |
| `time` | 2 |  |  no | Time integrator: 1 = forward Euler, 2 = SSP-RK2, 3 = SSP-RK3 |
| `drag` | None |  |  no | Drag model, global for all families: Newton, Stokes, Schlichting, Schiller-Naumann, Wen-Yu, Putnam, Clift-Gauvin, Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, Hermsen |
| `heat` | None |  |  no | Heat transfer model, global for all families: Stokes, JAXA1, JAXA2, JAXA3, Chang, Ranz-Marshall, Kavanau-Drake |

## [ICE-Multigrid]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `levels` | 1 | > 0 |  no | Number of grid levels. Each coarse level halves every block dimension, so every block must be divisible by 2^(levels-1) |
| `level1-iter` | 1000000000 | > 0 |  no | Max iterations at multigrid level 1 |
| `level2-iter` | 1000000000 | > 0 |  no | Max iterations at multigrid level 2 |

## [ICE-Family1]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `model` |  |  | yes | Closure for this family: MK (monokinetic, 6 variables), IG (isotropic Gaussian, 7) or AG (anisotropic Gaussian, 12) |

## [ICE-Physics]

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| `rho` | 2700.0 | > 0 |  no | Condensed-material density [kg/m^3], used when no property table is given |
| `cs` | 1598.0 | > 0 |  no | Condensed-material specific heat [J/(kg K)], used when no property table is given |
| `lv` | 1.08e7 | > 0 |  no | Latent heat of vaporisation [J/kg]; only acts through the mass-transfer term |
| `q` | 9.53e6 | > 0 |  no | Heat of combustion [J/kg]; only acts through the mass-transfer term |
| `emiss` | 1.0 | >= 0 |  no | Particle surface emissivity; 0 switches radiative exchange off |
