# Code Structure

## Repository layout

```
ICE/
├── bin/                        # ICE and DocGen (generated)
├── cmake/                      # Compiler-flag probing and Find modules
├── docs/                       # This documentation (MkDocs)
├── lib/
│   ├── ORION/                  # Submodule: mesh and solution I/O
│   └── third_party/FiNeR/      # Submodule: INI parser
├── src/
│   ├── app/
│   │   ├── main.f90            # The solver executable
│   │   └── docgen.f90          # Regenerates docs/user/registry.md
│   └── lib/                    # Everything else, built into libICEL.a
├── test/                       # Cases, wired into CTest
├── .githooks/pre-push          # Runs CTest before a push
├── CMakeLists.txt
├── install.sh
└── mkdocs.yml
```

## Library layout

```
src/lib/
├── base/
│   ├── Parameters_m.f90        # Kinds, string lengths, pi, Stefan-Boltzmann
│   ├── Base_Types_m.f90        # Vectors, faces, metric containers
│   ├── Advanced_Types_m.f90    # Block, domain, BC and simulation types
│   ├── Global_m.f90            # ngroups, ncond, ghost-layer count, phase prefix
│   └── Series_Data_m.f90
├── config/
│   ├── Registry.f90            # The input registry: add, validate, emit markdown
│   ├── Register_*.f90          # One per section group; the source of truth for inputs
│   ├── Backend_INI.f90         # FiNeR wrapper; scans for families, probes, MG levels
│   ├── Read_Ini.f90            # Scan, register, load, validate
│   ├── Assign_Setup.f90        # Post-read derivations: closures, ncond, coupling
│   └── Config_Types_m.f90      # The obj_* configuration objects
├── driver/
│   ├── Procedures_m.f90        # The ICE_type facade: setup / solve / postprocess
│   ├── Wrap_Setup.f90          # Read, allocate, metrics, BCs, partition, gas
│   ├── Wrap_Solve.f90          # One step, plus the grid-level transition
│   ├── Wrap_Postprocess.f90    # Output frequencies, shell reporting, restart files
│   ├── Mod_Allocate_Data.f90
│   └── Mod_Phase.f90           # Maps the ORION arrays onto the phase state
├── io/
│   ├── IO_Solution.f90         # Solution and restart read/write via ORION
│   ├── IO_BC.f90               # Parses the ATLAS boundary table
│   ├── IO_Probes.f90
│   └── Load_Table.f90          # Reads and checks the optional property table
├── numerics/
│   ├── space/
│   │   ├── Mod_Metrics.f90         # Areas, normals, volumes, dimensionality
│   │   ├── Lib_Ghost.f90           # Ghost fill for every boundary type
│   │   ├── Lib_Reconstruction.f90  # Limited piecewise-linear states
│   │   ├── Lib_Limiters.f90        # Eleven limiters
│   │   └── Lib_Shock_Detector.f90  # The Jameson shock detector
│   ├── fluxes/
│   │   ├── Mod_Fluxes.f90          # Interior faces, two-pass to avoid races
│   │   ├── bc/Mod_BC_Fluxes.f90    # Boundary faces, from the ghost values
│   │   └── riemann/Lib_Riemann.f90 # Saurel, Rusanov, HLLE
│   ├── time/
│   │   ├── Mod_dt.f90              # CFL step, local or global
│   │   └── explicit/               # RK stages, residual, IRS, state update
│   └── multigrid/                  # Coarse-grid construction and prolongation
├── parallel/
│   ├── Mod_MPI.f90             # Environment, block partition, collectives
│   └── Mod_GhostExchange.f90   # Halo schedule and persistent requests
├── physics/
│   ├── Lib_Model_MK/IG/AG.f90  # The three closures
│   ├── Lib_Model.f90           # Binds the model procedures to the active family
│   ├── Lib_Drag.f90            # Twelve drag correlations
│   ├── Lib_Heat.f90            # Seven Nusselt correlations
│   ├── Lib_Evaporation.f90     # Five evaporation models and their interface options
│   ├── Lib_Properties.f90      # Density, cp and energy of the material as functions of T
│   └── Mod_Sources.f90         # Assembles the source vector
└── diagnostic/
    └── Mod_Diagnostic.f90      # Residual norms and their output
```

## Where to change what

| To change | Edit |
|---|---|
| An input parameter, its default or its validation | the matching `config/Register_*.f90`, then run `bin/DocGen` |
| A drag or Nusselt correlation | `physics/Lib_Drag.f90` / `Lib_Heat.f90`, and the Python mirror in `test/verification/common.py` |
| An evaporation model | `physics/Lib_Evaporation.f90`, and the Python mirror in `test/verification/F-evaporation/run.py` |
| A closure | the three `physics/Lib_Model_*.f90` plus `ncond` in `config/Assign_Setup.f90` |
| A boundary type | `io/IO_BC.f90` to parse it and `numerics/space/Lib_Ghost.f90` to apply it |
| The time integrator | `numerics/time/explicit/Lib_RK.f90` |

## Module naming

A file `<dir>/<Name>.f90` defines module `ICE_<Name>`, with the `Mod_`/`Lib_` prefix
kept and `_m` suffixes preserved:

| File | Module |
|------|--------|
| `base/Parameters_m.f90` | `ICE_Parameters_m` |
| `parallel/Mod_MPI.f90` | `ICE_Mod_MPI` |
| `physics/Lib_Drag.f90` | `ICE_Lib_Drag` |
| `numerics/fluxes/bc/Mod_BC_Fluxes.f90` | `ICE_Mod_BC_Fluxes` |

Sources are collected by a recursive glob, so a new `.f90` anywhere under `src/lib/` is
picked up on the next configure without editing any `CMakeLists.txt`.

## Selecting behaviour at runtime

Two patterns coexist.

**Integer selectors.** `Lib_Drag`, `Lib_Heat` and `Lib_Evaporation` are collections of
`pure` procedures dispatched by an integer carried in the configuration (`dragSelect`,
`heatSelect`, `evapSelect`, `intfSelect`, `blowSelect`).
Nothing mutable is shared, so the source loops are safe to thread and the functions can
be called from a `pure` context.

**Procedure pointers.** `Lib_Model`, `Lib_Riemann` and `Lib_Limiters` hold module-level
pointers bound once per family, before the loop that uses them
(`assign_all`, `assign_riemann`, `assign_limiter`). Because the binding is per family
and the loops over families are serial, this is safe as written, but it is why
`Lib_Ghost` — which walks every family in one loop — dispatches on the closure name
explicitly instead of using the pointers.

## Build system

CMake 3.23 or newer. The library target is `ICEL` (alias `ICEL::ICEL`); the executables
are `ICE` and `DocGen`, both written to `bin/`. Options are listed under
[Installation](../getting-started/installation.md).

MPI support is a compile-time definition, `USE_MPI`. Everything MPI-specific is behind
it, and `Mod_MPI` provides no-op versions of its interface when it is off, so the rest
of the code has no conditionals in it.
