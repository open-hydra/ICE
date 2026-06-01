# Code Structure

## Repository layout

```
ICE/
├── bin/                    # Compiled executables (generated)
├── build/                  # CMake build directory (generated)
├── cmake/                  # CMake helper files
│   └── ICEConfig.cmake.in
├── docs/                   # MkDocs documentation source
├── lib/                    # Third-party dependencies (submodules)
│   ├── ORION/
│   └── third_party/FiNeR/
├── src/
│   ├── app/                # Entry-point programs
│   │   ├── main.f90        # ICE executable
│   │   └── docgen.f90      # Input registry documentation generator
│   └── lib/                # ICE library source
│       ├── base/           # Fundamental types and parameters
│       ├── config/         # INI file reading and configuration
│       ├── diagnostic/     # Residual monitoring and diagnostics
│       ├── driver/         # Top-level solve/setup/postprocess routines
│       ├── io/             # Solution and restart I/O
│       ├── numerics/       # Flux computation and time integration
│       ├── parallel/       # MPI and ghost-cell exchange
│       └── physics/        # Drag, heat transfer, and phase-change models
├── CMakeLists.txt
├── install.sh
├── LICENSE
├── CITATION.cff
└── mkdocs.yml
```

## Module naming convention

All public modules follow the pattern `ICE_<Subdirectory>_<Name>_m` or `ICE_<Name>_m` for cross-cutting modules. Examples:

| File | Module |
|------|--------|
| `base/Parameters_m.f90` | `ICE_Parameters_m` |
| `base/Advanced_Types_m.f90` | `ICE_Advanced_Types_m` |
| `parallel/Mod_MPI.f90` | `ICE_Mod_MPI` |
| `numerics/fluxes/bc/Mod_BC_Fluxes.f90` | `ICE_Mod_BC_Fluxes` |

## Build system

ICE uses CMake (≥ 3.23). The library target is `ICEL` with alias `ICEL::ICEL`. The main executable target is `ICE`; the documentation generator is `DocGen`.

Key CMake options:

| Option | Default | Description |
|--------|---------|-------------|
| `USE_OPENMP` | `OFF` | Enable OpenMP |
| `USE_MPI` | `OFF` | Enable MPI |
| `USE_TECIO` | `OFF` | Enable TecIO binary output |
| `ORION_PATH` | `lib/ORION/` | Path to ORION submodule |
| `FINER_PATH` | `lib/third_party/FiNeR/` | Path to FiNeR submodule |

---
