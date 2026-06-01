# Installation

This document describes how to obtain and build **ICE**. The instructions cover the `install.sh` script, CMake configuration, and Git submodule layout.

!!! note
    ICE has a dual nature: it is both a library and an executable. The installation process produces both the static library `libICEL.a` and the main executable `bin/ICE`. If you are only interested in using ICE as a library, you can link against `libICEL.a` without caring about the executable.

## Prerequisites

Before attempting to build ICE make sure your system provides the following external tools and compilers:

- **CMake** – 3.23 or newer.
- **Fortran compiler** – either the GNU toolchain (`gfortran`) or Intel/oneAPI (`ifort`/`ifx`) are supported.
- **C / C++ compiler** – required by the ORION I/O library and the optional TecIO component.
- **OpenMP** – needed for optional shared-memory parallelisation.
- **MPI** – needed for optional distributed-memory parallelisation.

### Git submodules

ICE depends on two repositories included as Git submodules.

| Path | Repository URL | Purpose |
|------|----------------|---------|
| `lib/ORION` | `https://github.com/MarcoGrossi92/ORION.git` | I/O routines (TecIO, VTK) |
| `lib/third_party/FiNeR` | `https://github.com/szaghi/FiNeR.git` | INI file parser |

## Build methods

First clone the repository with submodules:

```bash
git clone https://github.com/open-hydra/ICE.git
cd ICE
git submodule update --init --recursive
```

### Build with `install.sh` (recommended)

The script exposes three commands: `build`, `compile`, and `update`. It also maintains a `CMakePresets.json` file that records the configuration used for the most recent `build` invocation.

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**`build` command**

```bash
# minimal GNU build with OpenMP enabled
./install.sh build --compilers=gnu --use-openmp

# full configuration with MPI and all optional features
./install.sh build --compilers=gnu --use-openmp --use-mpi --use-tecio
```

Options accepted by `build`:

* `--compilers=<gnu|intel>` – select the compiler family (default: `gnu`).
* `--use-openmp` – enable OpenMP parallelisation.
* `--use-mpi` – enable MPI parallelisation.
* `--use-tecio` – enable TecIO support.
* `--include-orion=PATH` – use an external ORION tree instead of the submodule.
* `--include-finer=PATH` – same for FiNeR.

**`compile` command**

Re-runs CMake using the previously generated preset. Useful during development when only source files have changed.

```bash
./install.sh compile
```

**`update` command**

Synchronises the Git submodules.

```bash
./install.sh update           # sync to recorded commit
./install.sh update --remote  # update to newest remote commit
```

### Build with CMake

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DUSE_OPENMP=ON \
    -DUSE_MPI=OFF \
    -DUSE_TECIO=ON
cmake --build . --parallel
```

## CMake presets

After a successful `build`, a `CMakePresets.json` is written in the source root. Subsequent builds can reuse it:

```bash
cmake --preset default
cmake --build build
```

## Library linking (advanced)

To use ICE from an external CMake project:

```cmake
find_package(ICE REQUIRED)
add_executable(myapp main.f90)
target_link_libraries(myapp ICE::ICEL)
```

## Next steps

* **[Quick Start](quick-start.md)** – build and verify the installation.

---
