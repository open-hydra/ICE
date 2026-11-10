# Installation

ICE builds a static library, `libICEL.a`, and two executables: `bin/ICE`, the solver,
and `bin/DocGen`, which regenerates the [parameter reference](../user/registry.md)
from the source.

## Prerequisites

| | |
|---|---|
| **CMake** | 3.23 or newer |
| **Fortran compiler** | GNU `gfortran` or Intel oneAPI `ifx` |
| **C / C++ compiler** | For ORION, and for TecIO if enabled |
| **MPI** | Optional, for distributed-memory runs |
| **OpenMP** | Optional, for shared-memory runs |

Both parallel modes are off by default and can be combined. The OpenMP build needs
Linux: on macOS, build without it.

## Get the code

```bash
git clone --recurse-submodules https://github.com/open-hydra/ICE.git
cd ICE
```

ICE uses two libraries:

| Library | Purpose | Where it comes from |
|---------|---------|---------------------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Mesh and solution I/O (Tecplot, VTK, optional TecIO) | The newest 2.x release from 2.1.0, cloned into `lib/ORION` by the first configure |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI file parser | The submodule `lib/FiNeR` (v2.3.0) |

The first configure needs network access: besides ORION, it clones FiNeR's five
dependencies (PENF, FACE, FLAP, BeFoR64, StringiFor) into
`lib/FiNeR/src/third_party/`. The ORION release it picks is recorded in `deps.lock`,
so later configures work offline. `lib/ORION` and `deps.lock` are not tracked by git.
[Versioning](../development/versioning.md) explains how ICE chooses the releases.

## Build with `install.sh`

```bash
# GNU compilers, OpenMP
./install.sh build --compilers=gnu --use-openmp

# Intel compilers, MPI and binary Tecplot files
./install.sh build --compilers=intel --use-mpi --use-tecio
```

`build` starts from scratch (it removes `build/`), compiles, and records the settings
in `CMakePresets.json`.

| Option | Effect |
|---|---|
| `--compilers=gnu` | `gfortran`, `gcc`, `g++` |
| `--compilers=intel` | `ifx`, `icx`, `icpx` |
| `--use-openmp` | OpenMP threads |
| `--use-mpi` | MPI ranks; CMake finds the compiler wrappers |
| `--use-tecio` | Binary Tecplot files (`.szplt`); TecIO is built on the first configure, which takes a few minutes |
| `--include-orion=PATH` | Build this ORION clone instead of the release |
| `--include-finer=PATH` | Build this FiNeR clone instead of the submodule |
| `--update-deps` | Pick the newest ORION release again |

Two more commands:

```bash
./install.sh compile   # rebuild with the recorded settings, after editing the sources
./install.sh update    # put lib/FiNeR back at the commit ICE records
```

During development, use `compile`: it skips the configure step.

## Build with CMake

```bash
cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DUSE_OPENMP=ON
cmake --build build --parallel
```

| Variable | Default | Meaning |
|---|---|---|
| `USE_OPENMP` | `OFF` | OpenMP threads |
| `USE_MPI` | `OFF` | MPI ranks |
| `USE_TECIO` | `OFF` | Binary Tecplot files; without it, `tecplot binary` is refused when a run starts |
| `ICE_BIN_DIR` | `bin/` | Where the executables go: give each build folder its own, or they overwrite each other's `bin/ICE` |
| `FETCHCONTENT_SOURCE_DIR_ORION` | unset | An ORION clone to build instead of the release |
| `FETCHCONTENT_SOURCE_DIR_FINER` | unset | A FiNeR clone to build instead of the submodule |
| `SEMVER_UPDATE` | `OFF` | `ALL` picks the newest ORION release again, once |
| `ICE_PIN_CHECK` | `FATAL_ERROR` | `WARNING` accepts a library release that ICE does not require |
| `ICE_FINER_DEPS` | `FREE` | `BACKUP` takes FiNeR's dependencies at fixed releases instead of their default branches |

The compiler flags come from `cmake/SetFortranFlags.cmake`, which checks that the
compiler accepts each flag before using it.

### External clones

A clone given with `FETCHCONTENT_SOURCE_DIR_ORION` or `FETCHCONTENT_SOURCE_DIR_FINER`
is built as it is, uncommitted changes included. It must be a release ICE accepts,
or the configure stops: ORION 2.1.0 or a later 2.x, and FiNeR v2.3.0. A FiNeR clone
must already hold its dependencies in `src/third_party` (`fobis fetch --no-build`).
To build another release on purpose, such as an ORION branch under test, add
`-DICE_PIN_CHECK=WARNING`.

## Run the tests

When ICE is the top-level project, every build also registers the test suite:

```bash
ctest --test-dir build -j 5 --output-on-failure
```

The MPI cases are registered only with `USE_MPI=ON`. See [Testing](../development/testing.md).

## Use ICE as a library

The library target is `ICEL`, also called `ICEL::ICEL`. Add the source tree to your
project:

```cmake
add_subdirectory(path/to/ICE)
target_link_libraries(myapp PRIVATE ICEL::ICEL)
```

Inside another project, ICE builds the library only: no executables and no tests.
The parent project provides the `ORION` and `FiNeR::FiNeR` targets, as Hydra does.
If the parent selects ORION with `cmake/SemVerDeps.cmake`, ICE checks that the release
satisfies its own requirement.

## Troubleshooting

**`undefined reference` to C++ symbols when linking with TecIO.** TecIO is a C++
library linked into Fortran executables. A standalone build with `USE_TECIO=ON`
handles this by itself; a project that embeds ICE must enable C++ itself.

**`lib/FiNeR is checked out at ..., but this ICE commit records ...`.** The submodule
moved, usually after a `git checkout` of another branch. Run `./install.sh update`
and configure again.

**`lib/ORION ... is on a branch or has changes`.** The configure moves `lib/ORION` to
another release only if you have not touched it. Commit or stash your changes, or
build the folder as it is with `-DFETCHCONTENT_SOURCE_DIR_ORION=lib/ORION`.

**A clone from before ORION became a release.** Delete the folder `lib/third_party`
(FiNeR is now in `lib/FiNeR`) and run `./install.sh build` again.

## Next steps

* **[Quick Start](quick-start.md)** — run a shipped case and check its result.
