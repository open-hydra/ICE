# Installation

ICE builds into a static library `libICEL.a` and two executables: `bin/ICE`, the
solver, and `bin/DocGen`, which regenerates the
[parameter reference](../user/registry.md) from the input registry in the source.

## Prerequisites

| | |
|---|---|
| **CMake** | 3.23 or newer |
| **Fortran compiler** | GNU `gfortran` or Intel oneAPI `ifx` |
| **C / C++ compiler** | Required by ORION, and by TecIO if enabled |
| **MPI** | Optional, for distributed-memory runs |
| **OpenMP** | Optional, for shared-memory runs |

Both parallel modes are off by default and can be combined.

### Submodules

ICE bundles two dependencies as Git submodules. `install.sh` initialises them, and so
does CMake if they are missing.

| Path | Repository | Purpose |
|------|------------|---------|
| `lib/ORION` | `github.com/MarcoGrossi92/ORION` | Mesh and solution I/O (Tecplot, VTK, optional TecIO) |
| `lib/third_party/FiNeR` | `github.com/szaghi/FiNeR` | INI file parser |

```bash
git clone https://github.com/open-hydra/ICE.git
cd ICE
git submodule update --init
```

FiNeR's own five dependencies (PENF, FACE, FLAP, BeFoR64, StringiFor) are not submodules, so
a recursive clone does not fetch them. The CMake configure clones them from
`github.com/szaghi/` into `lib/third_party/FiNeR/src/third_party/`, which FiNeR's
`.gitignore` excludes, so the submodule stays clean. A bare `cmake -B build` does the same, and it needs network access
the first time. The same applies to a FiNeR tree given with `--include-finer`.

## Building with `install.sh`

```bash
# GNU build with OpenMP
./install.sh build --compilers=gnu --use-openmp

# Intel build with MPI and binary Tecplot support
./install.sh build --compilers=intel --use-mpi --use-tecio
```

`build` configures from scratch — it removes `build/` first — compiles, and writes a
`CMakePresets.json` recording what it used.

| Option | Effect |
|---|---|
| `--compilers=gnu` | `gfortran` / `gcc` / `g++` |
| `--compilers=intel` | `ifx` / `icx` / `icpx` |
| `--use-openmp` | Thread parallelism |
| `--use-mpi` | Rank parallelism; the compiler wrappers are found by CMake |
| `--use-tecio` | Binary Tecplot support in ORION |
| `--include-orion=PATH` | Use an existing ORION tree instead of the submodule |
| `--include-finer=PATH` | The same for FiNeR |

Two further commands:

```bash
./install.sh compile          # rebuild from the recorded preset, after editing sources
./install.sh update           # sync the submodules to their recorded commits
./install.sh update --remote  # ... or to the newest remote commit
```

`compile` is the one to use during development: it skips the configure step.

## Building with CMake directly

```bash
cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DUSE_OPENMP=ON \
    -DUSE_MPI=OFF \
    -DUSE_TECIO=OFF
cmake --build build --parallel
```

| Variable | Default | Meaning |
|---|---|---|
| `USE_OPENMP` | `OFF` | OpenMP |
| `USE_MPI` | `OFF` | MPI; also defines `USE_MPI` for the preprocessor |
| `USE_TECIO` | `OFF` | Passed to ORION for binary Tecplot |
| `ORION_PATH` | `lib/ORION/` | ORION source tree |
| `FINER_PATH` | `lib/third_party/FiNeR/` | FiNeR source tree |

Both executables land in `bin/` whatever the build directory, so two build directories
— one serial, one MPI — will overwrite each other's `bin/ICE`. Keep them apart if you
need both at once.

Compiler flags come from `cmake/SetFortranFlags.cmake`, which probes the compiler for
each flag rather than assuming it.

## Building the test suite

`enable_testing()` and the cases are added automatically when ICE is the top-level
project, so any of the builds above gives a working `ctest`:

```bash
ctest --test-dir build -j 5 --output-on-failure
```

The MPI-specific cases register only when `USE_MPI=ON`. See
[Testing](../development/testing.md).

## Using ICE as a library

The library target is `ICEL`, aliased `ICEL::ICEL`. ICE does not install an exported
CMake package, so an external project consumes it by adding the source tree:

```cmake
add_subdirectory(path/to/ICE)
target_link_libraries(myapp PRIVATE ICEL::ICEL)
```

Everything ICE's top-level `CMakeLists.txt` does is guarded on being the top-level
project, so as a subdirectory it builds only the library — no solver executable and no
tests — and the parent project must already define the `ORION` and `FiNeR::FiNeR`
targets it links against. This is how the Hydra suite embeds it.

## Troubleshooting

**`undefined reference` to C++ symbols when linking with TecIO.** TecIO is C++ and the
executable is linked by the Fortran driver, which does not pull in the C++ runtime.
Add it explicitly:

```bash
cmake -B build ... -DCMAKE_EXE_LINKER_FLAGS=-lstdc++
```

**`cannot find -ltecio::tecio`.** The ORION submodule is older than the ICE revision
that expects it. `./install.sh update --remote` and rebuild.

## Next steps

* **[Quick Start](quick-start.md)** — run a shipped case and check it.
