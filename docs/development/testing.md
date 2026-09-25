# Testing

ICE's test cases are wired into **CTest**, so the whole suite runs from the build
directory and can guard a `git push`.

## Test organisation

```
test/
├── CMakeLists.txt             # CTest definitions (cases, labels, timeouts)
├── Doisneau/                  # Crossing-jets cases, one per closure and topology
│   ├── MK/  IG/  AG/          # Single block, 100 x 100
│   ├── IG-chimera/            # Two overlapping non-matching blocks (ATLAS 102)
│   ├── IG-split4/             # Four blocks joined by connections (ATLAS 101)
│   └── plot_vv.py             # Regenerates the V&V figures from the case outputs
├── NoExchange/MK/             # Coupled slab with drag = NoDrag, heat-transfer = NoHeat: state kept bit for bit
├── Axis/                      # One-degree wedge about x: side faces 200, axis face 300
│   ├── MK/                    # Pressureless cloud along x: stationary field, boundary census
│   └── IG/  AG/               # Uniform cloud at rest in a closed wedge: stays at rest (hoop pressure)
├── Thermal/                   # Expanding blob in a closed box at uniform T
│   └── IG/  AG/               # T stays uniform: the pressure work in the energy flux
├── Riemann/AG/                # Anisotropic Sod tube against the exact gamma = 3 solution
├── Reflect/AG/                # Sheared blob on a code-200 symmetry plane: x-momentum conserved
├── Berthon/                   # 1D Riemann problems for the AG closure
│   ├── SCS/  RCS/  RCR/       # Shock-contact-shock, rarefaction-contact-shock, ...
│   ├── SCS-hlle/              # SCS again through the HLLE flux
│   ├── Results/               # The analytical wave patterns each case is checked against
│   ├── berthon.py             # Exact-solution reader and the L1 comparison
│   └── plot_vv.py             # Regenerates the V&V figures from the case outputs
├── verification/              # Exact-solution cases, inputs generated on the fly
│   ├── common.py              # Case generation, reference correlations and ODE solver
│   ├── plot_vv.py             # Regenerates the code-verification figures
│   ├── A-drag-relaxation/     # Stokes drag against its closed form + RK order study
│   ├── B-thermal-relaxation/  # Nu = 2 heat exchange against its closed form
│   ├── C-advection-sine/      # Periodic transport + grid-refinement order study
│   ├── D-drag-matrix/         # All 12 drag laws against an independent integration
│   ├── E-heat-matrix/         # All 7 Nusselt laws, likewise
│   ├── F-evaporation/         # All 5 evaporation models against the d-squared law
│   ├── G-cloud-translation/   # Transport + drag together, against exact translation
│   ├── H-linear-strain/       # Straining gas, exact affine map + small-St asymptote
│   ├── I-vortex-cloud/        # 2D cloud in a prescribed vortex, exact conformal map
│   └── K-properties-table/    # Property table read at a fixed T: interpolation, range, saturation
├── fast/                      # Short invariant checks, no stored references
│   ├── common.sh              # Shared helpers (short run, compare byte for byte)
│   ├── openmp-equiv/          # 1 vs 4 threads bit-identical
│   ├── refusals/              # Broken inputs and a diverging run must exit non-zero
│   └── mpi-equiv/             # 1 vs 2 ranks bit-identical (connection and chimera)
└── unit/                      # Programs linked against the library, no solver run
    └── test_properties.f90    # Property table: grammar, checks, loading, lookup, energy, Psat
```

Each case under `Doisneau/` is self-contained: `input.ini`, `INPUT/` (initial and
boundary conditions, particle properties), `MESH/`, a stored `reference/` solution, a
`verify.py` that compares the run against it (L2 norm of the density, tolerance $10^{-4}$, and
of the number density relative to its own scale, held to the same tolerance over the density's RMS),
and an `ICE.sh` run script. The cases and what they verify are described in
[Verification & Validation](../vv/index.md).

Each case under `Berthon/` is laid out the same way, but has no stored reference of its
own: its `verify.py` compares six fields against the analytical wave pattern in
`Berthon/Results/` (L1 norm over the domain), and `IBCB.f90` is the small program that
wrote its initial and boundary conditions. See
[Berthon Riemann Problems](../vv/berthon.md).

The cases under `Axis/` and `NoExchange/` carry no stored reference: each keeps the
generator of its inputs (`make_case.py`) and a `verify.py` that checks an exact property
of the run. The `Axis/` run scripts keep the solver's log in `logfile`, which their
`verify.py` reads.

The fast tests own no data: they copy one of the `Doisneau` cases, shorten it to 200
iterations, and run it twice under different parallel settings.

The unit tests call the library's routines directly and run no case. Each writes the
small files it reads into its own directory in the build tree, `build/test/unit/<name>/`,
and prints one line per assertion.

The verification cases own no data either: each `run.py` writes its own mesh, initial
condition, boundary conditions and `input.ini` into a scratch `work/` directory, runs
the solver there and compares against a closed form or an independently integrated
reference. What each one verifies is described in
[Code Verification](../vv/verification.md). Set `ICE_KEEP_WORK=1` to keep the generated
cases for inspection after a pass.

## Running the tests

```bash
ctest --test-dir build -j 5 --output-on-failure   # everything, about a minute
ctest --test-dir build -L fast                    # the fast tier only, ~15 s
ctest --test-dir build -R DoisneauIG -V           # a single case, with live output
```

CTest runs each case in its own source directory and leaves the output there, exactly as
running `./ICE.sh solve` by hand does.

### Labels

| Label | Meaning |
|-------|---------|
| `fast` | Short runs of hard invariants (no stored reference); a few seconds each |
| `validation` | Full case compared against its stored reference solution |
| `verification` | Compared against a solution that does not come from ICE (closed form or independent integration) |
| `unit` | A program linked against the library that asserts what its routines return |
| `needs-mpi` | Needs an MPI build; registered only when `USE_MPI=ON` |
| `sources`, `transport`, `implementation` | What part of the solver a case covers |
| `MK`, `IG`, `AG`, `chimera`, `connection`, `1D`, `2D` | What configuration it covers |

Labels combine, so `ctest -L AG` runs every case touching the anisotropic Gaussian
closure whatever its tier, and `ctest -L fast` is what the pre-push hook runs.

### Threads and ranks

The cases run on one thread and one rank by default. To change that, configure with:

```bash
cmake -B build -DICE_TEST_THREADS=4 -DICE_TEST_RANKS=2
```

`ICE_TEST_RANKS` above 1 needs an MPI build. Both builds write `bin/ICE`, so the binary
is whichever one was built last; the MPI tests check that they really got the MPI build
and fail with a clear message if not.

## The pre-push hook

`.githooks/pre-push` runs the suite before every push. Enable it once per clone:

```bash
git config core.hooksPath .githooks
```

It skips itself if there is no configured build directory, and otherwise picks `build/`,
or the first `build-*/`. Useful overrides:

```bash
ICE_PREPUSH_LABEL=fast git push   # only the fast tier
ICE_BUILD_DIR=/path/to/build git push
ICE_PREPUSH_JOBS=2 git push
```

## Adding a case

Decide first which kind it is, because the kinds are registered differently and carry
different weight.

### A validation case

It compares a full run against a stored solution of ICE's own, so it can only detect
*change*. Use one when the physics has no closed form but a regression would matter.

1. Create a directory under `test/`, with `input.ini`, `INPUT/` and `ICE.sh` — copy an
   existing case rather than writing `ICE.sh` from scratch.
2. Write a `verify.py` that exits 0 on success and prints enough to diagnose a failure.
3. Produce the reference from a run you have reason to trust, and say in the commit why
   you trust it:
   ```bash
   cd test/<case> && ./ICE.sh solve && mkdir -p reference && cp OUTPUT/part-field.tec reference/
   ```
4. Register it with `ice_add_case`, labelled `validation` plus its coverage.

### A verification case

It compares against a solution that does not come from ICE — a closed form, or an
integration done independently — so it can fail on its first run, and a failure means
something. Prefer one whenever the problem admits an exact answer.

1. Add a directory under `test/verification/` with a `run.py` that writes its own mesh,
   initial condition, boundary table and `input.ini` into `work/`, runs the solver there
   and compares. `common.py` has the helpers for all of that.
2. Register it with `ice_add_verification`, labelled `verification` plus its coverage,
   and add `fast` if it runs in a few seconds.

### A unit test

It asserts a routine's contract directly, without a run: a parser, a check, a lookup.
Use one when the contract is a routine's and a case would reach it only indirectly.

1. Add `test/unit/test_<name>.f90`, a program that `use`s the library modules, prints
   `OK` or `FAIL` with a description for every assertion, and ends with `error stop 1`
   if any failed.
2. Register it with `ice_add_unit_test(<Name> unit/test_<name>.f90 "fast;unit;implementation")`.

For a validation or verification case, document what it checks and what tolerance it
uses on a [V&V page](../vv/index.md) — a tolerance with no recorded reasoning is one
nobody can tighten later.

---
