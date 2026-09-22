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
└── fast/                      # Short invariant checks, no stored references
    ├── common.sh              # Shared helpers (short run, compare byte for byte)
    ├── openmp-equiv/          # 1 vs 4 threads bit-identical
    └── mpi-equiv/             # 1 vs 2 ranks bit-identical (connection and chimera)
```

Each case under `Doisneau/` is self-contained: `input.ini`, `INPUT/` (initial and
boundary conditions, particle properties), `MESH/`, a stored `reference/` solution, a
`verify.py` that compares the run against it (L2 norm of density, tolerance $10^{-4}$),
and an `ICE.sh` run script. The cases and what they verify are described in
[Verification & Validation](../vv/index.md).

The fast tests own no data: they copy one of the `Doisneau` cases, shorten it to 200
iterations, and run it twice under different parallel settings.

## Running the tests

```bash
ctest --test-dir build -j 5 --output-on-failure   # everything, about a minute
ctest --test-dir build -L fast                    # the fast tier only, ~12 s
ctest --test-dir build -R DoisneauIG -V           # a single case, with live output
```

CTest runs each case in its own source directory and leaves the output there, exactly as
running `./ICE.sh solve` by hand does.

### Labels

| Label | Meaning |
|-------|---------|
| `fast` | Short runs of hard invariants (no stored reference); a few seconds each |
| `validation` | Full case compared against its stored reference solution |
| `needs-mpi` | Needs an MPI build; registered only when `USE_MPI=ON` |
| `MK` / `IG` / `AG`, `chimera`, `connection`, `2D` | Coverage area |

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

## Adding a new test case

1. Create a directory under `test/`, with `input.ini`, `INPUT/` and `ICE.sh` (copy an
   existing case).
2. Add a `verify.py` that exits with status 0 on success.
3. Register it in `test/CMakeLists.txt` with `ice_add_case`, giving it labels.
4. Create its reference solution from a run you trust:
   ```bash
   cd test/<case> && ./ICE.sh solve && cp OUTPUT/part-field.tec reference/
   ```
5. Document the case and its expected solution on a V&V page.

---
