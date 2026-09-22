# Testing

## Running tests

Test cases are located under `test/`. Each case directory is self-contained: `input.ini`,
`INPUT/` (initial and boundary conditions, particle properties), a `reference/` solution and
a `verify.py` that compares the run against it.

```bash
cd test
./test.sh -p 8 check all            # run and verify every case
./test.sh -p 8 check Doisneau/IG    # a single case
./test.sh -p 8 update Doisneau/IG   # rerun and overwrite the stored reference
./test.sh -m 2 -p 4 check all       # 2 MPI ranks x 4 threads (needs an MPI build)
./test.sh clean                     # remove run outputs (run from test/ only)
```

The cases and what they verify are described in [Verification & Validation](../vv/index.md).

## Adding a new test case

1. Create a new directory under `test/`, with `input.ini`, `INPUT/` and `ICE.sh` (copy an
   existing case).
2. Add a `verify.py` that exits with status 0 on success.
3. Add the case to `ALL_TESTS` in `test/test.sh`, then create its reference with
   `./test.sh update <case>`.
4. Document the case and its expected solution on a V&V page.

---
