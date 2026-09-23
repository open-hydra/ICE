# Contributing

## Workflow

1. Branch from the default branch.
2. Make the change, and add or extend a test case for it under `test/`.
3. Run the suite: `ctest --test-dir build -j 5 --output-on-failure`.
4. Open a pull request describing what changed and why.

Installing the pre-push hook runs the suite before every push:

```bash
git config core.hooksPath .githooks
```

See [Testing](testing.md) for the tiers and how to add a case.

## Coding conventions

- **Language**: modern Fortran, free form, `implicit none` everywhere, modules rather
  than common blocks.
- **Modules**: `ICE_<Name>`, matching the file name; see
  [Code Structure](structure.md#module-naming).
- **Indentation**: two spaces. No tabs.
- **Precision**: `real(R8)` throughout, from `iso_fortran_env`; write literals as
  `1._R8`, never as bare `1.0`.
- **Comments**: explain *why*, not *what*. A comment that restates the line below it is
  noise; one that records why a loop is split, or why a guard exists, is not.
- **Purity**: new leaf procedures — correlations, limiters, flux functions — should be
  `pure` and take everything they need as arguments, rather than reaching for module
  state.

## Changing an input parameter

The registry in `src/lib/config/Register_*.f90` is the single source of truth for
inputs: it holds the default, the validation rule and the description. After changing a
`reg%add` call, regenerate the reference page:

```bash
./bin/DocGen
```

and commit `docs/user/registry.md` with the change.

## Changing a physical model

The verification suite mirrors ICE's drag and Nusselt correlations in Python, in
`test/verification/common.py`, so that the comparison is genuinely independent of the
Fortran. A change to one has to be made in both, or case D or E will fail — which is
the point.

## Documentation

The docs are MkDocs with the Material theme, built from `docs/` and `mkdocs.yml`:

```bash
mkdocs serve          # preview at localhost:8000
mkdocs build --strict # what CI checks: fails on a broken internal link
```

The figures on the V&V pages are regenerated from the case outputs by the `plot_vv.py`
scripts under `test/`, after the corresponding cases have run.

## Submitting issues

Use the GitHub issue tracker, with a minimal case that reproduces the problem — the
`input.ini` and the smallest mesh that still shows it.
