# Contributing

## Workflow

1. Fork the repository and create a feature branch.
2. Make your changes following the coding conventions below.
3. Test your changes against the existing test cases.
4. Open a pull request with a clear description of what was changed and why.

## Coding conventions

- **Language**: Modern Fortran (F90/F95 modules, `implicit none` everywhere).
- **Module names**: `ICE_<Description>_m` for library modules.
- **File names**: match the module name without the `ICE_` prefix (e.g., `Advanced_Types_m.f90`).
- **Comments**: add comments only where the *why* is non-obvious. No docstrings for trivial procedures.
- **Indentation**: 2 spaces.

## Submitting issues

Use the GitHub issue tracker. Include a minimal reproducible case whenever possible.

---
