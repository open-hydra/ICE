# Changelog

All notable changes to ICE are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). An entry that begins with **Breaking.** asks a user to
change an input or a build.

## [Unreleased]

### Added

- Drop internal heating and boiling, as in IGLOO: `liquid-conduction = P2T` (parabolic temperature profile) and
  `boiling = ZGR` (Zuo-Gomes-Rutland), with `k-liq` and `mu-liq`.
- `boiling-pressure`: the pressure at which `boiling-temperature` is given (default one atmosphere).
- The property table takes an `Emissivity` column by material, and a liquid cp that varies with temperature.
- The banner prints the compiler options of the build.
- CI builds and tests with TecIO (`USE_TECIO=ON`) too.

### Changed

- **Breaking.** An evaporating material must give its vapour data (latent heat, vapour molar mass and cp, boiling
  temperature, ...): the run stops naming what is missing, instead of using aluminium's values.
- A negative density after the update, floored silently before, is now counted and warned.
- Near-empty IG/AG cells have their dispersion capped, so an empty domain fills without diverging.

### Fixed

- Probes accept the pressure names of the solution files, and as many names as the key holds.
- A build without TecIO refuses `tecplot binary` at set-up, not at the first read or write.
- Warnings: forward Euler with MUSCL above its limiter's CFL bound, and a ZGR step that removes over half a cell's
  drop mass.

## [0.1.0-alpha] - 2026-11-04

First public release. See the [README](README.md) for the features.

### Added

- MK, IG and AG closures; drag, convective heat, radiation, evaporation and solidification.
- Multi-block grids with connections and chimera overset; MPI over blocks, OpenMP with thread groups.
- Materials and inlets from ATLAS's phase files, including the 401 mixture inlet over several families.
- Restart that resumes the counters, and residual history per family.
- Tecplot (ASCII and binary) and VTK output; point probes.

[Unreleased]: https://github.com/open-hydra/ICE/compare/v0.1.0-alpha...HEAD
[0.1.0-alpha]: https://github.com/open-hydra/ICE/releases/tag/v0.1.0-alpha
