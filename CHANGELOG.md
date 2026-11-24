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
- `install.sh build --update-deps` selects the newest ORION 2.x again; a [Versioning](docs/development/versioning.md)
  page.
- The documentation is built with Zensical (strict in CI), with a new home page, and its pages were rewritten to
  read more easily; the Overview moved into Getting Started.
- The Code Verification page is split into five topic pages, and the V&V section is grouped by subject.
- `combustion = direct` (Gallier and Godfroy 2009), in `[ICE-Physics]` or as the material's token: drops burn at
  `T-burn` (else `boiling-temperature`) down to `d-stop`, `q-comb` going to the gas; one gas species. Each burning
  material is checked on its own. Its latent heat is the `Lv` column of the property table at the burning
  temperature, else `latent-heat`.
- `combustion = Beckstead`, the d^n law (`K-burn`, `n-burn`, `X-eff`, `beta-part`, `T-ign`), as IGLOO's.
- The constants of a burning model are rows of `INPUT/<prefix>models.txt` (written by ATLAS GPB), one per material
  and kind: `<material> combustion <model> key=value ...`; the row gives `q-comb` (`q-comb=0` allowed).
- A burning material without `d-stop` burns out. Time-accurate, in a run with a gas, the burning drops of a cell at
  or below `d-burnout` (else a drop mass of 1e-20 kg), or whose step burns what is left, leave the cell empty (its
  velocity kept); their mass, momentum and energy are kept as an amount that a coupler takes (`take_removed`) and
  hands to its gas.
  Steady, a burning cell's exit burns a fraction of what leaves it, rising from 0 at the removal size to all of it
  at 0.90 times it (`burnout_band`, printed by the setup). A drop that does not burn is not removed.
- The generated mass of a burning material: a coupler that declares the material its condensed oxide is born as
  (`set_generation`, with the kg generated per kg burnt) gets it born in that material's family, with the burning
  drops' velocity and temperature, each drop at `d-gen` (a constant of the row). The births are the step's rates,
  weighted over the stages as the update weighs them, and the gas does not receive them: time-accurate, the step
  times them is added after every family's step; steady, they are the births in force on the generation family over
  its next step. ICE alone generates no mass.
- The coupling interface, `ICE_Mod_Exchange` (a [Coupling Interface](docs/development/coupling.md) page): the gas
  state a coupler sets (`set_gas`); what a family gives the gas over its last step (`exchange_rate`) and what its
  burn-out removed (`take_removed`), each as the gas's mass, momentum and energy, the generated mass and the heat of
  reaction; where mass is born (`set_generation`, `birth_family`); births from any solver, at once (`add_birth`) or
  as rates in force (`add_birth_rate`, `clear_birth_rates`), each a packet of mass, momentum, energy and number
  (`birth_packet`). ICE's own generated mass is born through the same calls.
- Steady coupled runs: a burning cell steps at most half its bulk density over its burn. The empty cells a source
  fills or a burning family empties step no further than what fills them (time-accurate: the global step at most
  the source's).
- `ICE_gas_by_coupler` (`ICE_Global_m`): a coupler that owns the gas phase sets it, and `gas-phase` is then refused.

### Changed

- **Breaking.** `d-res` is renamed `d-stop`. A models-file row with `d-res` stops naming `d-stop`; in
  `[ICE-Physics]` it is a key ICE does not know, which stops the run.
- **Breaking.** The constants of the burning models leave the material line for the models file: a constant on the
  line stops and names the file, and the line's real keys are `alpha-e`, `k-liq`, `mu-liq`, `T-melt`, `h-fus`,
  `T-nuc` and `cp-solid`. An oxide cap (`cap=` on the line, a cap row, a material of 0 groups) stops: ICE does
  not carry it.
- A breakup model (`breakup=` on the material line, a `breakup` row of the models file) is ignored with a warning:
  ICE has no breakup model (IGLOO has). It stopped the run before, as an unknown key or kind.
- **Breaking.** `q-comb`, `d-stop` and `d-burnout` are no longer keys of `[ICE-Physics]`: there they stop the run,
  naming the models file. A burning material's row must give `q-comb`.
- **Breaking.** `combustion = direct` without `T-burn` stops when `boiling-pressure` is not one atmosphere: the
  burning temperature `boiling-temperature` stands for is the one at 101325 Pa.
- **Breaking.** The burning models are named `COMB_BECKSTEAD` = 1 and `COMB_DIRECT` = 2 (`ICE_Config_Types_m`),
  IGLOO's numbering: a coupler that read `combSelect` = 1 as `direct` reads the names.
- **Breaking.** An evaporating material must give its vapour data (latent heat, vapour molar mass and cp, boiling
  temperature, ...): the run stops naming what is missing, instead of using aluminium's values.
- A negative density after the update, floored silently before, is now counted and warned.
- Near-empty IG/AG cells have their dispersion capped, so an empty domain fills without diverging.
- Two-way coupled with RK2 or RK3: after a step, `source` and `heat` hold the step's rates, the stages' weighted as
  the update weighs them, instead of the last stage's, so that the gas receives what the phase exchanged. Euler, and
  every run that is not two-way coupled, are unchanged.
- **Breaking.** ORION 2.2 or a later 2.x is required, in place of the submodule `lib/ORION`. The configure checks it
  out in `lib/ORION`. Results are unchanged.
- **Breaking.** `ORION_PATH` and `FINER_PATH` are removed: use `FETCHCONTENT_SOURCE_DIR_ORION` and
  `FETCHCONTENT_SOURCE_DIR_FINER` (`install.sh --include-*` as before). `install.sh update --remote` is removed.
- FiNeR (v2.3.0) moves to `lib/FiNeR`; delete an old `lib/third_party`. Its dependencies follow their default
  branches (`ICE_FINER_DEPS=BACKUP` for fixed releases).
- The property table is read by ORION's Tecplot reader: BLOCK packing works too, and ORION's message says why a
  file cannot be read.

### Fixed

- MUSCL face values stay between the two cells of the face on cells of different lengths (as in MOSE): before, a
  longer cell could overshoot its neighbour.
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
