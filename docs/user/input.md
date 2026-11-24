# Input File

`input.ini`, in the case folder, holds every setting of a run. It is an INI file:
sections in square brackets, with `key = value` lines inside them.

- A comment line starts with `;`, `#` or `!`. A `;` also ends a value early, so the
  rest of the line is a comment.
- An unknown key in a section of ICE stops the run and is named: a misspelt key would
  otherwise be ignored without a word. Sections that are not ICE's are left alone, so
  the file can also hold the settings of other tools.
- The section names and most keys are the same as [MOSE](https://github.com/open-hydra/MOSE)'s.

Every key, with its default and allowed values, is listed in the
[parameter reference](registry.md). This page explains how they fit together.

## A minimal file

```ini
[ICE-Parameters]
time-threshold = 0.125

[ICE-Numerics]
time-scheme          = RK2
cfl                  = 0.5
time-accurate        = true
space-reconstruction = MUSCL
flux-limiter         = vanleer

[ICE-Physics]
drag          = Stokes
heat-transfer = Ranz-Marshall

[ICE-IO]
ic-format   = tecplot ascii
sol-format  = tecplot ascii
shell-diter = 100

[ICE-Family1]
closure = AG
```

Five keys are required: `time-scheme`, `cfl`, `time-accurate`,
`space-reconstruction` and `closure`. Everything else has a default.

## Sections

| Section | Contents |
|---|---|
| `[ICE-Parameters]` | Restart, and when to stop |
| `[ICE-Numerics]` | Time scheme, CFL number, reconstruction, limiter, shock detector, Riemann solver, residual smoothing |
| `[ICE-Physics]` | Exchange models, and the properties of the materials and their vapour |
| `[ICE-IO]` | File formats, output frequencies, the folder of the gas field |
| `[ICE-Multigrid]` | Grid levels and the iterations on each |
| `[ICE-Probes]` | The names of the probe sections |
| `[<probe name>]` | One per probe: where, what, how often ([Probes](using.md#probes)) |
| `[ICE-FamilyN]` | One per particle family: its closure |

## Families

ICE counts the `[ICE-FamilyN]` sections from `[ICE-Family1]` up, and stops at the
first missing number. Two families:

```ini
[ICE-Family1]
closure = MK

[ICE-Family2]
closure = AG
```

The families share the mesh, the gas field and the scheme, but each has its own state,
and they do not interact. A file without any family section stops the run.

## Materials

Each family is made of one material. The phase file, `INPUT/part-phase.txt`, written
by ATLAS, lists the materials after its first line, one per line:
`<name> <groups> [key=value ...]`.

```
condensed-dispersed phase
A 1 evaporation=CEM
B 2
```

The families take the materials in order: here family 1 is material A, families 2
and 3 are material B. The groups must add up to the number of families. Without the
file, every family is made of one material. With several materials, the
[property table](initial-conditions.md#property-table) needs one zone per material.

### Properties given per material

These keys of `[ICE-Physics]` take one value per material, in the order of the phase
file (`emissivity = 0 0` for two materials); another count stops the run:

`density`, `specific-heat`, `latent-heat`, `emissivity`, `vapour-molar-mass`,
`boiling-temperature` (or `Tboil`), `boiling-pressure`, `vapour-specific-heat`,
`lewis-number`, `vapour-mass-fraction`, `evaporation-coefficient`.

When the property table has a column for a quantity (`Lv`, `Mv`, `Cpv` or
`Emissivity`), give the key 0 for that material: the table sets it.

### Models chosen on the material line

The tokens are the same as [IGLOO](https://github.com/open-hydra/IGLOO)'s. A token that
also exists in `[ICE-Physics]` overrides it for that material. At set-up, ICE prints the
models of each material.

| Token | Meaning | More |
|---|---|---|
| `evaporation`, `interface`, `alpha-e` | The evaporation model, the interface model and the accommodation coefficient | [Evaporation](#evaporation) |
| `liquid-conduction` | `ITC` (uniform drop temperature) or `P2T` (parabolic profile) | [Conduction and boiling](#conduction-and-boiling) |
| `boiling` | `clamp` or `ZGR`, the behaviour above the boiling temperature | [Conduction and boiling](#conduction-and-boiling) |
| `solidification` | `on`: the material freezes and melts | [Solidification](#solidification) |
| `combustion` | `direct` or `Beckstead`: the drop burns | [Combustion](#combustion) |
| `carry` | A material a burning drop carries (IGLOO's model): ICE does not model it, so the token is refused | |
| `breakup` | The breakup of a drop: ICE has no such model, so the token is ignored with a warning | |

`evaporation-blowing` has no token: it applies to every material.

## Exchange models

The exchange with the gas happens only in a coupled run, when `INPUT/gas.tec` exists.
ICE reads these keys only then. The physics of each model is on the
[Particle Physics](../theory/physics.md) page.

### Drag and heat transfer

`drag` and `heat-transfer` apply to every family, and a coupled run needs both: if
either is left at `none`, the run stops and lists the valid names. `NoDrag` and
`NoHeat` switch the exchange off on purpose. The word `Chang` is refused with a
pointer to `JAXA3`, which is its formula.

### Radiation

Radiation is set by `emissivity`, one value per material, not by a model key.

| Value | Effect |
|---|---|
| 0 (default) | no radiation |
| between 0 and 1 | a grey body of that emissivity |
| above 1 | cut to 1, with a warning |
| negative | the run stops |

An `Emissivity` column in the property table replaces the key for that material.

### Evaporation

`evaporation` is `none` by default; while it is, the vapour keys are not read. An
evaporating material needs six quantities, and the run stops naming the one missing:

| Key | Or the table column |
|---|---|
| `latent-heat` | `Lv`, a function of temperature |
| `vapour-molar-mass` | `Mv` |
| `vapour-specific-heat` (0: the gas cp) | `Cpv`, a function of temperature |
| `boiling-temperature` (or `Tboil`, not both) | |
| `lewis-number` | |
| `vapour-mass-fraction` | |

The saturation pressure follows the Clausius-Clapeyron curve through
`boiling-temperature` at `boiling-pressure` (one atmosphere unless given), or the
`Psat` column of the property table when it has one.

`evaporation-interface` and `evaporation-blowing` refine the chosen model:

```ini
[ICE-Physics]
evaporation           = CEM
evaporation-interface = LK
evaporation-blowing   = LK
latent-heat           = 2.26e6
vapour-molar-mass     = 18.015
vapour-specific-heat  = 0
boiling-temperature   = 373.15
lewis-number          = 1
vapour-mass-fraction  = 0
```

The run stops on `interface = LK` with the `d2-law`.

### Conduction and boiling

Both apply to an evaporating material only: the run stops on `P2T` or `ZGR` for a
material that does not evaporate.

| Token | Values | Needs |
|---|---|---|
| `liquid-conduction` | `ITC` (default), `P2T` | `P2T`: `k-liq` [W/(m K)] and `mu-liq` [Pa s] above 0, or the `Conductivity` and `Viscosity` columns |
| `boiling` | `clamp` (default), `ZGR` | |

### Solidification

All on the material's line of the phase file:

| Token | Meaning | Default |
|---|---|---|
| `solidification` | `on` or `off` | `off` |
| `T-melt` | melting temperature [K] | 2327 |
| `h-fus` | heat of fusion [J/kg] | required, > 0 |
| `T-nuc` | nucleation temperature [K], below `T-melt` | 0.8 `T-melt` |
| `cp-solid` | specific heat of the solid [J/(kg K)] | required, > 0 |

The liquid's specific heat and the density are the material's. The specific heat may
vary with temperature (a `Cp` column); the density must be constant. The run stops on
a `Density` column that varies, and on a material that also evaporates. The set-up
prints each solidifying material's values.

### Combustion

A burning material needs two things:

| What | Where |
|---|---|
| The model, `direct` or `Beckstead` | the `combustion` token on the material's line, or the key in `[ICE-Physics]` for every material whose line names none (the token wins) |
| Its constants | a row of the models file, `INPUT/<prefix>models.txt`, beside the phase file |

ATLAS writes both from `[GPB-Phase*]`. For example:

```
A 1 combustion=direct
```

```
A combustion direct q-comb=9.53e6 d-stop=6e-5 T-burn=2791
```

**The models file.** One row per burning material: `<material> combustion <model> key=value ...`,
the keys in any order. `-` in place of the model stands for the one the line or
`[ICE-Physics]` names. A line starting with `#` is a comment.

**The constants:**

| Key | Meaning | Model | Default |
|---|---|---|---|
| `q-comb` | heat of reaction per unit burnt mass [J/kg] | both | required (0 allowed) |
| `d-stop` | diameter of the inert residue [m] | both | none: the drop burns out |
| `d-burnout` | diameter at which the last drops of a cell are removed [m] | both | a drop mass of 1e-20 kg |
| `d-gen` | diameter of the drops of the generated mass [m], > 0 | both | none; read only when a coupler declares a generated mass |
| `T-burn` | burning temperature [K] | `direct` | `boiling-temperature` |
| `K-burn` | burning coefficient [m$^n$/s] | `Beckstead` | required |
| `n-burn` | exponent of the $d^n$ law, between 0 and 3 | `Beckstead` | 1.8 |
| `X-eff` | oxidiser fraction | `Beckstead` | 1 |
| `T-ign` | ignition temperature [K] | `Beckstead` | 2350 |
| `beta-part` | share of the heat of reaction the drop keeps | `Beckstead` | 0 |

`direct` also needs the latent heat at the burning temperature: the `Lv` column of the
property table, else `latent-heat`. With the column, leave `latent-heat` out or give 0
for that material.

`boiling-temperature` stands for the burning temperature at one atmosphere only: with
`boiling-pressure` at another value, the row must give `T-burn`.

**Residue or burn-out.** With `d-stop` the burning stops at that diameter and the
residue stays. Without it the drops burn out, and the last drops of a cell are removed
([burn-out](../theory/physics.md#burn-out)). A drop that does not burn is never removed.

!!! tip "Comparing with IGLOO"
    To compare a material with its IGLOO run, give `d-burnout` in both decks, at about
    0.1 of the injected diameter: the removal then hands over a thousandth of the drop's
    mass.

**The gas.** ICE running alone takes the burnt vapour to be the gas itself, so the gas
phase file must declare one species: `INPUT/phase.txt`, or `INPUT/<gas-phase>-phase.txt`
with `gas-phase` in `[ICE-Parameters]`. Written by hand:

```
ideal-gas phase
GAS 30.14
```

**What stops the set-up:**

- a burning material that also evaporates or solidifies;
- a constant of the models file left on the material line, or given in `[ICE-Physics]`
  (`q-comb`, `d-stop`, `d-burnout`);
- a row without `q-comb`; `d-stop` beside `d-burnout`; the former key `d-res`;
- a row of another kind (`carry`, IGLOO's carried material), of another model than the
  material's, with a key its model does not read, a second row for one material, or a
  row for no material of the phase file;
- `direct` without its latent heat or its burning temperature;
- `Beckstead` without `K-burn`;
- with ICE running alone, a gas phase file that does not declare one species.

A row of the kind `breakup` is ignored with a warning: ICE has no breakup model.

## Grid levels

```ini
[ICE-Multigrid]
levels      = 2
level1-iter = 5000
level2-iter = 500
```

`levels = 2` adds a coarse grid with half the cells in every direction, so every block
dimension must be even. The run starts on the coarsest level. See
[Grid sequencing](../theory/time-integration.md#grid-sequencing).

Each level needs **its own boundary file**, because the records are per boundary cell:

| Level | File |
|---|---|
| 1 (the fine grid) | `INPUT/part-bc.txt` |
| 2 | `INPUT/part-bc2.txt` |
| *n* | `INPUT/part-bc<n>.txt` |

ATLAS BCB writes them all when `MG-levels` is set in its `[ATLAS-Parameters]`, so the
two settings must agree. ICE checks the number of records of each file against its
level and stops if they differ.

## The time step: `cfl`

`cfl` is the Courant number of the signal speeds **summed over the directions**, on
the width of the cell ([theory](../theory/time-integration.md#the-time-step)). The
largest stable values:

| `time-scheme` | `cfl` up to |
|---|---|
| `RK2` | 1 |
| `RK3` | 1.17 (1.256 with first-order reconstruction) |
| `euler` | 1 with first-order reconstruction; with `MUSCL`, the limiter's bound (below) |

With `euler` and `MUSCL` the bound is 0.5 for `vanleer`, `mc`, `superbee`, `koren`,
`umist` and `ospre`; 0.586 for `vanalbada`; 0.571 for `sweby` and `osher`; 0.667 for
`minmod`. Above it ICE warns.

Keep a margin where shocks and limiters act, for example `cfl = 0.8` with `RK2`.

## Stopping conditions

The run stops at the first of three thresholds:

| Key | Stops when | Default |
|---|---|---|
| `iter-threshold` | the iterations reach it | never, in practice |
| `time-threshold` | the simulation time reaches it | never, in practice |
| `res-threshold` | the density residual has fallen by this factor (steady runs only) | $10^{-10}$ |

So a steady run with no threshold set stops on the residual. `res-threshold = 0`
switches the residual test off. A time-accurate run ignores it: the change of one time
step says nothing about convergence.

`res-threshold` is **relative**: it compares the density residual (column 3 of the
[residual history](using.md#the-residual-history)) with the first non-zero one, on the
finest grid. The same flow therefore stops at the same iteration whatever its
loading. The first *non-zero* residual is the reference because a cloud at rest that
the gas starts moving may not change at its first step.

- A steady restart takes back the reference of the run it continues
  ([Restarting](using.md#restarting)).
- A `res-threshold` switched on while the run goes on (from 0, see below) measures from
  the first non-zero residual after it.
- A residual that first rises and then falls is compared with the first one, not with
  its peak.

## Reading again during the run

With `ini-diter` set, ICE reads `input.ini` again every that many iterations:

- the thresholds and the output frequencies take effect at once. A threshold lowered
  below the iterations already done ends the run (or the grid level) at the next
  iteration;
- what is set up once does not change: the scheme, and the materials with their
  properties, models and tables. `irs = true` in a time-accurate run is ignored, with a
  warning.

An error in the new file — a value that breaks its rule, a changed `levels` or
`closure`, `riemann-solver = Saurel` with an IG or AG family, a file that cannot be
read or holds no key — is handled as `ini-error` says:

| `ini-error` | What happens |
|---|---|
| `warn` (default) | ICE warns, with the iteration and the error, and keeps **all** its previous settings, including the keys the new file had right. Keys that go together, such as `cfl` and `dt-max`, never change one without the other. |
| `stop` | The run stops with the error and exit code 1. |

The value of `ini-error` in force when the reading starts decides. At setup, every
error stops the run.
