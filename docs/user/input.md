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

### Keys with one value per material

These keys of `[ICE-Physics]` take one value per material, in the order of the phase
file (`emissivity = 0 0` for two materials); another count stops the run:

`density`, `specific-heat`, `latent-heat`, `emissivity`, `vapour-molar-mass`,
`boiling-temperature` (or `Tboil`), `boiling-pressure`, `vapour-specific-heat`,
`lewis-number`, `vapour-mass-fraction`, `evaporation-coefficient`.

When the property table has a column for a quantity (`Lv`, `Mv`, `Cpv` or
`Emissivity`), give the key 0 for that material: the table sets it.

### Tokens on the material line

The tokens are the same as [IGLOO](https://github.com/open-hydra/IGLOO)'s. The ones
that also exist in `[ICE-Physics]` override it for that material.

| Token | Meaning |
|---|---|
| `evaporation`, `interface`, `alpha-e` | The evaporation model, the interface model and the accommodation coefficient |
| `liquid-conduction` | `ITC` (uniform drop temperature) or `P2T` (parabolic profile, with `k-liq` [W/(m K)] and `mu-liq` [Pa s], or the `Conductivity` and `Viscosity` columns) |
| `boiling` | `clamp` or `ZGR`, the behaviour above the boiling temperature |
| `solidification = on` | The material freezes and melts, with `T-melt`, `h-fus`, `T-nuc` and `cp-solid` |
| `combustion = direct` | The drop burns on the d^2 law at its burning temperature; its constants are a row of the models file ([Combustion](#combustion)) |
| `combustion = Beckstead` | The drop burns on the d^n law from its ignition temperature; its constants are a row of the models file |
| `K-burn`, `n-burn`, `X-eff`, `beta-part`, `T-ign`, `T-burn`, `q-comb`, `d-stop`, `d-burnout`, `d-gen` | The constants of the burning models: refused on the line, their place is the models file |
| `cap` | The oxide cap of a burning drop, which ICE does not carry: refused, and so is a `cap` row of the models file |
| `breakup` | The breakup model of a drop, which ICE does not have: ignored with a warning, and so is a `breakup` row of the models file |

The physics of each is in [Particle Physics](../theory/physics.md):
[conduction inside the drop](../theory/physics.md#conduction-inside-the-drop),
[boiling](../theory/physics.md#boiling-above-the-boiling-temperature),
[solidification](../theory/physics.md#solidification) and
[direct combustion](../theory/physics.md#direct-combustion). At setup, ICE prints the models of each material.

The setup stops on:

- an unknown key, a value that is not a number, or `interface = LK` with the d2-law;
- solidification without a positive `h-fus` and `cp-solid`, with `T-nuc` at or above
  `T-melt`, together with evaporation, or with a `Density` column that varies;
- in a coupled run: `P2T` or `ZGR` on a material that does not evaporate, and `P2T`
  without a positive `k-liq` and `mu-liq` (or their columns);
- a burning material that also evaporates or solidifies, that gives `d-stop` beside `d-burnout`,
  the former key `d-res`, or a key of the oxide cap, or whose row of the models file has no
  `q-comb`; `combustion = direct` without its latent heat and its burning temperature, or without
  `T-burn` when `boiling-pressure` is not one atmosphere; `combustion = Beckstead` without `K-burn`;
  `q-comb`, `d-stop` or `d-burnout` in `[ICE-Physics]`; with ICE running alone, a gas phase file
  that does not declare one species.

`evaporation-blowing` has no token: it applies to every material.

## Exchange models

`drag` and `heat-transfer` apply to every family; `evaporation` is the default of every
material. ICE reads them only in a coupled run, when `INPUT/gas.tec` exists.

- A coupled run needs `drag` and `heat-transfer`: if either is left at `none`, the run
  stops and lists the valid names. `NoDrag` and `NoHeat` switch the exchange off on
  purpose.
- Radiation is set by `emissivity`, not by these keys.
- `evaporation` is `none` by default, and while it is, the vapour keys are not read.

### Evaporation

An evaporating material needs six quantities. The run stops and names the one missing.

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

### Combustion

`combustion = direct` burns the drops of a material at its burning temperature, `combustion =
Beckstead` on the d^n law from `T-ign` ([Particle Physics](../theory/physics.md#direct-combustion)).
Each material burns on its own model: the word is a token of the material's line of the phase
file, and the constants of the model are a row of the models file, `<prefix>models.txt` beside the
phase file. ATLAS writes both from `[GPB-Phase*]`:

```
A 1 combustion=direct
```

```
A combustion direct q-comb=9.53e6 d-stop=6e-5 T-burn=2791
```

A row is `<material> combustion <model> key=value ...`, the keys in any order; `-` in place of the
model stands for the one the line or `[ICE-Physics]` names; a line that starts with `#` is a comment.
`direct` reads `T-burn` [K], `q-comb` [J/kg], `d-stop` [m], `d-burnout` [m] and `d-gen` [m];
`Beckstead` reads `K-burn`, `n-burn`, `X-eff`, `beta-part`, `T-ign` [K], `q-comb`, `d-stop`,
`d-burnout` and `d-gen`. The setup
stops on one of these constants left on the material line, on a row of another kind (`cap`, the oxide
cap, is IGLOO's), of another model than the material's, with a key its model does not read, on a
second row for one material, and on a row for no material of the phase file. A row of the kind
`breakup` is the exception: ICE has no breakup model (IGLOO has), so it ignores the row with a warning,
its model and keys unread, as it ignores `breakup=` on the material line.

`combustion` has no default: given in `[ICE-Physics]`, it holds for every material whose line
names no other word. The constants are the row's alone: `q-comb`, `d-stop` and `d-burnout` in
`[ICE-Physics]` stop the run and name the models file. The row of every burning material gives
`q-comb`, the heat of reaction per kg burnt (`q-comb=0` written is allowed). The burning
temperature of `direct` is `T-burn` of its row, else `boiling-temperature`, which then stands for
the boiling temperature at one atmosphere: with `boiling-pressure` at another value, the row must
give `T-burn`. Its latent heat
is the `Lv` column of the material's property table, read at the burning temperature, else
`latent-heat`. The table owns the property: with the column, leave `latent-heat` out or give 0 for
that material (a value beside the column is refused).

With `d-stop` the burning stops at that diameter and the residue stays. Without it the drops burn
out. In a time-accurate run, while they burn, the drops of a cell are removed when they are at
or below `d-burnout` (absent: a drop mass of 1e-20 kg) or when the step burnt as much as is left,
and what they still hold is kept as an amount that a coupler hands to its gas
([Particle Physics](../theory/physics.md#without-a-residue)). A steady run burns, at the exit of
a cell whose drops burn, a fraction of what leaves it that rises from 0 at the removal size to all
of it at 0.90 times it; the setup prints this band, `burnout_band`, for such a run. A drop that
does not burn is not removed, whatever its size. To compare a material with its run in IGLOO, give `d-burnout` in both
decks, at about 0.1 of the injected diameter: the removal then hands a thousandth of the
drop's mass.

`d-gen` (> 0) is the diameter the drops of the mass the material generates are born with. The generated mass is the
condensed oxide of the burnt metal, born in another ICE family: a coupler that knows the reaction
(hydra-MI2) names that family and the kg of generated mass per kg burnt. ICE running alone generates no mass
([Particle Physics](../theory/physics.md#generated-mass)).

The model of every material in `[ICE-Physics]`, its constants in a row under `-`:

```ini
[ICE-Physics]
combustion          = direct
latent-heat         = 10.8e6
boiling-temperature = 2791
```

```
A combustion - q-comb=9.53e6 d-stop=6e-5
```

The gas phase file must declare one species: `INPUT/phase.txt`, or
`INPUT/<gas-phase>-phase.txt` with `gas-phase` in `[ICE-Parameters]`. Written by hand:

```
ideal-gas phase
GAS 30.14
```

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
