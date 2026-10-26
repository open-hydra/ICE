# Input File

`input.ini` is the only configuration file ICE reads, and it is read from the working
directory. It is an INI file parsed by [FiNeR](https://github.com/szaghi/FiNeR):
sections in square brackets, `key = value` inside them, and everything optional unless
noted below.

A comment line starts with `;`, `#` or `!`; `;` also works inline, trimming the rest of
the value.

A key ICE does not recognise, inside a section it does, stops the run and is named —
a misspelled or renamed key would otherwise take its default and change nothing
silently. Sections ICE does not own are left alone, so a case may carry the
mesh-generator sections that produced it.

The section names and most of the keys mirror
[MOSE](https://github.com/open-hydra/MOSE), so a case set up for one solver reads the
same way as a case set up for the other.

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

`time-scheme`, `cfl`, `time-accurate`, `space-reconstruction` and `closure` are
required; everything else has a default.

## Sections

| Section | Contents |
|---|---|
| `[ICE-Parameters]` | Restart flag and the three stopping thresholds |
| `[ICE-Numerics]` | Time scheme, CFL number and step ceiling, reconstruction, limiter, shock detector, Riemann solver, residual smoothing |
| `[ICE-Physics]` | Drag, heat-transfer and evaporation models, and the condensed-material and vapour properties |
| `[ICE-IO]` | Input and output formats and frequencies, gas-file directory |
| `[ICE-Multigrid]` | Number of grid levels and the iteration budget of each |
| `[ICE-Probes]` | Names the section that configures each probe |
| `[<probe name>]` | One per probe: location, variables, sampling frequency |
| `[ICE-FamilyN]` | One per particle family: its closure |

Every parameter, with its default and its validation rule, is listed in the
[parameter reference](registry.md), which is generated from the registry in the source.

### Families

The number of families is determined by counting `[ICE-FamilyN]` sections from 1
upwards, stopping at the first gap. Two families:

```ini
[ICE-Family1]
closure = MK

[ICE-Family2]
closure = AG
```

Each carries its own state; they share the mesh, the gas field and the numerical
scheme, and do not interact. `closure` is required — a file with no family section
stops the run.

### Materials

A family is made of one condensed material. `INPUT/<prefix>phase.txt`, which ATLAS writes,
names the materials after its type line, one line each, `<name> <groups> [key=value ...]`:

```
condensed-dispersed phase
A 1 evaporation=CEM
B 2
```

The families map onto the materials in that order: here family 1 is material A and
families 2 and 3 are material B. Without the file every family uses one material. With
several, the groups must add up to the `[ICE-FamilyN]` count and the
[property table](initial-conditions.md#property-table) must give one zone per material.

`density`, `specific-heat`, `latent-heat`, `emissivity`, `vapour-molar-mass`,
`boiling-temperature` (alias `Tboil`), `vapour-specific-heat`, `lewis-number`, `vapour-mass-fraction` and
`evaporation-coefficient` take one value per material, in the phase file's order
(`emissivity = 0 0` for two materials); any other count stops the run. `latent-heat`,
`vapour-molar-mass` and `vapour-specific-heat` are 0 for a material whose property table
carries the `Lv`, `Mv` or `Cpv` column of the same quantity: the table then owns it.

The tokens are the ones IGLOO reads. `evaporation`, `interface` and `alpha-e` set that
material's evaporation model, interface and accommodation coefficient; `[ICE-Physics]
evaporation`, `evaporation-interface` and `evaporation-coefficient` are the default of every
material without the token. `solidification = on` makes the material freeze and melt, with
`T-melt`, `h-fus`, `T-nuc` and `cp-solid` (see [Solidification](../theory/physics.md#solidification)); it
stops the run without a positive `h-fus` and `cp-solid`, with `T-nuc` at or above `T-melt`, together
with evaporation, and with a `Cp` or `Density` column that varies. `liquid-conduction`, `boiling`
and `combustion` accept only the value ICE implements (`ITC`, `clamp`, `none`), and IGLOO's other
numeric keys (`k-liq`, `mu-liq`, `K-burn`, …) are read and ignored. An unknown key, a value that is
not a number, or `interface = LK` with the d2-law stops the run. `evaporation-blowing` stays global.
The setup prints each material's evaporation model and each solidifying material's values.

### Choosing the exchange models

`drag` and `heat-transfer` are set once in `[ICE-Physics]` and apply to every family;
`evaporation` is the default of every material (see [Materials](#materials)). They are only consulted when the run is coupled (that is, when
`INPUT/gas.tec` exists); an uncoupled run ignores them, and a name none of them
recognises stops the solver with the list of valid ones. In a coupled run `drag` and
`heat-transfer` are required: leaving either at its default `none` stops the solver with
that list. `drag = NoDrag` and `heat-transfer = NoHeat` switch the momentum and the
convective heat exchange off explicitly; radiation stays under `emissivity`.

`evaporation` defaults to `none`, and while it is `none` the vapour keys beside it are
never read. Selecting a model makes `latent-heat`, `vapour-molar-mass` and
`boiling-temperature` matter — those three set the saturation curve unless the property
table carries a `Psat` column, and their defaults describe aluminium. A table with `Lv`,
`Mv` and `Cpv` columns supplies the latent heat, the molar mass and the vapour specific
heat in place of `latent-heat`, `vapour-molar-mass` and `vapour-specific-heat`, the first
and the last as functions of the particle temperature (see the
[property table](initial-conditions.md#property-table)). `boiling-temperature`
may be given as `Tboil`, the name IGLOO also accepts; giving both stops the run. `evaporation-interface` and `evaporation-blowing` are refinements of
the selected model rather than models of their own:

```ini
[ICE-Physics]
evaporation           = CEM
evaporation-interface = LK
evaporation-blowing   = LK
latent-heat           = 2.26e6
vapour-molar-mass     = 18.015
boiling-temperature   = 373.15
```

### Grid levels

```ini
[ICE-Multigrid]
levels      = 2
level1-iter = 5000
level2-iter = 500
```

`levels = 2` adds one coarse grid at half the resolution in every direction, so every
block dimension must be divisible by 2. The run starts on the coarsest level; see
[Grid sequencing](../theory/time-integration.md#grid-sequencing).

Each level is a mesh in its own right and needs **its own boundary file**, because the
records are per boundary cell and a coarse cell is not a fine one:

| Level | File |
|---|---|
| 1 (fine) | `INPUT/<prefix>bc.txt` |
| 2 | `INPUT/<prefix>bc2.txt` |
| *n* | `INPUT/<prefix>bc<n>.txt` |

ATLAS BCB writes all of them in one pass when `MG-levels` is set in its own
`[ATLAS-Parameters]`, so the two settings have to agree. ICE checks each file's record
count against that level's mesh and stops if they disagree, so a stale or missing file
is reported rather than absorbed.

## Stopping conditions

`iter-threshold`, `time-threshold` and `res-threshold` are independent, and the run
stops at whichever is met first. `res-threshold` ends steady runs only
(`time-accurate = false`): the density change of a time step says nothing about
convergence, so a time-accurate run ends on the other two. Their defaults are
effectively infinite for the first two and $10^{-10}$ for the third, so a steady run
that sets none of them stops on the residual. Set `res-threshold = 0` to disable it.

`res-threshold` is relative: a steady run ends when its density residual (the L2 norm
over the cells of the change of density in one iteration, the third column of
`residual-history.dat`, after the iteration and the time) has fallen by that factor from the first non-zero one, on the
finest level. The residual scales with the loading and with the square root of the
number of cells, so a relative test stops the same flow at the same iteration whatever
its density, and a dilute cloud does not stop at its first iteration. The reference is
the first non-zero residual rather than the first one, because a cloud at rest that the
gas starts to move can have no density change at its first step (none with the one-stage
`euler` scheme); while every residual is zero, the residual ends nothing, and
`iter-threshold` ends a run in which nothing moves. A `res-threshold` switched on by a
re-read of `input.ini` (from 0) measures from the first non-zero residual after it. If a
residual rises before it falls, the run ends when it is `res-threshold` times the first
one, later than a test relative to the peak would. A restart takes the first non-zero
residual of its own run as the reference (the solution files do not carry the first
run's) and says so. An input written for the former absolute test keeps its stop with
`res-threshold` divided by that first residual.

## Re-reading at runtime

If `ini-diter` is set, `input.ini` is re-read every that many iterations. Thresholds
and output frequencies then take effect during a run: an `iter-threshold` or
`level<n>-iter` lowered below the iterations already done ends that level at the next
iteration. Anything consumed once at setup does not: the scheme, and the materials with
their properties, models and property table. A change of `levels` or of a family's
`closure`, which setup has sized the run for, stops the run, as does a value that breaks
its rule; `irs = true` in a time-accurate run is ignored, with a warning.
