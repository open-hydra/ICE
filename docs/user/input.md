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
| `[ICE-Physics]` | Drag and heat-transfer models, and the condensed-material properties |
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

### Choosing the drag and heat correlations

`drag` and `heat-transfer` are set once in `[ICE-Physics]` and apply to every family.
They are only consulted when the run is coupled (that is, when `INPUT/gas.tec` exists);
an uncoupled run ignores them, and leaving them at `none` in a coupled run stops the
solver with the list of valid names.

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

## Stopping conditions

`iter-threshold`, `time-threshold` and `res-threshold` are independent, and the run
stops at whichever is met first. Their defaults are effectively infinite for the first
two and $10^{-10}$ for the third, so a file that sets none of them will stop on the
residual — which for a genuinely unsteady problem may be immediately. Set
`res-threshold = 0` to disable it.

## Re-reading at runtime

If `ini-diter` is set, `input.ini` is re-read every that many iterations. Thresholds
and output frequencies then take effect during a run; anything consumed once at setup
does not.
