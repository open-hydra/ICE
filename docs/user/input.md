# Input File

`input.ini` is the only configuration file ICE reads, and it is read from the working
directory. It is an INI file parsed by [FiNeR](https://github.com/szaghi/FiNeR):
sections in square brackets, `key = value` inside them, and everything optional unless
noted below.

!!! warning "Comments must start with `;`"
    Only a semicolon is recognised as a comment. A line starting with `#` or `!` is not
    treated as a comment and will corrupt the section it sits in — usually surfacing as
    a confusing complaint about an unrelated parameter.

Unknown sections and unknown keys are ignored in silence, so a misspelled key takes its
default rather than raising an error. The startup report is what to check against: it
prints the closures, the scheme and the boundary types that were actually selected.

## A minimal file

```ini
[ICE-Parameters]
cfl            = 0.5
time-accurate  = true
time-threshold = 0.125

[ICE-Scheme]
space-reconstruction = MUSCL
flux-limiter         = VANLEER
time                 = 2

[ICE-Family1]
model = AG

[ICE-IO]
bck-format  = tecplot ascii
sol-format  = tecplot ascii
shell-diter = 100
```

## Sections

| Section | Contents |
|---|---|
| `[ICE-Parameters]` | CFL number, time-step ceiling, time-accurate flag, residual smoothing, and the three stopping thresholds |
| `[ICE-Scheme]` | Reconstruction, limiter, Runge-Kutta stages, drag and heat correlations |
| `[ICE-FamilyN]` | One per particle family: its closure. **The only required key in the file** |
| `[ICE-IO]` | Output formats and frequencies, restart flag, gas-file directory |
| `[ICE-Multigrid]` | Number of grid levels and the iteration budget of each |
| `[ICE-Probes]` | Names the section that configures each probe |
| `[<probe name>]` | One per probe: location, variables, sampling frequency |
| `[ICE-Physics]` | Condensed-material density, specific heat, emissivity |

Every parameter, with its default and its validation rule, is listed in the
[parameter reference](registry.md), which is generated from the registry in the source.

### Families

The number of families is determined by counting `[ICE-FamilyN]` sections from 1
upwards, stopping at the first gap. Two families:

```ini
[ICE-Family1]
model = MK

[ICE-Family2]
model = AG
```

Each carries its own state; they share the mesh, the gas field and the numerical
scheme, and do not interact. `model` is required — a file with no family section stops
the run.

### Choosing the drag and heat correlations

`drag` and `heat` are set once in `[ICE-Scheme]` and apply to every family. They are
only consulted when the run is coupled (that is, when `INPUT/gas.tec` exists); an
uncoupled run ignores them, and leaving them at `None` in a coupled run stops the
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
