# Boundary Conditions

ICE reads its boundary conditions from `INPUT/part-bc.txt`, written by
[ATLAS BCB](https://github.com/open-hydra/ATLAS). Boundary types are chosen per block
face in the ATLAS input; ICE only consumes the resulting table.

Every boundary type reaches the solver through the same path: `compute_ghost` fills the
first ghost layer, `fill_second_ghost` the second, and the boundary flux is then built
from the ghost values by the same reconstruction and Riemann solver as an interior
face. Each block has two ghost layers on every face.

## File format

One record per boundary **cell**, in ATLAS's layout — a header line

```
block  i  j  k  face  code
```

followed by zero or more payload lines depending on `code`. Faces are numbered 1 and 2
for the low and high $i$ sides, 3 and 4 for $j$, 5 and 6 for $k$; the cell index is the
interior cell against the face, so face 2 records carry $i = n_x$ and face 4 records
$j = n_y$.

ICE expects exactly

$$
2\,(n_y n_z + n_x n_z + n_x n_y)
$$

records per block per **copy** of the boundary table, covering every boundary cell of
every face once.

A header whose sixth column is neither `0` nor a three-digit code stops the run with
the record number — which is what a file in an older, wider format looks like when it
is read with this schema.

### Populations

A dispersed phase can carry several materials, each with several populations, declared
in the ATLAS phase file as `<material> <npCP>`. ATLAS writes the whole boundary table
once per (material, population) pair — block by block, the copies of a block together —
and the copies differ only in the payload of their `401`–`403` records. ICE binds them
to its `[ICE-Family*]` sections:

| Copies in the file | What ICE does |
|---|---|
| 1 | The record applies to every family, which is what a single-population case and every case written before populations look like. The file may be in any order. |
| one per family | Copy *c* belongs to family *c*. The order is what identifies a copy, so the file must keep ATLAS's: blocks in mesh order, the copies of a block together, each copy over the same cells in the same order. |
| anything else | The setup stops. Neither reading is right, and binding the wrong copy to a family would inject the wrong size class with no other symptom. |

So a phase with two materials, one of them in two populations, needs three
`[ICE-Family*]` sections, and family 3 is the second material's single population.

## Supported codes

| Code | Type | Payload | ICE treatment |
|------|------|---------|---------------|
| `0` | Null | none | The ghost copies the interior cell; no boundary flux. This is the code for the faces of a direction the mesh does not resolve — faces 5 and 6 of a planar 2-D case |
| `101` | Connection | `b2 i2 j2 k2 f2 c1 c2 c3 c4` | Both ghost layers copied from the interior cells of the neighbouring block; the face's flux is computed once and the neighbouring block takes its negative |
| `201` | Periodic | as `101` | As `101` |
| `102` | Chimera (overset) | donor counts per layer, then `b i j k weight` per donor | Each ghost layer set to the weighted blend of its own donors, in conservative variables |
| `200` | Axisymmetry | none | The side faces of an axisymmetric wedge: see below |
| `300` | Symmetry | none | See below |
| `400` | Extrapolation | none | Ghost copies the interior cell: zero gradient. Particles that move into the domain through it enter with the state of the boundary cell, and the run says so: a line the first time a cell that carries the phase (density above the empty-cell density, $10^{-6}$ kg/m³) moves into the domain faster than round-off (above $10^{-12}$ of the signal speed $\sqrt{P/\rho}$ for IG and AG, of its own speed for MK, which has no signal speed) through the `400` faces of a face of a block (per family and grid level), with the iteration, the faces and the largest $\lvert w\rvert$ ($w$ the velocity along the normal out of the domain; $\lvert w\rvert/\sqrt{P/\rho}$ for IG, $\lvert w\rvert/\sqrt{P_{nn}/\rho}$ for AG), another each time that largest value is more than 10 times that of the last line, and at the end of a standalone run the stages each such face let particles in (under MPI the rank that owns the block prints its lines, at once, and rank 0 the end lines; inside hydra only the first-inflow and 10-times lines print, not the end lines, since hydra never runs ICE's end of run) |
| `401` | Inlet, ratios to the gas | `krho kV alpha beta kT rp` | See below |
| `402` | Inlet, absolute | `gp v alpha beta Tp rp` | |
| `403` | Inlet, speed as a ratio | `gp kV alpha beta Tp rp` | |

Any other code stops ICE at setup.

The second ghost layer is copied from the neighbour for `101`/`201` and taken from its own
donors for `102`. On a `200`, `300`/`301` or inlet face it follows the rule of the first:
the mirror image of the second interior cell where the first ghost mirrors the first, a copy
of it where the first is a copy, and the injected state again where an inlet injects. `400`
builds it by a second-order extrapolation $P_{g2} = 3P_{g1} - 3P_m + P_{m+1}$, which from its
copied first ghost is the second interior cell too; `0` leaves it unfilled, since no stencil
reaches it.

## Symmetry (`300`)

The ghost copies the interior cell and its velocity is mirrored about the face normal (for
AG the dispersion tensor too, $\mathsf{P}' = H\mathsf{P}H$ with $H = I - 2\hat{\mathbf n}\hat{\mathbf n}^T$) —
except when the particles are moving *out* through the face, in which case the copy is
left unmirrored and the face behaves as an outflow. A cloud reaching a symmetry plane
is therefore reflected, while one leaving through it is allowed to go.

Face 3 is the exception: it always mirrors, because it is conventionally the axis of an
axisymmetric case, where letting the cloud leave would be wrong.

## Axisymmetry (`200`)

A 2-D axisymmetric mesh is one layer of cells rotated about the $x$ axis, and its two side
faces (5 and 6) carry `200`. The ghost is the mirror image of the interior cell, with the
velocity (and for AG the dispersion tensor) reflected about the face; for a state without swirl this is exactly the
neighbouring wedge rotated into place. The second ghost mirrors the second interior cell.
For IG and AG the boundary flux is the Riemann flux between the cell and its mirror image:
no mass crosses, and the pressure on the two side faces supplies the hoop term of the
radial momentum. A pressureless MK cloud has no flux through these faces.

The axis face (face 3) of such a mesh mirrors whether it carries `200` or `300`.

## Faces of a direction with one cell

Faces 5 and 6 of a 2-D mesh, and faces 3 to 6 of a 1-D one, bound a direction the mesh does
not resolve. Where both faces of such a direction mirror the flow (`200`, `300` or `301`), the
components of the state that are odd across them have a mode the time step does not follow:
between two mirror planes $\Delta z$ apart, $dw/dt = -2a\,w/\Delta z$ ($a$ the signal speed of the
closure), which on a layer thin against its in-plane cells makes round-off in $w$ grow until the
run fails. The flow between two mirror planes has those components zero, so ICE sets them to zero
at set-up and does not integrate them, and says so in a warning:

| Mesh | Faces | Not integrated (with AG also) |
|---|---|---|
| 2-D plane | 5 and 6 all `200`, `300` or `301`, each $k$ node plane normal to $z$ | $w$ ($P_{13}$, $P_{23}$) |
| 2-D axisymmetric wedge | 5 and 6 all `200`, the $k$ node planes at $\mp\theta/2$ about $z = 0$ | $w$, the swirl ($P_{13}$, $P_{23}$) |
| 1-D | 3 and 4 as above, each $j$ node plane normal to $y$ | $v$ ($P_{12}$, $P_{23}$) |
| 1-D | 5 and 6 as above | $w$ ($P_{13}$, $P_{23}$) |

A plane is checked to 1 % of the thickness of its block. A wedge whose side faces carry any other
code, or that is not symmetric about $z = 0$, stops the setup. Mirror faces that are not planes
normal to their axis, as on a layer between parallel planes inclined to $z$, ask the cloud to
follow them across the layer: the component is then integrated, with a warning, and a layer thin
against its other cells can still fail. So can a layer with a mirror on one side only and `0` or
`400` on the other, where the mode is half as fast: give both faces the same code.

## Inlets (`401`–`403`)

All three prescribe a mass flux, a direction, a temperature and a particle radius; they
differ in whether those are absolute or ratios to the local gas.

| | Mass term | Velocity | Temperature |
|---|---|---|---|
| `401` | `krho`, the family's share of the mixture flow: $\rho_{p,i} = \frac{krho_i}{1-\sum_j krho_j}\,\rho_g / kV$ | `kV` times the gas velocity normal to the face | `kT` times the gas temperature |
| `402` | `gp`, a mass flux in kg s⁻¹ m⁻²: $\rho_p = gp/|v_n|$ | `v`, an absolute speed | `Tp`, absolute |
| `403` | `gp`, as `402` | `kV` times the gas speed | `Tp`, absolute |

In the `401` formula the sum runs over the families of the cell, as in IGLOO; with one copy of the table, its
record counts once for every family. A cell whose `krho` sum reaches 1, or a negative `krho`, stops the setup.

`alpha` and `beta` are the injection angles about the $x$ axis; written as the literal
`normal`, they are replaced by the angles of the face normal, which is the usual choice.
The number density follows from the prescribed radius,
$n = \rho_p / \big(\rho_{al}(T_p)\tfrac43\pi r_p^3\big)$, with the condensed density at the
inlet temperature: the property table's when the case has one, the constant `density` of
`[ICE-Physics]` otherwise. For the Gaussian closures the dispersion in
the ghost cell is set to $10^{-6}$, so an inlet injects an effectively monokinetic
stream.

Two cases are handled before any of that:

- if the particles at the face are already **leaving** the domain, the inlet becomes an
  extrapolation, so a wrongly-placed inlet does not inject against an outflow;
- if the mass flux is **zero**, the inlet becomes a symmetry, which is how a partly
  active inlet face is expressed.

`401` takes its density, velocity and temperature from the gas, and `403` its velocity: a run without a gas
phase (0-way) refuses an injecting `401` or `403` record at set-up; use `402` there.

## Chimera (`102`) { #chimera }

ATLAS BCB finds the donor cells overlapping each of the two ghost layers and writes
their volume weights, normalised to 1 per layer. ICE blends the donors in conservative
variables, converting back to primitives afterwards, so the interpolation conserves
mass, momentum and energy.

ICE stops at setup if any ghost layer's weights sum to nearly zero, which means ATLAS
found no donor for it. Partial coverage is *not* caught: ATLAS normalises over the
donors it did find, so a ghost cell covered by half a donor block still gets weights
summing to 1. Check `chimera.log` from the BCB run.

See [Chimera Overset](../vv/chimera.md) for the accuracy this reaches.

## File naming

All of ICE's input files carry the phase prefix `part-`, which is fixed for the
standalone solver; it is a variable only so that a coupled Hydra run can give each
phase its own set of files.
