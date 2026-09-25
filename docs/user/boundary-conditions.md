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

records per block, covering every boundary cell of every face once, in any order. A
record applies to all particle families.

A header whose sixth column is neither `0` nor a three-digit code stops the run with
the record number — which is what a file in an older, wider format looks like when it
is read with this schema.

## Supported codes

| Code | Type | Payload | ICE treatment |
|------|------|---------|---------------|
| `0` | Null | none | The ghost copies the interior cell; no boundary flux. This is the code for the faces of a direction the mesh does not resolve — faces 5 and 6 of a planar 2-D case |
| `101` | Connection | `b2 i2 j2 k2 f2 c1 c2 c3 c4` | Both ghost layers copied from the interior cells of the neighbouring block |
| `201` | Periodic | as `101` | As `101` |
| `102` | Chimera (overset) | donor counts per layer, then `b i j k weight` per donor | Each ghost layer set to the weighted blend of its own donors, in conservative variables |
| `200` | Axisymmetry | none | The side faces of an axisymmetric wedge: see below |
| `300` | Symmetry | none | See below |
| `400` | Extrapolation | none | Ghost copies the interior cell: zero gradient |
| `401` | Inlet, ratios to the gas | `krho kV alpha beta kT rp` | See below |
| `402` | Inlet, absolute | `gp v alpha beta Tp rp` | |
| `403` | Inlet, speed as a ratio | `gp kV alpha beta Tp rp` | |

Any other code stops ICE at setup.

The second ghost layer is copied from the neighbour for `101`/`201`, taken from its own
donors for `102`, the mirror image of the second interior cell for `200`, and otherwise
built by a second-order extrapolation $P_{g2} = 3P_{g1} - 3P_m + P_{m+1}$; `0` leaves it
unfilled, since no stencil reaches it.

## Symmetry (`300`)

The ghost copies the interior cell and its velocity is mirrored about the face normal —
except when the particles are moving *out* through the face, in which case the copy is
left unmirrored and the face behaves as an outflow. A cloud reaching a symmetry plane
is therefore reflected, while one leaving through it is allowed to go.

Face 3 is the exception: it always mirrors, because it is conventionally the axis of an
axisymmetric case, where letting the cloud leave would be wrong.

## Axisymmetry (`200`)

A 2-D axisymmetric mesh is one layer of cells rotated about the $x$ axis, and its two side
faces (5 and 6) carry `200`. The ghost is the mirror image of the interior cell, with the
velocity reflected about the face; for a state without swirl this is exactly the
neighbouring wedge rotated into place. The second ghost mirrors the second interior cell.
For IG and AG the boundary flux is the Riemann flux between the cell and its mirror image:
no mass crosses, and the pressure on the two side faces supplies the hoop term of the
radial momentum. A pressureless MK cloud has no flux through these faces.

The axis face (face 3) of such a mesh mirrors whether it carries `200` or `300`.

## Inlets (`401`–`403`)

All three prescribe a mass flux, a direction, a temperature and a particle radius; they
differ in whether those are absolute or ratios to the local gas.

| | Mass term | Velocity | Temperature |
|---|---|---|---|
| `401` | `krho`, a loading ratio: $\rho_p = \frac{krho}{1-krho}\,\rho_g / kV$ | `kV` times the gas velocity normal to the face | `kT` times the gas temperature |
| `402` | `gp`, a mass flux in kg s⁻¹ m⁻²: $\rho_p = gp/|v_n|$ | `v`, an absolute speed | `Tp`, absolute |
| `403` | `gp`, as `402` | `kV` times the gas speed | `Tp`, absolute |

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

`401` needs a gas field to be meaningful, and `403` needs one for its velocity.

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
