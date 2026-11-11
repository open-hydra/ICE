# Boundary Conditions

ICE reads its boundary conditions from `INPUT/part-bc.txt`, a table written by
[ATLAS BCB](https://open-hydra.github.io/ATLAS/user-guide/bcb/). The type of each face
is chosen in ATLAS; ICE applies the table.

Every block has two layers of ghost cells around it. A boundary condition sets the
ghost cells, and the flux through the boundary face is then computed from them with
the same reconstruction and Riemann solver as an interior face.

## The codes

| Code | Type | What ICE does |
|------|------|---------------|
| `0` | None | The ghost copies the interior cell, and no flux crosses. For the faces of a direction the mesh does not resolve, such as faces 5 and 6 of a planar 2-D case |
| `101` | Connection | The ghost cells come from the neighbouring block; the flux through the face is computed once, for both blocks |
| `201` | Periodic | As `101` |
| `102` | Chimera | Each ghost layer is a weighted blend of donor cells of an overlapping block ([below](#chimera)) |
| `200` | Axisymmetry | The side faces of an axisymmetric wedge ([below](#axisymmetry-200)) |
| `300` | Symmetry | Reflects the cloud, or lets it leave ([below](#symmetry-and-walls-300-301)) |
| `301` | Wall | As `300` |
| `400` | Extrapolation | The ghost copies the interior cell: a zero-gradient outflow ([below](#extrapolation-400)) |
| `401` | Inlet, as ratios to the gas | [Inlets](#inlets) |
| `402` | Inlet, absolute values | [Inlets](#inlets) |
| `403` | Inlet, the speed as a ratio to the gas | [Inlets](#inlets) |

Any other code stops ICE at setup.

## File format

There is one record per boundary **cell**. A record starts with the line

```
block  i  j  k  face  code
```

followed by the data lines its code needs. The faces are numbered 1 and 2 for the low
and high $i$ sides, 3 and 4 for $j$, 5 and 6 for $k$. The cell is the interior cell
against the face, so a face-2 record has $i = n_x$, a face-4 record $j = n_y$.

A block has $2\,(n_y n_z + n_x n_z + n_x n_y)$ boundary cells, and the table holds one
record for each. A header line whose sixth column is not `0` or a three-digit code
stops the run, naming the record: that is what a file in an older format looks like.

### Several families

A phase can have several materials, each with several populations
(`<material> <populations>` in the ATLAS phase file). ATLAS then writes the whole table
once per (material, population) pair: the copies differ only in the data of their inlet
records. ICE matches them with its families:

| Copies in the file | What ICE does |
|---|---|
| 1 | Every family uses it. This is the single-population case. |
| One per family | Copy *c* belongs to family *c*. ICE identifies the copies by their order, so the file must keep ATLAS's: blocks in mesh order, the copies of a block together, each copy over the same cells. |
| Any other number | The setup stops: binding a wrong copy would inject the wrong particles without any other sign. |

For example, a phase with two materials, the first in two populations, needs three
`[ICE-FamilyN]` sections; family 3 is the second material.

## Symmetry and walls (`300`, `301`)

The ghost copies the interior cell with its velocity mirrored about the face (for AG,
the dispersion tensor too, $\mathsf{P}' = H\mathsf{P}H$ with
$H = I - 2\hat{\mathbf n}\hat{\mathbf n}^T$). A cloud that reaches the face is
reflected.

A cloud moving **out** through the face is let go instead: the ghost is not mirrored,
and the face acts as an outflow. The exception is face 3, the axis of an axisymmetric
case, which always mirrors.

## Axisymmetry (`200`)

A 2-D axisymmetric mesh is one layer of cells rotated about the $x$ axis; its two side
faces (5 and 6) carry `200`. The ghost is the mirror image of the interior cell: for a
flow without swirl, this is exactly the neighbouring wedge rotated into place. For IG
and AG, no mass crosses these faces, and the pressure on them supplies the hoop term of
the radial momentum. A pressureless MK cloud has no flux through them.

The axis (face 3) mirrors whether it carries `200` or `300`.

## Extrapolation (`400`)

The ghost copies the interior cell, so the gradient across the face is zero: the cloud
leaves freely.

Particles can also *enter* through such a face, with the state of the boundary cell.
Since that is rarely intended, ICE reports it:

- a line the first time particles enter through the `400` faces of a block face (per
  family and grid level), with the iteration and how fast they enter compared with the
  signal speed of the cloud;
- another line each time that speed grows more than tenfold;
- at the end of a standalone run, which faces let particles in.

Only cells that hold particles (density above $10^{-6}$ kg/m³) moving in faster than
round-off count.

## Inlets (`401`–`403`) { #inlets }

An inlet prescribes a mass flux, a direction, a temperature and a particle radius.
The three codes differ in what is absolute and what is a ratio to the local gas:

| Code | Data | Mass flux | Speed | Temperature |
|---|---|---|---|---|
| `401` | `k kV alpha beta kT rp`, then `G` | ratios to the gas flow, plus `G` (below) | `kV` × the gas speed | `kT` × the gas temperature |
| `402` | `gp v alpha beta Tp rp` | `gp`, in kg s⁻¹ m⁻² | `v`, absolute | `Tp`, absolute |
| `403` | `gp kV alpha beta Tp rp` | `gp`, in kg s⁻¹ m⁻² | `kV` × the gas speed | `Tp`, absolute |

- **Direction.** `alpha beta` written as `normal normal` (the usual choice) inject
  along the inward normal of the face. Written as numbers, they are angles in the global
  frame: $\alpha$ in the $x$–$y$ plane from $x$, and $\beta$ out of it, giving the
  direction $(\cos\beta\cos\alpha, \cos\beta\sin\alpha, \sin\beta)$. One `normal` and
  one number is refused.
- **Number density.** It follows from the radius `rp`:
  $n = \rho_p / \big(\rho_{al}(T_p)\tfrac43\pi r_p^3\big)$, with the material density
  at the inlet temperature (from the property table, or the constant `density`).
- **Dispersion.** The Gaussian closures inject a dispersion of $10^{-6}$: an
  effectively monokinetic stream.
- **A reversed inlet.** Where the particles at the face are already leaving, the inlet
  acts as an extrapolation, so it never injects against an outflow.
- **A zero mass flux.** The face acts as a symmetry: this is how a partly active inlet
  face is written.
- **Without gas.** `401` and `403` take values from the gas, so an uncoupled run
  refuses them; use `402`.

### The mass flux of `401`

`401` shares the gas flow among all the dispersed phases of the cell, as ATLAS BCB
computes it. With $krho_{tot}$ the sum of `krho` over every `401` population of the
cell and $G$ the sum of `gp` over its `402` and `403` ones, ATLAS writes
$k = krho/(1-krho_{tot})$ as the first value and $G$ as the tenth. Family $i$ then
injects the mass flux

$$
k_i\,(\rho_g u_{g,n} + G)
$$

along the normal. With the direction given as angles, the gas part
carries the cosine between the direction and the normal, and $G$ does not. $G$ only
counts where it is positive and the injection moves into the domain. A `401` line without the tenth value comes from an
older ATLAS: regenerate the file with ATLAS BCB. A negative $k$, or a $G$ that is
negative or not finite, stops the setup.

## Chimera (`102`) { #chimera }

ATLAS BCB finds the cells of an overlapping block that cover each ghost layer, and
writes their volume weights, which add up to 1 per layer. ICE blends these donor cells
in conservative variables, so the interpolation conserves mass, momentum and energy.

ICE stops at setup if the weights of a layer add up to nearly zero: ATLAS found no
donor for it. Partial coverage is **not** caught, because ATLAS normalises the weights
over the donors it did find. Check `chimera.log` from the BCB run.

[Chimera Overset](../vv/chimera.md) shows the accuracy this reaches.

## Directions with a single cell

Faces 5 and 6 of a 2-D mesh, and faces 3 to 6 of a 1-D one, close a direction the mesh
does not resolve. When both faces of such a direction mirror the flow (`200`, `300` or
`301`), the velocity across them can grow from round-off on a thin layer until the run
fails. Since that component is zero anyway between two mirror planes, ICE sets it to
zero at setup, does not integrate it, and says so:

| Mesh | Faces | Not integrated (AG also) |
|---|---|---|
| 2-D plane | 5 and 6, each a plane normal to $z$ | $w$ ($P_{13}$, $P_{23}$) |
| 2-D axisymmetric wedge | 5 and 6 all `200`, symmetric about $z = 0$ | $w$, the swirl ($P_{13}$, $P_{23}$) |
| 1-D | 3 and 4, planes normal to $y$ | $v$ ($P_{12}$, $P_{23}$) |
| 1-D | 5 and 6 | $w$ ($P_{13}$, $P_{23}$) |

A wedge whose side faces carry another code, or that is not symmetric about $z = 0$,
stops the setup. Mirror faces that are not planes normal to their axis, or a mirror on
one side only, leave the component integrated, with a warning: give both faces the
same code.

??? note "The second ghost layer"
    `101` and `201` copy it from the neighbour, `102` blends its own donors. On `200`,
    `300`, `301` and inlet faces it follows the first layer's rule: the mirror image of
    the second interior cell, a copy of it, or the injected state again. `400` builds it
    by second-order extrapolation, $P_{g2} = 3P_{g1} - 3P_m + P_{m+1}$, which gives the
    second interior cell again. `0` leaves it unset, since no stencil reaches it.

## File names

All of ICE's input files carry the prefix `part-`. A coupled Hydra run gives each
phase its own prefix, through the `phase` key of `[ICE-Parameters]`.
