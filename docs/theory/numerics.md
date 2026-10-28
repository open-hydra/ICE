# Spatial Discretization

## Finite-volume framework

ICE uses a cell-centred finite-volume method on multi-block structured hexahedral
grids. For cell $i$,

$$
\frac{d \mathbf{U}_i}{d t}
= -\frac{1}{\Omega_i} \sum_{f \in \partial \Omega_i} \hat{\mathbf{F}}_f \cdot \hat{\mathbf{n}}_f \, A_f
  + \mathbf{S}_i ,
$$

with $\Omega_i$ the cell volume, $A_f$ the face area and $\hat{\mathbf{n}}_f$ the outward
unit normal. Face areas, normals, cell volumes and the three cell lengths are computed
once at setup from the node coordinates and stored. The volume of a cell is the integral
of the Jacobian determinant of the trilinear map of its eight nodes,
$\Omega = \frac{1}{6} \sum_d \mathbf{e}_d \cdot (\mathbf{S}_d^- + \mathbf{S}_d^+)$, with
$\mathbf{e}_d$ the mean edge of direction $d$ and $\mathbf{S}_d^\mp$ the area vectors (half
the cross product of the diagonals, both along $+d$) of its two faces of that direction:
exact on every cell, twisted faces included, so that the volumes of the cells add up to
the volume their boundary faces enclose. A cell whose volume is zero, or has not the sign
of the determinant of its mean edges (it is then turned over somewhere), stops the setup.

Each block carries **two ghost layers** on every face. The dimensionality is inferred
from the blocks: the mesh is 3-D if a block has more than one cell in $k$, 2-D if not and
a block has more than one cell in $j$, 1-D otherwise. Nothing has to be declared — the
flux loops over a direction with a single cell simply do no work. A 2-D mesh is an
axisymmetric wedge when, in every block, the $k = 1$ node plane is the $k = 0$ plane turned
about the $x$ axis by the same angle ($x$ and $r = \sqrt{y^2+z^2}$ kept), and plane when, in
every block, it is the $k = 0$ plane translated. Every node is checked against 1 % of the
thickness of its block, the largest move of a node from one plane to the other; a mesh
with blocks of both kinds, or a block that fits neither to 1 % (coordinates written with
too few digits, or $k$ planes neither parallel nor turned about $x$), stops at setup, and
the setup summary names the type. The ghost nodes, one layer beyond every face, are
extrapolated linearly along the index lines, $\mathbf{x}_{-1} = 2\mathbf{x}_0 - \mathbf{x}_1$,
so that a ghost cell has the length across the face of the first cell, whatever the
stretching (the length the reconstruction at the faces next to a boundary reads); beyond
the side faces of a wedge they are rotated by its angle instead.

The interior faces of a block are swept in two passes, odd faces then even, so that two
faces sharing a cell never accumulate into it at the same time. This is what makes the
OpenMP result independent of the thread count. Boundary faces are handled separately,
from the ghost values.

## Reconstruction

The left and right states at a face come from a piecewise-linear reconstruction of the
**primitive** variables, limited component by component:

$$
\mathbf{P}_L = \mathbf{P}_i + \beta\,\phi(s_{i}, s_{i-1})\,\delta_L, \qquad
\mathbf{P}_R = \mathbf{P}_{i+1} - \beta\,\phi(s_{i+1}, s_{i})\,\delta_R,
$$

where $s$ are the one-sided slopes on the non-uniform mesh, $\phi$ the limiter and $\beta$ the shock-detector weight described below.

A family of a solidifying material reconstructs, in place of $T_p$, its liquid-branch temperature $T_\ell = e/c_l$,
which is continuous where the content nucleates, and takes each face's $T_p$ and $f$ from the face's $T_\ell$ and
$\chi$ by the rule of a cell ([Solidification](physics.md#solidification)); a liquid face is unchanged, and a cold solid,
whose $T_\ell$ is negative, gets its physical temperature back.

### Limiters

Several options are available to compute $\phi$:

- minmod
- Van Albada
- Van Leer
- MC
- Superbee

### Shock detector

The shock-detector adds a sensor on density,

$$
s = \max_{d}\ \left|\frac{\rho_{d+1} - 2\rho + \rho_{d-1}}{\rho_{d+1} + 2\rho + \rho_{d-1}}\right|,
$$

from which the weight multiplying every limited slope is

$$
\beta = \begin{cases}
1 - \tanh\!\big(10\,\varphi^3\big),\quad \varphi = s\,\Delta, & s < 1/\Delta \\[2pt]
0, & \text{otherwise}
\end{cases}
\qquad \Delta = 20 .
$$

So $\beta = 1$ in smooth flow, where the scheme is the plain MUSCL one, falls continuously with the sensor, and is 0 at
a shock ($s \ge 1/\Delta$), where the scheme drops to first order. A face inside a block takes the smaller weight of its
two cells, $\min(\beta_L, \beta_R)$, and so does a face between two connected blocks (`101`, `201`), which lies inside the
domain just the same ([Boundary faces](#boundary-faces)). A face on a physical boundary, and a chimera face (`102`), whose
ghost is a blend of donor cells rather than a cell, take their own cell's weight.

The form is that of Tramel, Nichols and Buning (AIAA 2009-3988, Eqs. 29-30), with one change. There the weight is
$1 - \tanh(10\,\varphi^3)$ with $\varphi = s/\Delta$, and $\Delta$ between 1 and 10 was found adequate. With
$\varphi = s/\Delta$ and $\Delta = 20$, $\tanh(10\,\varphi^3)$ stays below $1.6\times10^{-7}$ for every $s$ below the
threshold $1/\Delta$, so $\beta$ is a switch: 1 below the threshold and 0 above it. ICE takes $\varphi = s\,\Delta$, which
keeps the threshold and spreads the fall of $\beta$ over $0 \le s < 1/\Delta$, from 1 at $s = 0$ to $1 - \tanh 10$
($4\times10^{-9}$) at the threshold. The paper also thickens the sensor (its Eq. 31 takes, for each cell, the largest $s$
of the 27 cells around it); the face weight $\min(\beta_L, \beta_R)$ gives each face the weight of the more shocked of the
two cells it separates, and treats a face the same from both sides. With the switch and the weight of the cell on the
face's low-index side, the crossing jets of MK and AG, mirror images of each other, were $1.0\times10^{-2}$ and
$1.7\times10^{-2}$ away from their own mirror image at iteration 5000; MK is now symmetric to round-off. AG keeps an
asymmetry of $3\times10^{-5}$: its steady state with the detector is one of a pair of steady states that are each other's
mirror image, and round-off decides which one a run reaches.

A switch also kept steady computations from converging: the cells near the threshold switched from one iteration to the
next, so the density residual levelled off (at about $2\times10^{-3}$ for MK and $4\times10^{-3}$ for AG on the crossing
jets) instead of falling. With the continuous weight both converge: MK reaches the default `res-threshold` at iteration
2493 and AG at iteration 6907. Below the threshold the weight is smaller than the switch's wherever the sensor is not
zero, so the scheme is more dissipative there: on the AG crossing jets the RMS distance to the exact solution grows from
$2.7\times10^{-2}$ to $2.8\times10^{-2}$ kg/m³, and the field's maximum, which lies on a jet (0.1 kg/m³ in the exact
solution), falls from 0.16, where the jets meet, to 0.12, just upstream of the crossing. A domain split into connected
blocks gets the weights of the same domain in
one block, so the two differ by round-off only — which the detector amplifies where the crossing jets collide: the
four-block crossing jets, cut at $y = 0$ through the collision, differ from the single block by $3\times10^{-10}$ of the
peak density at iteration 5000, less than a change of one inlet value in its last digit moves the single block
($1\times10^{-9}$; [Multi-block and MPI](../vv/multiblock-mpi.md)). Blocks joined by chimera faces still differ where
the detector acts at the interface.

On the [crossing jets](../vv/crossing-jets.md) the detector lets the MK and AG runs converge, keeps the AG pressure
tensor off its floor in the nearly empty cells around the jets, and keeps AG from gathering mass on either side of the
crossing. It costs accuracy where the flow is smooth: without it MK and AG are nearer their exact solutions, but neither
converges (MK's residual stalls and its jets are not each other's mirror image at the edges of the inlets; AG's
residual cycles). The crossing-jets tests run without it, for all three closures.

## Riemann solvers

The face flux is a two-state flux solved via three possible schemes.

| Name | Form | Applies to |
|---|---|---|
| **Saurel** | Sign of the mean normal velocity selects the donor state; there is no pressure and no sound speed to upwind against | MK only. It assembles the flux from the monokinetic variable layout, so ICE refuses it for IG and AG |
| **Rusanov** | $\tfrac12(\mathbf F_L + \mathbf F_R) - \tfrac12 A\,(\mathbf U_R - \mathbf U_L)$, with $A$ the largest of $|u_n \pm a|$ on the two sides and $a$ the signal speed across the face: $\sqrt{3P/\rho_p}$ for IG, $\sqrt{3P_{nn}/\rho_p}$ for AG | Any closure; the default for IG and AG |
| **HLLE** | Two-wave solver with Roe-averaged speed estimates, falling back to the upwind flux when both waves run the same way | Any closure with a sound speed |

The time step adds up the same directional speed over the directions the mesh resolves ([Time Integration](time-integration.md#the-time-step)).

HLLE is less dissipative than Rusanov on a contact and is worth trying when a contact is being smeared, at the cost of a Roe average per face.

## Boundary faces

Boundary fluxes are built from the ghost values, so every boundary type — connection, chimera, symmetry, extrapolation, inlet —
reaches the flux loop through the same path. The ghost fill is described under [Boundary Conditions](../user/boundary-conditions.md).
On the side faces of an axisymmetric wedge the ghost is the cell's mirror image, so the flux carries no mass, and for IG and
AG the pressure on those faces is the hoop term of the radial momentum; MK has no flux there.

A face between two blocks (`101`, `201`) has a record in each block. Its flux is computed once, from the owner record — the
one on the face's low-index side when the face is the high face of one block and the low face of the other, otherwise the
first in (block, face, k, j, i) order — with that block's stencil and, as for a face inside a block, the smaller of its two
cells' detector weights; the other block's cell takes its negative. The two blocks therefore exchange the same mass,
momentum and energy across the cut, which a reconstruction from each side would not give: each block would take the lengths
of the other's cells from its own extrapolated ghost nodes. With MPI the rank of the other block sends its cell's weight to
the owner's rank, and the owner's rank sends back the flux, once per stage.

The state an inlet injects is the state of the boundary face, not of the ghost cell's centre: the reconstructions take the
length of an injecting ghost as zero, so the first cell's slope toward the inlet is the one-sided slope from the face, on
the face itself and on the first interior face, and the scheme stays second order at the inlet. Where the particles leave
through an inlet face, or its mass flux is zero, the ghost is a copy or a mirror of the cell and keeps its length.

A reconstructed state that is unphysical has its slopes halved until it is not. If it is still unphysical at first order —
a ghost that nothing filled, a NaN in the stencil — ICE prints the stencil and stops with a non-zero status.

## Source terms

Drag, convective and radiative heat, and evaporation are evaluated **explicitly**, once per Runge-Kutta stage, from the current primitive state, and added to the residual before it is integrated. There
is no point-implicit or operator-split treatment.
