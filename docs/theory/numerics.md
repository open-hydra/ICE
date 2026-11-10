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
unit normal. Face areas, normals, volumes and the three cell lengths are computed once, at
setup, from the node coordinates.

The volume of a cell is the integral of the Jacobian determinant of the trilinear map of its
eight nodes,
$\Omega = \frac{1}{6} \sum_d \mathbf{e}_d \cdot (\mathbf{S}_d^- + \mathbf{S}_d^+)$, with
$\mathbf{e}_d$ the mean edge of direction $d$ and $\mathbf{S}_d^\mp$ the area vectors (half
the cross product of the diagonals) of its two faces of that direction. It is exact on every
cell, twisted faces included, so the cell volumes add up to the volume their boundary
encloses. A cell whose volume is zero, or turned over, stops the setup.

Each block carries **two ghost layers** on every face. Nothing about the dimensions has to
be declared:

- the mesh is 3-D if a block has more than one cell in $k$, 2-D if not and a block has more
  than one cell in $j$, and 1-D otherwise. The flux loops over a direction with one cell
  simply do no work;
- a 2-D mesh is an **axisymmetric wedge** when, in every block, the $k = 1$ node plane is the
  $k = 0$ plane turned about the $x$ axis by the same angle, and **plane** when it is the
  $k = 0$ plane translated. Every node is checked to 1 % of the thickness of its block; a
  mesh with blocks of both kinds, or a block that fits neither, stops at setup. The setup
  summary names the type.

The ghost nodes, one layer beyond every face, are extrapolated linearly along the index
lines, $\mathbf{x}_{-1} = 2\mathbf{x}_0 - \mathbf{x}_1$, so a ghost cell has the width of the
first cell across the face, whatever the stretching. Beyond the side faces of a wedge they
are rotated by its angle instead.

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

A side whose three cells include an empty one (a density at or below $10^{-6}$ kg/m³) is
first order: $\mathbf{P}_L$ is then $\mathbf{P}_i$, and $\mathbf{P}_R$ is $\mathbf{P}_{i+1}$.
Each side looks at its own three cells — $i-1$, $i$, $i+1$ for $\mathbf{P}_L$ and $i$, $i+1$,
$i+2$ for $\mathbf{P}_R$ — so the two sides of a face fall back as their mirror images do, and
a symmetric flow stays symmetric to round-off.

A family of a solidifying material reconstructs, in place of $T_p$, its liquid-branch temperature $T_\ell = e/c_l$,
which is continuous where the content nucleates, and takes each face's $T_p$ and $f$ from the face's $T_\ell$ and
$\chi$ by the rule of a cell ([Solidification](physics.md#solidification)); a liquid face is unchanged, and a cold solid,
whose $T_\ell$ is negative, gets its physical temperature back.

### Limiters

`flux-limiter` selects $\phi$ (first order is `space-reconstruction = first-order`, not a limiter). With
$r = y/x$, $\phi(x, y) = \psi(r)\,x$: for $\mathbf{P}_L$, $x = s_i$ is the forward and $y = s_{i-1}$ the backward
slope of cell $i$; for $\mathbf{P}_R$ the same of cell $i+1$.

| `flux-limiter` | $\psi(r)$ | $\nu_{\max}$ |
|---|---|---|
| `minmod` | $\max(0, \min(1, r))$ | 2/3 |
| `vanalbada` | $(r^2 + r)/(r^2 + 1)$ | $2 - \sqrt{2} \approx 0.586$ |
| `vanleer` | $(r + \lvert r\rvert)/(1 + \lvert r\rvert)$ | 1/2 |
| `ospre` | $1.5\,(r^2 + r)/(r^2 + r + 1)$ | 1/2 |
| `umist` | $\max(0, \min(2r, 0.25 + 0.75r, 0.75 + 0.25r, 2))$ | 1/2 |
| `osher` | $\max(0, \min(r, 1.5))$ | 4/7 $\approx 0.571$ |
| `sweby` | $\max(0, \min(1.5r, 1), \min(r, 1.5))$ | 4/7 |
| `mc` | $\max(0, \min(2r, (1 + r)/2, 2))$ | 1/2 |
| `koren` | $\max(0, \min(2r, (1 + 2r)/3, 2))$ | 1/2 |
| `superbee` | $\max(0, \min(2r, 1), \min(r, 2))$ | 1/2 |

$\nu_{\max}$ is the largest `cfl` at which explicit Euler (`time-scheme = euler`) keeps the
reconstruction TVD for 1-D scalar advection on a uniform mesh. It follows from Harten's
sufficient condition in incremental form, over every $r$ and both signs of the speed:
$\nu_{\max} = 1/(1 + A/2 - \psi_{\min}/2)$ for $a > 0$ and $1/(1 + \psi_{\max}/2 - B/2)$ for
$a < 0$, with $A$ and $B$ the largest and smallest $\psi(r)/r$. It bounds `cfl` because `cfl`
bounds the sum of the direction limits ([time step](time-integration.md#the-time-step)).

- Van Albada's $\psi$ is not clipped at $r < 0$ (it goes down to $(1 - \sqrt{2})/2$), so its
  bound is below the 0.624 of $r > 0$ alone.
- OSPRE's goes down to $-1/2$ at $r = -1/2$, which brings its bound from 4/7 to 1/2.
- $\mathbf{P}_R$ reads Osher's $\psi$ at the ratio of cell $i+1$, so $a < 0$ gives 4/7 where
  $a > 0$ gives 2/3.

Above $\nu_{\max}$ the run goes on with a warning,
`[WARNING] [ICE-Numerics] forward Euler with MUSCL ...`, at setup and whenever `input.ini`
is read again; with `irs = true` an `[INFO]` line takes its place. RK2 and RK3 are not
checked.

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

The form is that of Tramel, Nichols and Buning (AIAA 2009-3988, Eqs. 29-30), with one
change. There the weight is $1 - \tanh(10\,\varphi^3)$ with $\varphi = s/\Delta$; with
$\Delta = 20$ this is a switch, 1 below the threshold $1/\Delta$ and 0 above it. ICE takes
$\varphi = s\,\Delta$, which keeps the threshold but spreads the fall of $\beta$ over
$0 \le s < 1/\Delta$, from 1 at $s = 0$ to $1 - \tanh 10$ ($4\times10^{-9}$) at the threshold.

The paper also thickens the sensor (its Eq. 31 takes the largest $s$ of the 27 cells around a
cell). ICE gives each face the weight of the more shocked of its two cells,
$\min(\beta_L, \beta_R)$, which treats a face the same from both sides and keeps mirror-image
flows symmetric to round-off.

The continuous weight matters for steady computations: with a switch, the cells near the
threshold flip from one iteration to the next and the residual levels off instead of
falling. The price is a little more dissipation below the threshold: on the AG crossing
jets, the RMS distance to the exact solution grows from $2.7\times10^{-2}$ to
$2.8\times10^{-2}$ kg/m³.

A domain split into connected blocks gets the same weights as in one block, so the two
differ by round-off only, which the detector can amplify where the crossing jets collide:
$2\times10^{-9}$ of the peak density after 5000 iterations, less than a change of the last
digit of an inlet value ([Multi-block and MPI](../vv/multiblock-mpi.md)). Blocks joined by
chimera faces still differ where the detector acts at the interface.

On the [crossing jets](../vv/crossing-jets.md) the detector lets the MK and AG runs converge, keeps the AG pressure
tensor off its floor in the nearly empty cells around the jets, and keeps AG from gathering mass on either side of the
crossing. It costs accuracy where the flow is smooth: without it MK and AG are nearer their exact solutions, but neither
converges (MK's residual stalls; AG's residual cycles). The crossing-jets tests run without it, for all three
closures.

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

A face between two blocks (`101`, `201`) has a record in each block. Its flux is computed
once, by one of the two blocks — the one on the low-index side when the face is the high
face of one block and the low face of the other, otherwise the first in
(block, face, k, j, i) order — and the other block takes its negative. The two blocks
therefore exchange exactly the same mass, momentum and energy across the cut. With MPI,
the other block's rank sends its cell's detector weight, and gets the flux back, once per
stage.

The state an inlet injects is the state of the boundary face, not of the ghost cell's centre: the reconstructions take the
length of an injecting ghost as zero, so the first cell's slope toward the inlet is the one-sided slope from the face, on
the face itself and on the first interior face, and the scheme stays second order at the inlet. Where the particles leave
through an inlet face, or its mass flux is zero, the ghost is a copy or a mirror of the cell and keeps its length.

A reconstructed state that is unphysical has its slopes halved until it is not. If it is still unphysical at first order —
a ghost that nothing filled, a NaN in the stencil — ICE prints the stencil and stops with a non-zero status.

## Source terms

Drag, convective and radiative heat, and evaporation are evaluated **explicitly**, once per Runge-Kutta stage, from the current primitive state, and added to the residual before it is integrated. There
is no point-implicit or operator-split treatment.
