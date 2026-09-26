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
once at setup from the node coordinates and stored.

Each block carries **two ghost layers** on every face. The dimensionality is inferred
from the first block: a mesh with one cell in $k$ is 2-D, one cell in both $j$ and $k$
is 1-D. Nothing has to be declared — the flux loops over a direction with a single
cell simply do no work. A 2-D mesh whose two node planes make the same angle about the $x$
axis at both ends of its first node line is an axisymmetric wedge of that angle: the ghost
nodes beyond its side faces are then rotated by the angle, not extrapolated in a straight
line.

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
1 - \tanh\!\big(10\,(s/\Delta)^3\big), & s < 1/\Delta \\[2pt]
0, & \text{otherwise}
\end{cases}
\qquad \Delta = 20 .
$$

So $\beta = 1$ in smooth flow, where the scheme is the plain MUSCL one, and falls to 0 at a shock, where it drops to first order.

The weight is close to a switch, and that has two consequences. In a steady computation the
cells near the threshold can keep switching from one iteration to the next, so the residual
levels off instead of converging. And a face on a block boundary takes the weight of its own
block's cell, while inside a block every face takes the weight of the cell on its low-index
side, so a domain split into blocks switches some faces differently from the same domain in
one block. The detector earns its place where a collision would otherwise go wrong: it keeps
the MK delta shock of the [crossing jets](../vv/crossing-jets.md) from creeping upstream and
the AG pressure tensor from reaching its floor where the jets cross. The IG crossing-jets
cases run without it.

## Riemann solvers

The face flux is a two-state flux solved via three possible schemes.

| Name | Form | Applies to |
|---|---|---|
| **Saurel** | Sign of the mean normal velocity selects the donor state; there is no pressure and no sound speed to upwind against | MK only. It assembles the flux from the monokinetic variable layout, so ICE refuses it for IG and AG |
| **Rusanov** | $\tfrac12(\mathbf F_L + \mathbf F_R) - \tfrac12 A\,(\mathbf U_R - \mathbf U_L)$, with $A$ the largest of $|u_n \pm a|$ on the two sides and $a$ the signal speed across the face: $\sqrt{3P/\rho_p}$ for IG, $\sqrt{3P_{nn}/\rho_p}$ for AG | Any closure; the default for IG and AG |
| **HLLE** | Two-wave solver with Roe-averaged speed estimates, falling back to the upwind flux when both waves run the same way | Any closure with a sound speed |

The time step uses the same directional speed, direction by direction.

HLLE is less dissipative than Rusanov on a contact and is worth trying when a contact is being smeared, at the cost of a Roe average per face.

## Boundary faces

Boundary fluxes are built from the ghost values, so every boundary type — connection, chimera, symmetry, extrapolation, inlet —
reaches the flux loop through the same path. The ghost fill is described under [Boundary Conditions](../user/boundary-conditions.md).
On the side faces of an axisymmetric wedge the ghost is the cell's mirror image, so the flux carries no mass, and for IG and
AG the pressure on those faces is the hoop term of the radial momentum; MK has no flux there.

A reconstructed state that is unphysical has its slopes halved until it is not. If it is still unphysical at first order —
a ghost that nothing filled, a NaN in the stencil — ICE prints the stencil and stops with a non-zero status.

## Source terms

Drag, convective and radiative heat, and evaporation are evaluated **explicitly**, once per Runge-Kutta stage, from the current primitive state, and added to the residual before it is integrated. There
is no point-implicit or operator-split treatment.
