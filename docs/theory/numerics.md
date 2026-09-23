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
cell simply do no work.

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

## Riemann solvers

The face flux is a two-state flux solved via three possible schemes.

| Name | Form | Applies to |
|---|---|---|
| **Saurel** | Sign of the mean normal velocity selects the donor state; there is no pressure and no sound speed to upwind against | MK only. It assembles the flux from the monokinetic variable layout, so ICE refuses it for IG and AG |
| **Rusanov** | $\tfrac12(\mathbf F_L + \mathbf F_R) - \tfrac12 A\,(\mathbf U_R - \mathbf U_L)$, with $A$ the largest of $|u_n \pm a|$ on the two sides | Any closure; the default for IG and AG |
| **HLLE** | Two-wave solver with Roe-averaged speed estimates, falling back to the upwind flux when both waves run the same way | Any closure with a sound speed |

HLLE is less dissipative than Rusanov on a contact and is worth trying when a contact is being smeared, at the cost of a Roe average per face.

## Boundary faces

Boundary fluxes are built from the ghost values, so every boundary type — connection, chimera, symmetry, extrapolation, inlet —
reaches the flux loop through the same path. The ghost fill is described under [Boundary Conditions](../user/boundary-conditions.md).

## Source terms

Drag, convective and radiative heat, and evaporation are evaluated **explicitly**, once per Runge-Kutta stage, from the current primitive state, and added to the residual before it is integrated. There
is no point-implicit or operator-split treatment.
