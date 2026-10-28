# Time Integration

## Runge-Kutta schemes

An explicit Runge-Kutta time scheme is used to advance the flow in time. Each stage is written in the low-storage form

$$
\mathbf{U}^{(k)} = \mathbf{U}^n + c_k\Big(\mathbf{U}^{(k-1)} - \mathbf{U}^n
                    + \Delta t\,\mathcal{R}\big(\mathbf{U}^{(k-1)}\big)\Big),
\qquad \mathbf{U}^{(0)} = \mathbf{U}^n,
$$

with

| Scheme | $c_1, c_2, c_3$ |
|---|---|
| Forward Euler | 1 |
| SSP-RK2 (Heun) | 1, 1/2 |
| SSP-RK3 (Shu-Osher) | 1, 1/4, 2/3 |

Expanding the last row recovers the familiar Shu-Osher form
$\mathbf{U}^{(2)} = \tfrac34\mathbf{U}^n + \tfrac14(\mathbf{U}^{(1)} + \Delta t \mathcal R^{(1)})$,
$\mathbf{U}^{n+1} = \tfrac13\mathbf{U}^n + \tfrac23(\mathbf{U}^{(2)} + \Delta t \mathcal R^{(2)})$.

Every stage is a full pass: sources, ghost fill (including the MPI halo exchange), shock-detector weights, boundary fluxes, interior fluxes, residual, optional smoothing, update. The stage update is done in conservative variables — the primitives are converted in, integrated, and converted back — and the converted-back state is checked for positivity.

The measured orders of these schemes on a problem with an exact solution are in [Code Verification](../vv/verification.md).

## The time step

For each cell ICE adds up the signal speeds of the directions the mesh resolves — all three on a
3-D mesh, $i$ and $j$ on a 2-D one, $i$ on a 1-D one:

$$
\frac{1}{\Delta t_c} = \sum_{d=1}^{n_\mathrm{dir}} \Big(|\mathbf{u}\cdot\nabla\xi_d| + a_d\,|\nabla\xi_d|\Big),
$$

with $\nabla\xi_d$ the row $d$ of the cell's metric tensor, normal to its $d$ faces and of modulus one over the width of the cell across them, and $a_d$ the closure's signal speed along $\hat{\mathbf e}_d = \nabla\xi_d/|\nabla\xi_d|$ — zero for MK, so the MK step is set by convection alone; $\sqrt{3P/\rho_p}$ for IG; $\sqrt{3P_{dd}/\rho_p}$ for AG, with $P_{dd} = \hat{\mathbf e}_d\cdot\mathsf P\,\hat{\mathbf e}_d$. A cell with nothing to carry, an MK cell at rest, takes $10^5$ s. The cell step is then

$$
\Delta t_i = \min\Big(\texttt{dt-max},\ \ \mathrm{CFL} \cdot \Delta t_c,\ \ \texttt{tau-factor}\cdot\tau_{p,i}\Big),
$$

so `dt-max` bounds the step that is actually taken. If a CFL ramp is set, the CFL factor is additionally scaled by $\min(1, \text{iteration}/N)$, ramping it linearly from zero over the first $N$ iterations; the ceiling is applied after that too.

A scheme that is not split by direction needs the sum of the direction limits, not the smallest of
them (W. Hundsdorfer, *Numerical solution of advection-diffusion-reaction equations*, lecture
notes, CWI, 2000, sec. 9.3), and the width of the cell rather than its edge, which is longer on a
skewed cell. For the linearised MUSCL scheme with upwind fluxes (the van Leer, van Albada and MC
slopes, in smooth flow) the rule is stable for $\mathrm{CFL} \le 1$ with RK2 (Heun) and
$\mathrm{CFL} \le 1.17$ with RK3, $1.256$ with first-order fluxes (von Neumann limits of the
linearised scheme); explicit Euler is stable with first-order fluxes only, up to 1. Keep a margin
below these where limiters and shocks act. The rule used to take the smallest of
$\ell_d/(|\mathbf u\cdot\hat{\mathbf e}_d| + a_d)$ over the three index directions, with $\ell_d$
the length of the cell edge: the same CFL now gives a step down to one half of that on square 2-D
cells and to one third on cubes.

The directions the mesh does not resolve add nothing. Across them a cell has boundary faces only,
and the flux between the cell and its ghost gives a mode only to the components that are odd
across a mirror face: $w$ between the $k$ faces of a 2-D layer of thickness $\Delta z$,
$dw/dt = -2a\,w/\Delta z$, which on a layer thin against its in-plane cells no explicit scheme
follows at the in-plane step. Where those faces mirror the flow, ICE sets those components to zero
and does not integrate them
(see [Boundary Conditions](../user/boundary-conditions.md#faces-of-a-direction-with-one-cell)). The
old rule bounded the step by that width too, so cells thin across those faces, as near the axis of
a wedge, now take a larger step.

`dt-max` and `tau-factor` exist because the source terms are explicit: the relaxation time $\tau_p$ does not appear in the CFL condition at all, so nothing else stops the step from overshooting it. `dt-max` is a fixed ceiling with a finite default, and on a coarse mesh or a slow flow it — rather than the CFL condition — is what sets the pace. `tau-factor` follows $\tau_p$ cell by cell. It applies to coupled runs only, and only in cells that carry the phase ($\rho_p > 10^{-6}$) and have a finite $\tau_{p,i}$; `tau-factor = 0` switches it off. $\tau_{p,i}$ is the value the source terms computed during the previous step, so the first step of a run is bounded by `dt-max` alone.

The explicit relaxation of the velocity is stable for $\Delta t < 2\tau_p$ with Euler and RK2 and $2.51\tau_p$ with RK3, which the default `tau-factor = 1` respects. The IG and AG dispersion relaxes twice as fast ($\Delta t < \tau_p$ with Euler and RK2, $1.25\tau_p$ with RK3), so a coupled IG or AG run with RK2 wants `tau-factor = 0.5`.

### Time-accurate mode

The smallest $\Delta t_i$ over every cell, every block, every
family and every MPI rank becomes the step for all of them, and the simulation clock advances by it. This is the mode to use whenever the answer at a given physical time matters. With `tau-factor` the global step follows the smallest $\tau_p$ in the domain, so a family of small or burning droplets slows the whole run.

### Steady-state mode

Every cell keeps its own $\Delta t_i$. Convergence to a steady
state is faster, since a cell is no longer held back by the most restrictive cell in the domain, but the intermediate fields have no physical time.

## Residual smoothing

This technique applies an implicit residual smoothing pass after the residual is formed and before the state is updated. It is a Jacobi solve of

$$
\tilde{\mathcal R}_i - \beta\big(\tilde{\mathcal R}_{i-1} - 2\tilde{\mathcal R}_i + \tilde{\mathcal R}_{i+1}\big) = \mathcal R_i,
$$

applied direction by direction, with $\beta$ a free coefficient and **two** Jacobi iterations per direction. The residual ghost cells are filled by zero-gradient extrapolation and held fixed across the iterations. Smoothing widens the stability region, so a larger CFL number becomes usable; it also changes the answer, and so belongs to steady-state runs rather than time-accurate ones.

## Grid sequencing

Coarse grids are built by halving every block dimension, so
each block must be divisible by $2^{\texttt{levels}-1}$;.

The run then **starts on the coarsest level** and works up. Each level is advanced for its own iterations, after which the solution is prolongated onto the next finer grid and the iteration counter restarts. There is no fine-to-coarse transfer during the solve: this is grid sequencing — a cheap way to get a good initial guess on the fine grid — not a multigrid cycle, and it does not accelerate the fine-level convergence once it is reached. A restart resumes on the level of the solution it reads (`ICE_GRID_LEVEL`, see [Restarting](../user/using.md#restarting)): the solution is written on the finest level, so a restart from it does not run the coarse levels again.

The one fine-to-coarse transfer happens at setup, where the initial solution is restricted onto every level in turn so that the coarsest one starts from the initial condition rather than from nothing. Restriction is volume-weighted and done in conserved variables, so it conserves mass; prolongation is trilinear.

A level is a complete grid and carries everything one needs: its own metrics, its own [boundary table](../user/input.md#grid-levels) and, under MPI, its own halo schedule. Block ownership is the exception — it is decided once on the fine grid and every level keeps it, since all levels have the same blocks in the same order.

