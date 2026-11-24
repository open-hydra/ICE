# Time Integration

ICE advances the state with an explicit Runge-Kutta scheme. Each cell takes the largest step its flow, its sources
and its mesh allow; a time-accurate run then uses the smallest of them everywhere, a steady run lets each cell keep its
own.

## Runge-Kutta schemes

Each stage is written in the low-storage form

$$
\mathbf{U}^{(k)} = \mathbf{U}^n + c_k\Big(\mathbf{U}^{(k-1)} - \mathbf{U}^n
                    + \Delta t\,\mathcal{R}\big(\mathbf{U}^{(k-1)}\big)\Big),
\qquad \mathbf{U}^{(0)} = \mathbf{U}^n,
$$

with

| Scheme | `time-scheme` | $c_1, c_2, c_3$ |
|---|---|---|
| Forward Euler | `euler` | 1 |
| SSP-RK2 (Heun) | `RK2` | 1, 1/2 |
| SSP-RK3 (Shu-Osher) | `RK3` | 1, 1/4, 2/3 |

Expanding the last row gives the familiar Shu-Osher form
$\mathbf{U}^{(2)} = \tfrac34\mathbf{U}^n + \tfrac14(\mathbf{U}^{(1)} + \Delta t \mathcal R^{(1)})$,
$\mathbf{U}^{n+1} = \tfrac13\mathbf{U}^n + \tfrac23(\mathbf{U}^{(2)} + \Delta t \mathcal R^{(2)})$.
The measured orders of these schemes, on a problem with an exact solution, are in
[Drag and Heat](../vv/verification/exchange.md#a-relaxation-under-stokes-drag).

### One stage

Every stage is a full pass:

1. the sources;
2. the ghost cells, with the MPI exchange;
3. the shock-detector weights;
4. the boundary fluxes, then the interior fluxes;
5. the residual, and its smoothing if switched on;
6. the update.

The update is done in conservative variables: the primitives are converted, updated and converted back.

### Empty and invalid cells

Converting back, ICE floors a density below $10^{-25}$. The density, velocity, temperature, number density and the
pressure (IG) or the diagonal of the tensor (AG) are set to $10^{-25}$; the off-diagonal terms, and the fractions $f$
and $\chi$ of a solidifying family, to 0.

| Density after the update | What ICE does |
|---|---|
| between 0 and $10^{-25}$ | floors it, silently |
| negative | floors it and counts it: the first one is reported with its family, block, cell, stage, level and iteration, then the count at every progress line and at the end of the run (a cell counts once per stage) |
| NaN or infinite (any variable) | stops the run |

## The time step

### The convective limit

For each cell ICE adds up the signal speeds of the directions the mesh resolves (three on a 3-D mesh, two on a 2-D
one, one on a 1-D one):

$$
\frac{1}{\Delta t_c} = \sum_{d=1}^{n_\mathrm{dir}} \Big(|\mathbf{u}\cdot\nabla\xi_d| + a_d\,|\nabla\xi_d|\Big).
$$

- $\nabla\xi_d$ is row $d$ of the cell's metric tensor: normal to the cell's faces of direction $d$, with modulus one
  over the width of the cell across them.
- $a_d$ is the closure's signal speed along that direction:

| Closure | $a_d$ |
|---|---|
| MK | 0: the step is set by convection alone |
| IG | $\sqrt{3P/\rho_p}$ |
| AG | $\sqrt{3P_{dd}/\rho_p}$, with $P_{dd}$ the tensor's component along the direction |

A cell with nothing to carry (an MK cell at rest) takes $10^5$ s.

The rule sums the direction limits rather than taking the smallest, because ICE does not split the scheme by
direction (W. Hundsdorfer, *Numerical solution of advection-diffusion-reaction equations*, CWI lecture notes, 2000,
sec. 9.3). It uses the width of the cell rather than its edge, which is longer on a skewed cell.

### The cell step

$$
\Delta t_i = \min\Big(\texttt{dt-max},\ \ \mathrm{CFL} \cdot \Delta t_c,\ \ \texttt{tau-factor}\cdot\tau_{p,i}\Big).
$$

- `dt-max` is a fixed ceiling. On a coarse mesh or in a slow flow it, rather than the CFL condition, sets the pace.
- A CFL ramp of $N$ iterations scales the CFL number by $\min(1, \text{iteration}/N)$.
- `tau-factor` follows the drops' relaxation time $\tau_p$, cell by cell (below).

### Stable CFL numbers

These are the von Neumann limits of the linearised MUSCL scheme with upwind fluxes (van Leer, van Albada and MC slopes,
in smooth flow):

| Scheme | Largest stable CFL |
|---|---|
| RK2 | 1 |
| RK3 | 1.17 (1.256 with first-order fluxes) |
| Euler | 1 with first-order fluxes; with MUSCL, the limiter's bound |

With Euler and MUSCL, the limited scheme stays TVD only up to the limiter's bound: 0.5 for van Leer on a uniform mesh,
and 1/2 on cells of different lengths ([Limiters](numerics.md#limiters)). Above it ICE warns (with `irs = true`, an
information line instead) and goes on. Keep a margin
below these limits where limiters and shocks act.

### Why `dt-max` and `tau-factor`

The source terms are explicit, and the relaxation time $\tau_p$ does not appear in the CFL condition at all. Without a
further bound, nothing would stop the step from overshooting the relaxation.

The explicit relaxation is stable for:

| Exchange | Euler, RK2 | RK3 |
|---|---|---|
| velocity (every closure) | $\Delta t < 2\tau_p$ | $\Delta t < 2.51\tau_p$ |
| IG and AG dispersion, which relaxes twice as fast | $\Delta t < \tau_p$ | $\Delta t < 1.25\tau_p$ |

The default `tau-factor = 1` respects the first row; a coupled IG or AG run with RK2 wants `tau-factor = 0.5`.

`tau-factor` acts only in coupled runs, in cells that carry drops ($\rho_p > 10^{-6}$) with a finite $\tau_p$, and
`tau-factor = 0` switches it off. $\tau_p$ is the value of the previous step, so the first step of a run is bounded by
`dt-max` alone.

### Directions the mesh does not resolve

A direction with a single cell adds nothing to the step. Across it a cell has boundary faces only, and the only motion
there is a velocity component that is odd across a mirror face, such as $w$ in a 2-D layer of thickness $\Delta z$,
which decays as $dw/dt = -2a\,w/\Delta z$. No explicit scheme can follow that at the in-plane step, so where those faces
mirror the flow ICE sets the component to zero and does not integrate it
([Boundary Conditions](../user/boundary-conditions.md#directions-with-a-single-cell)).

### Time-accurate or steady

| | Time-accurate | Steady |
|---|---|---|
| Step | the smallest $\Delta t_i$ over every cell, block, family and MPI rank, for all of them | each cell its own $\Delta t_i$ |
| Clock | advances by the step | none: the intermediate fields have no physical time |
| Use it when | the answer at a given time matters | only the converged field matters; it converges faster |

With `tau-factor`, a time-accurate step follows the smallest $\tau_p$ in the domain: a family of small or burning drops
slows the whole run.

### Burning and empty cells

The rule above does not see two things that change a cell, so ICE adds two more bounds.

- **A burning cell, in a steady coupled run**, steps at most half its bulk density over its burn rate. The burn per
  unit mass grows without bound as the drops burn down, and nothing else limits the step by it. A cell at or below
  $\rho_p = 10^{-6}$ that still holds drops also takes the `tau-factor` cap.
- **An empty cell that is about to fill** steps no further than what fills it. It fills from births in force (a steady
  run's generated mass, a coupler's births) or from a burning family's own cells emptying into it.
    - Steady: at most the family's smallest step over the cells it occupies and the cells the births fill, and the
      cell's own width at the family's largest speed.
    - Time-accurate: the global step at most the births' step in the empty cells they fill.

A family that nothing fills and that does not burn keeps the steps above, bit for bit.

## Residual smoothing

Implicit residual smoothing filters the residual after it is formed and before the state is updated. It is a Jacobi
solve of

$$
\tilde{\mathcal R}_i - \beta\big(\tilde{\mathcal R}_{i-1} - 2\tilde{\mathcal R}_i + \tilde{\mathcal R}_{i+1}\big) = \mathcal R_i,
$$

direction by direction, with a free coefficient $\beta$ and two Jacobi iterations per direction. The residual's ghost
cells take zero-gradient values, held fixed across the iterations.

Smoothing widens the stability region, so a larger CFL number becomes usable. It also changes the answer, so it belongs
to steady runs, not time-accurate ones.

## Grid sequencing

Grid sequencing gives the fine grid a good initial guess cheaply. It is not a multigrid cycle: there is no transfer
from fine to coarse during the solve, and it does not speed up convergence on the fine grid itself.

1. Coarse grids halve every block dimension, so each block must be divisible by $2^{\texttt{levels}-1}$.
2. At set-up the initial solution is restricted onto every level, so that the coarsest one starts from the initial
   condition. The restriction is volume-weighted, in conserved variables, so it conserves mass.
3. The run starts on the coarsest level. Each level runs its own iterations.
4. The solution is then interpolated (trilinearly) onto the next finer grid, and the iteration count starts again.

- While a coarse level runs, it writes its own files, `field-level<n>` and `residuals-level<n>`, with that level's
  grid, fields, time and iteration count.
- A restart resumes on the level of the solution it reads ([Restarting](../user/using.md#restarting)). The solution
  file is written on the finest level, so a restart from it does not run the coarse levels again.
- Each level is a complete grid with its own metrics, its own [boundary table](../user/numerics.md#grid-levels) and,
  under MPI, its own halo schedule. Block ownership is decided once on the fine grid and kept on every level.
