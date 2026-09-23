# Time Integration

## Runge-Kutta schemes

`time-scheme` in `[ICE-Numerics]` selects the integrator. Each stage is written in the
low-storage form

$$
\mathbf{U}^{(k)} = \mathbf{U}^n + c_k\Big(\mathbf{U}^{(k-1)} - \mathbf{U}^n
                    + \Delta t\,\mathcal{R}\big(\mathbf{U}^{(k-1)}\big)\Big),
\qquad \mathbf{U}^{(0)} = \mathbf{U}^n,
$$

with

| `time-scheme` | Scheme | $c_1, c_2, c_3$ |
|---|---|---|
| `euler` | forward Euler | 1 |
| `RK2` | SSP-RK2 (Heun) | 1, 1/2 |
| `RK3` | SSP-RK3 (Shu-Osher) | 1, 1/4, 2/3 |

Expanding the last row recovers the familiar Shu-Osher form
$\mathbf{U}^{(2)} = \tfrac34\mathbf{U}^n + \tfrac14(\mathbf{U}^{(1)} + \Delta t \mathcal R^{(1)})$,
$\mathbf{U}^{n+1} = \tfrac13\mathbf{U}^n + \tfrac23(\mathbf{U}^{(2)} + \Delta t \mathcal R^{(2)})$.

Every stage is a full pass: sources, ghost fill (including the MPI halo exchange),
boundary fluxes, interior fluxes, residual, optional smoothing, update. The stage
update is done in conservative variables — the primitives are converted in, integrated,
and converted back — and the converted-back state is checked for positivity and for
NaNs, which stops the run with the offending cell index if it fails.

The measured orders of these schemes on a problem with an exact solution are in
[Code Verification](../vv/verification.md).

## The time step

For each cell and each direction $d$ ICE forms

$$
\Delta t_d = \frac{\ell_d}{|\mathbf{u}\cdot\hat{\mathbf e}_d| + a},
$$

with $\ell_d$ the cell length along $d$, $\hat{\mathbf e}_d$ the unit vector of the
metric row and $a$ the closure's sound speed — zero for MK, so the MK step is set by
convection alone. The cell step is then

$$
\Delta t_i = \min\Big(\texttt{dt-max},\ \ \mathrm{CFL} \cdot \min_d \Delta t_d\Big),
$$

so `dt-max` bounds the step that is actually taken. If `cfl-rise-threshold` is set to
$N > 0$, the CFL factor is additionally scaled by $\min(1, \text{iteration}/N)$, ramping
it linearly from zero over the first $N$ iterations; the ceiling is applied after that
too.

`dt-max` exists because the source terms are explicit: the drag relaxation time
$\tau_p$ does not appear in the CFL condition at all, so nothing else stops the step
from overshooting it. It has a finite default for that reason, and on a coarse mesh or
a slow flow it — rather than the CFL condition — is what sets the pace. Compare the
reported `Delta t` against `dt-max` to see which is binding.

### Time-accurate mode

`time-accurate = true`: the smallest $\Delta t_i$ over every cell, every block, every
family and every MPI rank becomes the step for all of them, and the simulation clock
advances by it. This is the mode to use whenever the answer at a given physical time
matters, and the one the `time-threshold` stop condition refers to.

### Steady-state mode

`time-accurate = false`: every cell keeps its own $\Delta t_i$. Convergence to a steady
state is faster, since a cell is no longer held back by the most restrictive cell in the
domain, but the intermediate fields have no physical time. In this mode the solution
files record the iteration number in place of the time, and the shell prints the density
residual rather than the step.

## Residual smoothing

`irs = true` applies an implicit residual smoothing pass after the residual is formed
and before the state is updated. It is a Jacobi solve of

$$
\tilde{\mathcal R}_i - \beta\big(\tilde{\mathcal R}_{i-1} - 2\tilde{\mathcal R}_i + \tilde{\mathcal R}_{i+1}\big) = \mathcal R_i,
$$

applied direction by direction, with $\beta$ from `irs-beta` and **two** Jacobi
iterations per direction. The residual ghost cells are filled by zero-gradient
extrapolation and held fixed across the iterations. Smoothing widens the stability
region, so a larger CFL number becomes usable; it also changes the answer, and so
belongs to steady-state runs rather than time-accurate ones.

## Grid sequencing

`levels` in `[ICE-Multigrid]` builds coarse grids by halving every block dimension, so
each block must be divisible by $2^{\texttt{levels}-1}$; ICE checks this at setup and
stops if it is not.

The run then **starts on the coarsest level** and works up. Each level is advanced for
its own `levelN-iter` iterations, after which the solution is prolongated onto the next
finer grid and the iteration counter restarts. There is no fine-to-coarse transfer
during the solve: this is grid sequencing — a cheap way to get a good initial guess on
the fine grid — not a multigrid cycle, and it does not accelerate the fine-level
convergence once it is reached. Output from a coarse level is written to
`part-field-levelN`.

## Stopping

The run ends when any of these holds:

| Condition | Parameter |
|---|---|
| Iteration count on the current level reached | `iter-threshold`, or `levelN-iter` under grid sequencing |
| Simulation time reached | `time-threshold` (time-accurate mode) |
| Density residual fell below the threshold, on the finest level only | `res-threshold` |

The residual is the $L_2$ norm, over all cells and families, of the change in the
primitive variables across the step — $\rho_p$, the three velocity components and
$T_p$ — and it is the $\rho_p$ component that `res-threshold` tests.

!!! note "`res-threshold = 0` disables the residual test"
    A transient in which the density happens not to move — a cloud relaxing only in
    velocity or temperature — has an exactly zero density residual from the first step
    and would otherwise be declared converged immediately. Setting `res-threshold = 0`
    turns the test off and leaves the iteration and time limits in charge.
