# Drag and Heat

A cloud with no gradients in space: the transport operator contributes nothing, and each
source term reduces to an ODE with a known solution. A and B check the drag and the heat
exchange against a closed form, and J checks the time step near the relaxation time. D and E
run every drag and Nusselt correlation through the same box.

| Case | ctest | What it pins down | Reference |
|---|---|---|---|
| [A](#a-relaxation-under-stokes-drag) | `DragRelaxation` | The momentum source, and the order of the time scheme | Closed form |
| [J](#j-the-step-at-the-relaxation-time) | `TauLimit` | `tau-factor`: the step bounded by $\tau_p$ in both stepping modes; cells below the empty-cell density left out | Closed form |
| [B](#b-thermal-relaxation) | `ThermalRelaxation` | The energy source at $Nu = 2$ | Closed form |
| [D](#d-every-drag-correlation) | `DragMatrix` | All twelve drag laws, up to $Ma = 2.6$ | Independent RK4 |
| [E](#e-every-nusselt-correlation) | `HeatMatrix` | All seven heat laws, velocity and temperature relaxing together | Independent RK4 and closed form |

<figure>
  --8<-- "docs/vv/images/verification-relaxation.svg"
</figure>

## A. Relaxation under Stokes drag

A cloud at rest in a uniform gas stream, one-way coupled. Nothing varies in space, so the
momentum equation reduces to an ODE with an exact solution:

$$
u_p(t) = u_g + (u_{p,0} - u_g)\,e^{-t/\tau_p},
\qquad
\tau_p = \frac{\rho_{al} d_p^2}{18 \mu_g}.
$$

That $\tau_p$ is not assumed: with $C_d = 24/Re$, ICE's relaxation time
$8\rho_{al}R_p/(3\rho_g C_d |\Delta u|)$ reduces to it exactly, and the test confirms the
resulting curve. With $\rho_{al} = 1000$ kg/m³, $d_p = 10^{-4}$ m and
$\mu_g = 1.8\times10^{-5}$ Pa s, $\tau_p = 0.030864$ s.

| $t/\tau_p$ | $u_p/u_g$ (ICE) | $u_p/u_g$ (exact) | Relative error |
|---|---|---|---|
| 0.1 | 0.095560464 | 0.095560623 | $1.6\times10^{-7}$ |
| 0.5 | 0.394801708 | 0.394802241 | $5.3\times10^{-7}$ |
| 1.0 | 0.632546406 | 0.632547052 | $6.5\times10^{-7}$ |
| 2.0 | 0.864977857 | 0.864978331 | $4.7\times10^{-7}$ |
| 3.0 | 0.950224617 | 0.950224879 | $2.6\times10^{-7}$ |
| 5.0 | 0.993279221 | 0.993279280 | $5.9\times10^{-8}$ |

**What must not move.** The density, the number density and the particle temperature stay at
their initial values, and the field stays uniform to round-off. The temperature is a real
check: the drag work in the energy equation cancels exactly the kinetic energy it produces, so
an imbalance there would show up as spurious heating.

**Time-step refinement.** The field stays uniform, so the error is purely the time
integration of the source term. Halving the CFL number halves the step, and the observed
order is the scheme's:

| CFL | 0.8 | 0.4 | 0.2 | Observed order |
|---|---|---|---|---|
| RK2 | $2.47\times10^{-7}$ | $6.17\times10^{-8}$ | $1.54\times10^{-8}$ | 2.00 |
| RK3 | $1.32\times10^{-10}$ | $1.64\times10^{-11}$ | $2.05\times10^{-12}$ | 3.00 |

- The refinement runs on two cells in $y$, with a velocity along $y$ shared by the gas and the
  cloud, which the CFL limit then reads. A direction of a single cell sets no limit, so the run
  of the first table, on one cell in $y$, steps at `dt-max`.
- The case also checks the `dt-max` ceiling: a run with a ceiling well below its CFL limit
  must end at an exact multiple of it, since every step is `dt-max` itself, not
  `cfl * dt-max`.

## J. The step at the relaxation time

**Setup.** Case A's cloud with particles small enough that the relaxation time is below
`dt-max`: $d_p = 2$ µm, $\rho_{al} = 2000$ kg/m³, $\mu_g = 2\times10^{-5}$ Pa s, so
$\tau_p = \rho_{al}d_p^2/(18\mu_g) = 2.222\times10^{-5}$ s $= \mathrm{dt\text{-}max}/4.5$.

**Why it needs a limit.** A cloud at rest has no signal speed under MK, so the CFL condition
does not bound its step. With `dt-max` as the only bound, SSP-RK3 integrates the relaxation at
$z = -4.5$, outside its stability interval: $R(z) = 1 + z + z^2/2 + z^3/6 = -8.56$, and the
slip grows by that factor every step until the convective limit takes over. At $t = 20\tau_p$
the cloud is still $1.8\times10^3$ m/s off the gas velocity.

**What `tau-factor` does.** It bounds the step of every populated cell by
$\mathrm{tau\text{-}factor}\times\tau_p$ (default 1). Then $R(-1) = 1/3$ and the slip falls by
3 per step. The first step is the exception: $\tau_p$ is not known before the first source
evaluation, so that step is `dt-max` alone and overshoots to $u_p = 95.6$ m/s. Every step after
it is exactly $\tau_p$, and the end time shows it: $(t - \mathrm{dt\text{-}max})/\tau_p$ is an
integer.

| Leg | Setting | Check | Tolerance | Measured |
|---|---|---|---|---|
| J1 | time accurate, to $20\tau_p$ | $\max\lvert u_p/u_g - 1\rvert$ | $10^{-3}$ | $1.99\times10^{-7}$ ($= 85.6\cdot3^{-16}/u_g$) |
| | | $(t - \mathrm{dt\text{-}max})/\tau_p$ an integer | $10^{-6}$ | 16.000000 |
| J1a | `tau-factor` 0.2, `dt-max` $5\times10^{-6}$, to $2\tau_p$ | $\lvert u_p - u_g(1 - e^{-t/\tau_p})\rvert/u_g$ at the time reached | $5\times10^{-3}$ | $1.1\times10^{-4}$ |
| | | $(t - 5\times10^{-6})/(0.2\tau_p)$ an integer | $10^{-6}$ | 9.000000 |
| J2 | two clouds, $d_p$ = 2 and 20 µm, time accurate | small: $\max\lvert u_p/u_g - 1\rvert$ | $10^{-3}$ | $1.99\times10^{-7}$ |
| | | large: $\lvert u_p - u_g(1 - e^{-t/\tau_L})\rvert/u_g$ | $10^{-4}$ | $1.5\times10^{-7}$ |
| | | $(t - \mathrm{dt\text{-}max})/\tau_S$ an integer | $10^{-6}$ | 16.000000 |
| J2l | the same, local steps, 20 iterations | small: $\max\lvert u_p/u_g - 1\rvert$ | $10^{-6}$ | $7.4\times10^{-9}$ ($= 85.6\cdot3^{-19}/u_g$) |
| | | large: on its exponential at $20\,\mathrm{dt\text{-}max}$ | $10^{-4}$ | $1.4\times10^{-6}$ |
| J3 | dilute half, $\rho_p = 10^{-7}$, $\tau_p = 0.611\tau_S$ | $(t - \mathrm{dt\text{-}max})/\tau_S$ an integer | $10^{-6}$ | 16.000000 |
| | | dilute $\rho_p$ and $n$ unchanged | bit for bit | unchanged |

What each leg shows:

- **J1** is the limit itself.
- **J1a** shows that `tau-factor` buys accuracy as well as stability: at $0.2\tau_p$ the RK3
  error per step is $8\times10^{-5}$ of the slip. The comparison is made at the time the run
  reached ($2.025\tau_p$), not at the requested end time.
- **J2** puts two clouds side by side, one per row of a two-row slab, so that no particle
  crosses from one to the other. (Split in $x$, the small particles would flow into the large
  ones' cells and change their number density, hence their $\tau_p$.)
    - In time-accurate mode the global step is the smallest limit in the domain, $\tau_S$. The
      large particles ($\tau_L = 100\,\tau_S$) simply take small steps, $\Delta t/\tau_L = 0.01$,
      which RK3 integrates to $10^{-7}$.
    - In local mode every cell keeps its own limit. The small particles take $\tau_S$-sized
      steps; the large ones keep `dt-max`, since $\tau_L > \mathrm{dt\text{-}max}$, so after 20
      iterations they sit on their exponential at $t = 20\,\mathrm{dt\text{-}max}$, not at
      $20\tau_S$. That check is what tells a per-cell limit from a global one.
- **J3** checks the population floor. A cell whose bulk density is below $10^{-6}$ counts as
  empty, and its $\tau_p$, the ratio of two state variables that mean nothing there, must not
  set the step.
    - The dilute half has $\rho_p = 10^{-7}$ and $n = 2.5\times10^7$: a radius of 0.78 µm and a
      $\tau_p$ of $0.611\tau_S$. Were it counted, the steps would be $0.611\tau_S$ and the end
      time would read 15.883 (26 such steps), which is what a build without the floor gives.
    - That $\tau_p$ is kept close to $\tau_S$ on purpose. A cell whose $\tau_p$ is orders of
      magnitude below the step is unstable under the explicit source whatever the step, so it
      cannot show which cells set the step. At $\Delta t/\tau_p = 1.64$ the dilute half relaxes
      stably and keeps $\rho_p$ and $n$ bit for bit.

## B. Thermal relaxation

The same idea for the energy equation. The particles travel with the gas, so there is no
slip, no drag force and no drag work. With `heat-transfer = Stokes` the Nusselt number is
exactly 2, and the convective exchange is again a linear relaxation:

$$
T_p(t) = T_g + (T_{p,0}-T_g)\,e^{-t/\tau_T},
\qquad
\tau_T = \frac{\rho_{al} c_s d_p^2}{12 k_g} = 0.028846\ \text{s}.
$$

Radiation, the only other term in the energy source, is switched off (`emiss = 0`).

| $t/\tau_T$ | $T_p$ [K] (ICE) | $T_p$ [K] (exact) | Relative error |
|---|---|---|---|
| 0.1 | 390.247598 | 390.247586 | $1.2\times10^{-7}$ |
| 1.0 | 336.744852 | 336.744805 | $4.7\times10^{-7}$ |
| 3.0 | 304.975010 | 304.974991 | $1.9\times10^{-7}$ |
| 5.0 | 300.673583 | 300.673579 | $4.3\times10^{-8}$ |

The errors are relative to the initial 100 K gap. The velocity, density and number density
must not move, and do not.

## D. Every drag correlation

ICE offers twelve drag laws. Case A verifies one of them against a closed form; this case
covers the rest, in three sweeps.

<figure>
  --8<-- "docs/vv/images/verification-correlations.svg"
</figure>

Every one of the nineteen correlations in D and E sits about three decades inside the
tolerance the cases assert.

**At vanishing slip** ($Re = 6.7\times10^{-3}$), the laws that are a Stokes law plus a
correction must reproduce the exact Stokes exponential. This check does not depend on how the
correlation is written in ICE:

| Law | Deviation from Stokes |
|---|---|
| Stokes | $4.1\times10^{-7}$ |
| Morsi-Alexander | $4.1\times10^{-7}$ |
| Schlichting | $2.9\times10^{-4}$ |
| Schiller-Naumann, Wen-Yu, Clift-Gauvin | $1.3\times10^{-3}$ |
| Putnam | $1.6\times10^{-3}$ |

Morsi-Alexander is piecewise, and its $Re \le 0.1$ branch *is* the Stokes law, which is why it
matches to round-off. The others differ by their own correction terms, as they should.

**At finite slip**, every law is compared against an RK4 integration of
$du/dt = (u_g-u)/\tau(u)$ with the same correlation, evaluated independently in Python. Two
sweeps run, at $Re = 67$ and at $Re = 1333$; the second reaches the branches that switch
formula above $Re = 1000$. All twelve agree with the reference to between $4\times10^{-10}$
and $4.4\times10^{-6}$ of $u_g$.

**At transonic and supersonic slip** ($Ma = 1.30$ and $Ma = 2.59$), the four laws with a
compressible branch (Carlson-Hoglund, Henderson, Crowe and Hermsen) are compared the same way,
and agree to between $1.3\times10^{-8}$ and $1.4\times10^{-6}$ of $u_g$. In this range
Henderson leaves its subsonic fit, crosses the linear bridge and reaches its supersonic fit, so
the case also checks that the branches meet: the coefficient of the bridge is fixed by
continuity at $Ma = 1.75$, and both joins are continuous to $2\times10^{-9}$.

!!! note "What this check can and cannot catch"
    It catches a law that is mis-wired, mis-selected, divergent, or integrated incorrectly by
    the solver. It cannot catch a correlation written wrongly in the same way in ICE and in the
    reference: only the low-slip limits above, and a reading against the original
    publications, can do that.

## E. Every Nusselt correlation

The counterpart of D for the convective heat exchange, over ICE's seven laws. The cloud starts
hot *and* slipping, so velocity and temperature relax together: $Re$ and $Ma$ change while the
particles cool, and the Nusselt number follows them. The reference integrates the coupled pair
with RK4.

| Leg | What it checks | Measured |
|---|---|---|
| Vanishing slip | The four laws that tend to $Nu = 2$ reproduce case B's exact relaxation | $7\times10^{-7}$ (Stokes) to $6\times10^{-3}$ (Ranz-Marshall) |
| $Re = 67$ | All seven against the RK4 reference | within $2.6\times10^{-6}$ of the initial gap |
| Constant slip | `JAXA4` against an exact exponential | $7\times10^{-6}$ of the gap at $t = 0.01$ s (tolerance $3\times10^{-4}$) |

- **Vanishing slip.** Ranz-Marshall's error is the largest because its $Re^{1/2}$ correction
  is the largest of the four at this slip. Three laws are left out of this leg and covered by
  the finite-slip one: JAXA1 tends to zero rather than 2, and the two laws with a Mach
  correction keep a finite $Ma/Re$ ratio there.
- **Constant slip.** Under `NoDrag` the particles stay at rest in the moving gas, so
  $Re = 67$, $Ma = 0.029$ and $Nu$ are frozen and the cooling is an exact exponential. The
  measured error is the RK2 error of the default step, and the tolerance is about 40 times
  that. This is the leg with margin on the law's constant: a 1.4 % change of it moves the
  result by $3.4\times10^{-3}$ of the gap, against $1.4\times10^{-3}$ in the finite-slip leg.
