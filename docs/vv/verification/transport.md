# Transport

The cases on the other pages hold the cloud uniform, so the transport operator has nothing to
do. These four give it work. C is transport alone; G adds a source term that is uniform in
space; H one that is not; I does it in two dimensions. All four compare against a closed form,
and none can be passed by a solver that transports correctly but couples its source terms
wrongly.

| Case | ctest | What it pins down | Reference |
|---|---|---|---|
| [C](#c-sinusoidal-advection-on-a-periodic-mesh) | `SineAdvection` | Transport alone, and the order of the space scheme | Closed form |
| [G](#g-a-cloud-released-into-a-uniform-gas) | `CloudTranslation` | Transport and drag together | Closed form |
| [H](#h-a-cloud-in-a-straining-gas) | `LinearStrain` | A non-trivial particle velocity field, and the small-Stokes limit | Closed form |
| [I](#i-a-cloud-in-a-prescribed-vortex) | `VortexCloud` | Two dimensions, over four decades of Stokes number | Closed form |

The mesh refinements of all four are collected in [Convergence](#convergence).

## C. Sinusoidal advection on a periodic mesh

**Setup.** No gas at all (0-way coupling). The cloud moves at a uniform $U = 1$ m/s, and its
density carries a sine on a periodic unit domain, so free transport gives

$$
\rho_p(x,t) = \rho_0 + A \sin\!\big[k(x - Ut)\big], \qquad k = 2\pi/L.
$$

**What is measured.**

- The comparison is against the exact *cell average*, so the measured order is not polluted
  by the difference between a cell average and a point value.
- The `dt-max` ceiling holds $\Delta t$ at the same value on all four meshes, so the time
  error is common to them, and the refinement measures the space discretisation.
- Two reconstructions are refined side by side, and both reach their design order (see
  [Convergence](#convergence)).
- Mass is conserved to round-off across the periodic interface, and the velocity stays
  uniform: a pressureless cloud moving at one speed has nothing that could change it.

## G. A cloud released into a uniform gas

Case A leaves the transport operator idle, and case C leaves the source terms idle. This case
is the product of the two: the smallest problem in which getting either one wrong, or coupling
them wrongly, changes the answer.

**Setup.** A frozen gas moves at $U_g = 10$ m/s, and a localized cloud starts at rest. The
Stokes relaxation time $\tau = \rho_{al} d_p^2 / 18\mu$ depends only on the diameter, which
stays uniform because $n$ is seeded proportional to $\rho_p$ and the two obey the same
transport equation.

**The exact solution.** Every particle feels the same drag, so the particle velocity stays
uniform in space for all time, and the cloud translates rigidly:

$$
u_p(t) = U_g\big(1 - e^{-t/\tau}\big), \qquad
X(t) = U_g t - U_g\tau\big(1 - e^{-t/\tau}\big), \qquad
\rho_p(x,t) = \rho_p(x - X(t), 0).
$$

At $t = 3\tau$ the cloud has moved 0.633 m, where the gas would have carried it 0.926 m.

<figure>
  --8<-- "docs/vv/images/verification-translation.svg"
</figure>

Grey is the exact profile, red is ICE, and the labels mark the elapsed time in units of
$\tau$. The cloud keeps its shape: only the top hat smears, and only because a discontinuity
must. What sets its position is the velocity curve on the right, case A's exponential
recovered from a cloud that is now moving.

**The uniformity of $u_p$ is the sharper check.** It is not exact at finite $\Delta t$: the
drag source is evaluated on the density that the flux divergence is about to change, which
leaves a spread wherever the density gradient is steep. That spread is mesh-independent,
survives a first-order reconstruction, and falls at second order in $\Delta t$, which shows
that MK's momentum flux and its drag source are consistent:

| `dt-max` | $2\times10^{-5}$ | $1\times10^{-5}$ | $5\times10^{-6}$ | Observed order |
|---|---|---|---|---|
| $(\max u_p - \min u_p)/U_g$ | $9.22\times10^{-6}$ | $2.32\times10^{-6}$ | $5.81\times10^{-7}$ | 1.99 – 2.00 |

The other checks:

- the mass-weighted $u_p$ matches the exponential to $3\times10^{-8}$ m/s throughout;
- a top hat runs alongside the Gaussian, for the invariants that survive any amount of
  smearing: mass is conserved to $2\times10^{-16}$ across the periodic interface, and the
  centroid sits within 0.007 cells of $X(t)$.

## H. A cloud in a straining gas

In case G the gas has one velocity, so the particle velocity stays uniform and the cloud only
translates. Here the frozen gas carries a uniform strain, $u_g(x) = U_0 + ax$ with
$U_0 = 2$ m/s and $a = 6\ \mathrm{s^{-1}}$, and the cloud can no longer move rigidly: the
particles lag the gas acceleration, the defining behaviour of a dispersed phase. It is the
first case with a non-trivial particle velocity *field*, so the first that can detect a
momentum flux that is wrong in a way a uniform velocity hides.

**The exact solution.** With $\xi = x + U_0/a$ the gas velocity is $a\xi$, and Stokes drag
with a uniform diameter gives one linear equation along a trajectory:

$$
\tau\ddot\xi + \dot\xi - a\xi = 0, \qquad
\lambda_\pm = \frac{-1 \pm \sqrt{1 + 4a\tau}}{2\tau}.
$$

Seeding the cloud at the local gas velocity, $\dot\xi(0) = a\xi(0)$, makes the solution
separable: every particle is displaced by the same factor, $\xi(t) = C(t)\,\xi_0$, with
$C = Ae^{\lambda_+t} + Be^{\lambda_-t}$, $A + B = 1$ and $A\lambda_+ + B\lambda_- = a$. The map
$x_0 \mapsto x(t)$ is affine, which gives the whole Eulerian solution in closed form:

$$
\rho_p(x,t) = \frac{1}{C}\,\rho_{p0}\!\left(\frac{\xi}{C}\right), \qquad
u_p(x,t) = \frac{\dot C}{C}\,\xi ,
$$

exact at every Stokes number, with no crossing at any time since $C > 0$. At
$St = a\tau = 0.185$ and $t = 0.1$ s the cloud is stretched by $C = 1.7096$.

<figure>
  --8<-- "docs/vv/images/verification-strain.svg"
</figure>

Grey is exact, red is ICE. The cloud is stretched and diluted by the same factor, and the
particles run slower than the gas everywhere (the shaded band is the slip). That slip
approaches the textbook asymptote from below as the Stokes number falls.

**Results.** On the finest mesh the centroid lands within 0.004 cells of $C\xi_0$, and the
mass-weighted slope of the velocity field matches $\dot C/C$ to $2\times10^{-6}$ relative.

**The small-Stokes limit.** For a small Stokes number the standard result is
$u_p \approx u_g - \tau\,Du_g/Dt$, here a slip of $a^2\tau\xi$. Expanding the exact slip
coefficient $a - \dot C/C$ in $St$ gives $a^2\tau(1 - 2St + O(St^2))$, so the departure from
the asymptote must be proportional to $St$, with coefficient 2:

| $St = a\tau$ | 0.200 | 0.100 | 0.050 | 0.025 |
|---|---|---|---|---|
| measured slip / $a^2\tau$ | 0.7167 | 0.8385 | 0.9109 | 0.9529 |
| $(1 - \text{ratio})/St$ | 1.42 | 1.62 | 1.78 | 1.89 |

The ratio approaches 1 monotonically and the last row climbs towards 2: the solver has the
right asymptotic behaviour, not merely a plausible one. Against the exact $\dot C/C$, which
holds at any Stokes number, the velocity field agrees to $10^{-5}$ relative at each of these
Stokes numbers.

## I. A cloud in a prescribed vortex

The frozen gas is in solid-body rotation, $u_g = -\Omega y$, $v_g = \Omega x$, and a compact
cloud is released into it at the local gas velocity. It is the only code-verification case
that exercises the $j$-direction fluxes: G and H leave them idle, so a momentum flux wrong
across directions, or a metric wrong on the $j$ faces, cannot show in either.

**The exact solution.** In the complex plane $z = x + iy$ the carrier field is
$u_g = i\Omega z$, so the trajectory equation is case H's over the complex field:

$$
\tau\ddot z + \dot z = i\Omega z, \qquad
\lambda_\pm = \frac{-1 \pm \sqrt{1 + 4i\Omega\tau}}{2\tau}.
$$

Seeding $\dot z(0) = i\Omega z(0)$ makes it separable again: $z(t) = C(t)\,z_0$, with
$C = Ae^{\lambda_+t} + Be^{\lambda_-t}$, $A + B = 1$, $A\lambda_+ + B\lambda_- = i\Omega$. The
map is a rotation composed with a dilation, conformal and non-singular for all time, so the
Eulerian solution is in closed form:

$$
\rho_p(z,t) = \frac{1}{|C|^2}\,\rho_{p0}\!\left(\frac{z}{C}\right), \qquad
u_p(z,t) = \frac{\dot C}{C}\,z .
$$

There is **no trajectory crossing at any Stokes number**, and therefore no $\delta$-shock: MK
is used inside its domain of validity throughout, which makes the case a clean measurement
rather than a test of how the closure fails.

Two numbers carry the physics: $\arg C$ lags $\Omega t$, and $|C|$ exceeds 1 because inertia
throws the particles outwards.

<figure>
  --8<-- "docs/vv/images/verification-vortex.svg"
</figure>

The cloud is released on the dotted gas orbit, and the thick outline is where the closed form
puts it a quarter turn later. At $St = 0.01$ it has ridden round with the gas and is still the
same size. By $St = 10$ it has covered only two thirds of the turn, and inertia has thrown it
1.83 times further from the axis, spreading it by the same factor. ICE sits inside the outline
in every panel.

<figure>
  --8<-- "docs/vv/images/verification-vortex-transition.svg"
</figure>

The same results as a path and as a curve; lines are the closed form and markers are ICE. On
the left every cloud leaves the same point and spirals outwards, the more so the heavier it
is. On the right both quantities are drawn against a continuum of Stokes numbers, not only the
four that were run, which makes those four measured points a test rather than a
demonstration.

**Results.** Four decades of Stokes number, on $96^2$ cells over a quarter turn (the gas turns
1.5708 rad):

| $St = \Omega\tau$ | 0.01 | 0.1 | 1 | 10 |
|---|---|---|---|---|
| dilation $\lvert C\rvert$, ICE | 1.015428 | 1.151011 | 1.616801 | 1.828395 |
| dilation $\lvert C\rvert$, exact | 1.015725 | 1.151573 | 1.617528 | 1.828992 |
| phase lag [rad], ICE | 0.002476 | 0.027484 | 0.345841 | 0.538749 |
| phase lag [rad], exact | 0.000310 | 0.025868 | 0.345502 | 0.538606 |

This is the transition the case exists to show: at $St \ll 1$ the cloud is a tracer, at
$St \sim 1$ it lags appreciably, and at $St \gg 1$ it barely responds to the vortex while being
flung outwards. Two more checks:

- the velocity field is checked directly, not only through where the cloud ended up: the
  mass-weighted complex slope of $u_p + iv_p$ against $z$ matches $\dot C/C$ to $10^{-4}$
  relative at every Stokes number;
- the total mass stays within $4\times10^{-7}$ of the exact solution, which conserves it
  identically because the $1/|C|^2$ factor cancels the area growth.

hydra's `test/vortex-cloud` solves this case with ICE and with IGLOO's Lagrangian parcels, and
sets both beside the closed form at the quarter turn.

!!! note "Why the low-$St$ lag is the loosest number in the table"
    At $St = 0.01$ the cloud barely dilates, so it stays as compact as it started and is the
    least resolved case on a fixed mesh. The measured lag of $2.5\times10^{-3}$ rad against an
    exact $3.1\times10^{-4}$ is numerical diffusion, not physics. It is still two orders of
    magnitude below the $St = 1$ lag, which is the statement being made.

## Convergence

Four of the cases refine the mesh against an exact solution. Together they show what the
space discretisation delivers.

<figure>
  --8<-- "docs/vv/images/verification-orders.svg"
</figure>

| Study | Meshes | $L_1$ density error | Observed order |
|---|---|---|---|
| C. sine, MUSCL / `vanleer` | 25 – 200 | $1.63\times10^{-2}$ → $2.61\times10^{-4}$ | 1.94 – 2.04 |
| C. sine, `first-order` | 25 – 200 | $1.12\times10^{-1}$ → $1.55\times10^{-2}$ | 0.92 – 0.98 |
| G. translation | 50 – 200 | $1.46\times10^{-1}$ → $1.43\times10^{-2}$ | 1.64 – 1.71 |
| H. strain | 100 – 400 | $6.08\times10^{-2}$ → $5.31\times10^{-3}$ | 1.69 – 1.83 |
| I. vortex | $48^2$ – $192^2$ | $3.55\times10^{-1}$ → $2.72\times10^{-2}$ | 1.77 – 1.94 |

Only C reaches the design order cleanly, because only C is smooth and well resolved
everywhere. The cloud cases carry a compact packet a few cells across at the coarse end, so
the limiter is active over much of the profile, and the measured order sits between 1.6 and
1.9, rising towards 2 as the packet becomes better resolved. The vortex *centroid*, which does
not depend on how well the peak is resolved, converges faster, at 1.84 – 2.28.
