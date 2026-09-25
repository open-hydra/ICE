# Code Verification

The cases on the other V&V pages compare ICE against a stored solution of its own, so
they detect change. The cases here compare it against solutions that do not come from
ICE at all — a closed form, or an accurate integration of the equation the solver is
supposed to be integrating — so they can be wrong on the first run.

They live in `test/verification/`. Each one writes its own mesh, initial condition,
boundary conditions and `input.ini` into a scratch directory, runs the solver and
compares: there is no case data in the repository and no reference to regenerate.

All of them use the MK closure. A to F and K hold the cloud uniform in space, so the exact
answer is the same in every cell and the mesh only has to be large enough to exercise
the flux loops. G to I transport a localized cloud through a prescribed frozen carrier
field, so the answer varies from cell to cell and the mesh resolution is part of what
is measured.

| | Case | What it pins down | Reference |
|---|---|---|---|
| **A** | [Stokes drag relaxation](#a-relaxation-under-stokes-drag) | The momentum source, and the order of the time scheme | Closed form |
| **B** | [Thermal relaxation](#b-thermal-relaxation) | The energy source at $Nu = 2$ | Closed form |
| **D** | [Every drag correlation](#d-every-drag-correlation) | All twelve drag laws, to $Ma = 2.6$ | Independent RK4 |
| **E** | [Every Nusselt correlation](#e-every-nusselt-correlation) | All seven heat laws, velocity and temperature relaxing together | Independent RK4 and closed form |
| **F** | [Evaporation](#f-evaporation) | All five evaporation models, the latent sink, two exact identities, the boiling clamp, the table's `Psat` | Closed form and RK4 |
| **K** | [The property table](#k-the-property-table-at-a-fixed-temperature) | Linear interpolation between the table's rows, a table that starts above 1 K, saturation past its ends | Closed form |
| **C** | [Sinusoidal advection](#c-sinusoidal-advection-on-a-periodic-mesh) | Transport alone, and the order of the space scheme | Closed form |
| **G** | [Cloud in a uniform gas](#g-a-cloud-released-into-a-uniform-gas) | Transport and drag together | Closed form |
| **H** | [Cloud in a straining gas](#h-a-cloud-in-a-straining-gas) | A non-trivial particle velocity field, and the small-Stokes limit | Closed form |
| **I** | [Cloud in a vortex](#i-a-cloud-in-a-prescribed-vortex) | Two dimensions, over four decades of Stokes number | Closed form |

The letters are the directory names under `test/verification/` and were assigned as the
cases were written; the page is grouped by what each one isolates instead, which is why
C appears after F.

## The source terms on their own

A and B strip the problem down as far as it goes — a cloud with no spatial gradients, so
the transport operator contributes nothing and each source term reduces to a single
ODE with a known solution. D and E then take the same box and run every correlation
through it.

<figure>
  {% include "vv/images/verification-relaxation.svg" %}
</figure>

### A. Relaxation under Stokes drag

A cloud at rest in a uniform gas stream, one-way coupled. Nothing varies in space, so
the momentum equation reduces to an ODE whose solution is exact:

$$
u_p(t) = u_g + (u_{p,0} - u_g)\,e^{-t/\tau_p},
\qquad
\tau_p = \frac{\rho_{al} d_p^2}{18 \mu_g}.
$$

That $\tau_p$ is not assumed: with $C_d = 24/Re$ ICE's relaxation time
$8\rho_{al}R_p/(3\rho_g C_d |\Delta u|)$ reduces to it exactly, and the test confirms
the resulting curve. With $\rho_{al} = 1000$ kg/m³, $d_p = 10^{-4}$ m and
$\mu_g = 1.8\times10^{-5}$ Pa s, $\tau_p = 0.030864$ s.

| $t/\tau_p$ | $u_p/u_g$ (ICE) | $u_p/u_g$ (exact) | Relative error |
|---|---|---|---|
| 0.1 | 0.096146407 | 0.096146510 | $1.0\times10^{-7}$ |
| 0.5 | 0.393624252 | 0.393624592 | $3.4\times10^{-7}$ |
| 1.0 | 0.632308452 | 0.632308865 | $4.1\times10^{-7}$ |
| 2.0 | 0.864802926 | 0.864803229 | $3.0\times10^{-7}$ |
| 3.0 | 0.950289178 | 0.950289346 | $1.7\times10^{-7}$ |
| 5.0 | 0.993279242 | 0.993279280 | $3.8\times10^{-8}$ |

The case also checks what must *not* move: the density, the number density and the
particle temperature stay at their initial values, and the field stays uniform to
round-off. The temperature is a real check rather than a formality — the drag work in
the energy equation exactly cancels the kinetic energy it produces, so any imbalance
there would show up as spurious heating.

**Time-step refinement.** Since the field stays uniform the space discretisation
contributes nothing, and the error is purely the time integration of the source term.
Halving the CFL number halves the step, and the observed order is the scheme's. The case
also checks the `dt-max` ceiling, by giving a run a ceiling well below its CFL limit and
confirming that its end time is an exact multiple of it — every step is `dt-max` itself,
not `cfl * dt-max`:

| CFL | 0.8 | 0.4 | 0.2 | Observed order |
|---|---|---|---|---|
| RK2 | $3.8\times10^{-7}$ | $9.5\times10^{-8}$ | $2.4\times10^{-8}$ | 2.00 |
| RK3 | $2.5\times10^{-10}$ | $3.1\times10^{-11}$ | $3.9\times10^{-12}$ | 3.00 |

### B. Thermal relaxation

The same idea for the energy equation. The particles travel with the gas, so there is
no slip, no drag force and no drag work; with `heat-transfer = Stokes` the Nusselt number is
exactly 2 and the convective exchange is again a linear relaxation:

$$
T_p(t) = T_g + (T_{p,0}-T_g)\,e^{-t/\tau_T},
\qquad
\tau_T = \frac{\rho_{al} c_s d_p^2}{12 k_g} = 0.028846\ \text{s}.
$$

Radiation is switched off (`emiss = 0`); it is the only other term in the energy
source.

| $t/\tau_T$ | $T_p$ [K] (ICE) | $T_p$ [K] (exact) | Relative error |
|---|---|---|---|
| 0.1 | 390.247598 | 390.247586 | $1.2\times10^{-7}$ |
| 1.0 | 336.744852 | 336.744805 | $4.7\times10^{-7}$ |
| 3.0 | 304.975010 | 304.974991 | $1.9\times10^{-7}$ |
| 5.0 | 300.673583 | 300.673579 | $4.3\times10^{-8}$ |

Errors are relative to the initial 100 K gap. The velocity, density and number density
must not move, and do not.

### D. Every drag correlation

ICE offers twelve drag laws. Case A verifies one of them against a closed form; this
case covers the rest, in three ways.

<figure>
  {% include "vv/images/verification-correlations.svg" %}
</figure>

Every one of the nineteen correlations in D and E sits about three decades inside
the tolerance the cases assert. The detail is what each sweep is for.

**At vanishing slip** ($Re = 6.7\times10^{-3}$) the laws that are a Stokes law plus a
correction must reproduce the exact Stokes exponential. This check is independent of
how the correlation is written in ICE:

| Law | Deviation from Stokes |
|---|---|
| Stokes | $4.1\times10^{-7}$ |
| Morsi-Alexander | $4.1\times10^{-7}$ |
| Schlichting | $2.9\times10^{-4}$ |
| Schiller-Naumann, Wen-Yu, Clift-Gauvin | $1.3\times10^{-3}$ |
| Putnam | $1.6\times10^{-3}$ |

Morsi-Alexander is piecewise and its $Re \le 0.1$ branch *is* the Stokes law, which is
why it matches to round-off; the others differ by their own correction terms, which is
the expected behaviour and not an error.

**At finite slip** every law is compared against an RK4 integration of
$du/dt = (u_g-u)/\tau(u)$ using the same correlation, evaluated independently in
Python. Two sweeps are run, at $Re = 67$ and at $Re = 1333$, the second one to reach
the branches that only switch formula above $Re = 1000$. All twelve agree with the
reference to between $5\times10^{-10}$ and $2\times10^{-6}$ of $u_g$.

**At transonic and supersonic slip** ($Ma = 1.30$ and $Ma = 2.59$) the four laws with a
compressible branch — Carlson-Hoglund, Henderson, Crowe and Hermsen — are compared the
same way, and agree to between $9\times10^{-8}$ and $2\times10^{-6}$ of $u_g$. This is
the range in which Henderson leaves its subsonic fit, crosses the linear bridge and
reaches its supersonic fit, so the case also checks that the two branches meet: the
coefficient of the bridge is fixed by requiring continuity at $Ma = 1.75$, and both
joins are continuous to $2\times10^{-9}$.

!!! note "What this check can and cannot catch"
    It catches a law that is mis-wired, mis-selected, divergent, or integrated
    incorrectly by the solver. It cannot catch a correlation that is written wrongly in
    the same way in both ICE and the reference — only the low-slip limits above, and a
    reading against the original publications, can do that.

### E. Every Nusselt correlation

The counterpart of D for convective heat exchange, over ICE's seven laws. The cloud
starts hot *and* slipping, so velocity and temperature relax together: $Re$ and $Ma$
change while the particles cool and the Nusselt number follows them. The reference
integrates the coupled pair with RK4.

At vanishing slip the four laws that tend to $Nu = 2$ must reproduce case B's exact
relaxation, and do, to between $7\times10^{-7}$ (Stokes) and $6\times10^{-3}$
(Ranz-Marshall, whose $Re^{1/2}$ correction is the largest of the four at this slip).
JAXA1 tends to zero rather than 2, and the two laws with a Mach correction keep a
finite $Ma/Re$ ratio there, so they are excluded from that check and covered by the
finite-slip comparison, `JAXA4` also by the constant-slip leg below. At $Re = 67$ all seven agree with the reference to
$2.6\times10^{-6}$ of the initial temperature gap.

A last leg holds the slip constant: under `NoDrag` the particles stay at rest in the
moving gas, so $Re = 67$, $Ma = 0.029$ and $Nu$ are frozen and the cooling is an exact
exponential. `JAXA4` meets it to $8\times10^{-6}$ of the gap at $t = 0.01$ s, which is the
RK2 error of the default step; the tolerance, $3\times10^{-4}$, is about 40 times that.
This is the leg with margin on the law's constant: a 1.4 % change of
it moves the result by $3.5\times10^{-3}$ of the gap, against $1.4\times10^{-3}$ in the
finite-slip leg.

## Phase change

### F. Evaporation

A uniform cloud of droplets at rest in a uniform, hotter gas at rest. There is no
slip, so $Re = 0$ and every correlation collapses onto its stagnant-film value.
The droplets are water, 100 µm across, at 300 K in air at 800 K and one atmosphere.

<figure>
  {% include "vv/images/verification-evaporation.svg" %}
</figure>

Lines are the reference and markers are ICE — the closed form on the left, the RK4
integration in the centre. Left: every model shrinks the droplet along a straight line
in $d^2$, at its own rate, and `ASM` lies exactly on `CEM`, dashed so that both stay
visible. Centre: released, the droplet climbs towards a wet-bulb temperature that
differs by model, because each resolves the gas-side heat differently. Right: one
combination is
rigorously isothermal, and the case uses that as its sharpest check.

**The $d^2$ law.** ICE transports $\rho_p$ and $n$; $n$ has no source term, so the
diameter comes back out of the state as $d = 2\,[0.75\rho_p/(n\pi\rho_l)]^{1/3}$.
Every model here gives $\dot m \propto d$, so

$$
\frac{\mathrm d (d^2)}{\mathrm d t} = \frac{4\dot m}{\pi d \rho_l} = -K,
\qquad K \ \text{constant},
$$

and $d^2(t) = d_0^2 - Kt$ exactly — provided the droplet temperature holds. The first
part of the case forces that with an enormous specific heat, which is the idealisation
the law is derived under, and compares $K$ against the published correlations evaluated
independently in Python:

| `evaporation` | $K$ [m²/s] | Relative error | $d^2/d_0^2$ left after 0.04 s |
|---|---|---|---|
| `d2-law` | $4.15539205\times10^{-8}$ | $7\times10^{-10}$ | 0.834 |
| `CEM`    | $5.40199730\times10^{-9}$ | $1.6\times10^{-8}$ | 0.978 |
| `CEM-B`  | $3.70421443\times10^{-9}$ | $1.6\times10^{-8}$ | 0.985 |
| `ASM`    | $5.40199730\times10^{-9}$ | $1.5\times10^{-8}$ | 0.978 |
| `TC`     | $7.81868682\times10^{-9}$ | $1.4\times10^{-8}$ | 0.969 |

`ASM` and `CEM` agree to round-off here, and must: the Frössling $Sh_0$ is 2 at
$Re = 0$, so the Abramzon-Sirignano film correction $Sh^\star = 2 + (Sh_0-2)/F(B_M)$
has nothing to correct. The case asserts that identity rather than treating the
coincidence as luck. It also checks that `evaporation = none` leaves the bulk density
untouched to round-off, that the field stays uniform, and that the number density does
not move — droplets shrink, they are never destroyed.

**The non-equilibrium interface.** `evaporation-interface = LK` depresses the surface
mole fraction by a Knudsen-layer term, and the initial slope drops 0.313 % below the
equilibrium one. It is *not* a $d^2$ law: the layer thickness goes as $1/d$ while the
rate driving it does not, so the depression deepens as the droplet shrinks. The
reference is the integrated ODE instead, and ICE matches it to $4\times10^{-13}$ of the
initial bulk density.

**The energy equation.** Freezing the temperature says nothing about the latent sink,
so the second part releases it — real specific heat, $\tau_T = 0.134$ s — and compares
against an RK4 integration of the coupled $(\rho_p, T_p)$ system. This is the part that
reaches the gas-side heat `ASM` and `TC` compute for themselves in place of the Nusselt
one. After 0.05 s:

| `evaporation` | $T_p$ [K] (ICE) | $T_p$ [K] (RK4) | $\rho_p$ error |
|---|---|---|---|
| `d2-law` | 319.519825 | 319.519825 | $2.4\times10^{-11}$ |
| `ASM`    | 337.352129 | 337.352129 | $6.5\times10^{-12}$ |
| `TC`     | 331.380563 | 331.380563 | $3.1\times10^{-12}$ |

**An exact check with nothing frozen.** Putting the $d^2$ law together with $Nu = 2$
and the Miller-Harstad-Bellan blowing factor makes the Stefan number collapse:
$b = -\tfrac32 Pr\,\tau_p\dot m/m_p$ reduces to $\ln(1+B_T)$ identically, every
property cancelling, so $f_2 = \ln(1+B_T)/B_T$ and the blown convective heat equals the
latent sink term for term. That combination is rigorously isothermal, for any gas, any
material and any droplet temperature. ICE holds $T_p$ to $8\times10^{-12}$ K over ten
thousand steps with a physical specific heat, and reproduces the $d^2$ slope to
$5\times10^{-11}$ — where the same run without the blowing factor heats by 19.5 K.

**Above the boiling point.** At 380 K the saturation pressure is 1.27 atm, above the gas
pressure, so the surface sits on the boiling clamp $X_s = 1 - 10^{-12}$ and
$B_M = 6.2\times10^{11}$. The temperature is frozen again and the run lasts until `CEM`
has lost 10 % of its mass ($1.2\times10^{-4}$ s). Every model reproduces its $d^2$ slope
to $4\times10^{-10}$, and `CEM` with the `LK` interface, which takes $X_s$ off the clamp
to 0.98, matches the integrated ODE to $10^{-13}$. The tolerance on the slope is
$10^{-5}$: one ulp of $Y_s$ at the clamp moves the rate by $2.5\times10^{-6}$.

**The saturation pressure from the property table.** The last part writes the droplet
material as a [property table](../user/initial-conditions.md#property-table) on 250 to
450 K with a `Psat` column, which replaces Clausius-Clapeyron. A specific heat of
$10^{16}$ holds $T_p$ on the 300 K row to $6\times10^{-11}$ K, so the interpolation
between rows never enters, and the reference is the integrated `CEM` ODE with $p_{sat}$
read from the table as written. A column 1.2 times the curve raises the rate by a factor
1.207, and ICE matches its reference to $5\times10^{-13}$ of the initial bulk density;
a column equal to the curve is the control, and gives what the curve gives, to
$3\times10^{-13}$. The tolerance is $10^{-8}$, that of the other integrated references;
reading the curve instead of the 1.2 column misses it by $6.6\times10^{-3}$.

!!! note "What this check can and cannot catch"
    As with D and E, a correlation written wrongly in the same way in both ICE and the
    reference would pass. What is not vulnerable to that is the $\dot m \propto d$
    scaling, the `ASM`–`CEM` identity at $Re = 0$ and the isothermal identity above:
    those follow from the published forms and would break under any transcription
    error.

## Material properties

### K. The property table at a fixed temperature

Case A's cloud with its density taken from `INPUT/part-properties.dat` instead of the
INI. With `heat-transfer = NoHeat` and no radiation nothing heats the particles, and the
drag work cancels the kinetic energy it produces, so $T_p$ keeps its initial value and
the table is read at one temperature for the whole run. The velocity is case A's
exponential with $\tau_p = \rho_{al}(T_p)\,d_p^2/(18\mu_g)$, and $\rho_{al}(T_p)$ is what
the case measures: ICE recovers the radius from $\rho_p$ and $n$ through the table's
density, so at fixed $(\rho_p, n)$ the relaxation time goes as $\rho_{al}^{1/3}$.

| Leg | Table | $T_p$ | $\rho_{al}(T_p)$ | Relative error on $u_p(\tau_p)$ |
|---|---|---|---|---|
| K1 | 1 to 5000 K; $\rho$ = 2000 up to 300 K, 1000 from 301 K | 300.4 K | 1600 (linear) | $9.8\times10^{-8}$ |
| K1c | the same rows, $\rho$ = 1600 everywhere | 300.4 K | 1600 | $9.7\times10^{-8}$ |
| K6 | 280 to 400 K, constant | 300.4 K | 1500 | $9.7\times10^{-8}$ |
| K7 | 280 to 400 K, $\rho = 1500 + 5\,(T - 280)$ | 250 K | 1500 (the end value) | $9.7\times10^{-8}$ |

The tolerance is $10^{-4}$. RK2 at $\Delta t = \tau_p/1000$ leaves
$e^{-1}(\Delta t/\tau_p)^2/6 = 6.1\times10^{-8}\,u_g$ on $u_p(\tau_p)$, that is
$9.7\times10^{-8}$ of $u_p(\tau_p)$ itself, which is what the four legs measure. The
nearest row would give K1 $\rho_{al} = 2000$ and put $u_p(\tau_p)$ 4.3 % low; a linear
extrapolation below $T_{min}$ would put K7 2.0 % high, and the row at $T_{max}$ 6.5 % low.
K1c is K1's control, the same run with nothing to interpolate. $T_p$ must stay within
$10^{-4}$ K of its initial value, well under the 0.1 K that would change the nearest row;
it moves by round-off.

## Clouds carried by the gas

The remaining four give the transport operator something to do. C is transport alone;
G adds a source term that is uniform in space; H one that is not; I does it in two
dimensions. All four compare against a closed form, and none of them can be passed by a
solver that transports correctly but couples its source terms wrongly.

### C. Sinusoidal advection on a periodic mesh

Cases A and B leave the transport operator idle. Here there is no gas at all (0-way
coupling), the cloud moves at a uniform $U = 1$ m/s and the density carries a sine on a
periodic unit domain, so free transport gives

$$
\rho_p(x,t) = \rho_0 + A \sin\!\big[k(x - Ut)\big], \qquad k = 2\pi/L.
$$

The comparison is against the exact *cell average*, so the measured order is not
polluted by the difference between a cell average and a point value. The `dt-max`
ceiling holds $\Delta t$ at the same value on all four meshes, so the time error is
common to them and what the refinement measures is the space discretisation.

Two reconstructions are refined side by side, and both reach their design order — the
numbers are in [Convergence](#convergence-in-one-picture) with the other studies. The
case also checks that mass is conserved to
round-off across the periodic interface and that the velocity stays uniform — a
pressureless cloud moving at one speed has nothing that could change it.

### G. A cloud released into a uniform gas

Cases A and B leave the transport operator idle; case C leaves the source terms idle.
This one is the product of the two, and it is the smallest problem in which getting
either one wrong — or coupling them wrongly — changes the answer.

A frozen gas moves at $U_g = 10$ m/s and a localized cloud starts at rest. The Stokes
relaxation time $\tau = \rho_{al} d_p^2 / 18\mu$ depends only on the diameter, which
stays uniform because $n$ is seeded proportional to $\rho_p$ and the two obey the same
transport equation. Every particle therefore feels the same drag, the particle velocity
stays uniform in space for all time, and the cloud translates rigidly:

$$
u_p(t) = U_g\big(1 - e^{-t/\tau}\big), \qquad
X(t) = U_g t - U_g\tau\big(1 - e^{-t/\tau}\big), \qquad
\rho_p(x,t) = \rho_p(x - X(t), 0).
$$

At $t = 3\tau$ the cloud has moved 0.633 m where the gas would have carried it 0.926 m.

<figure>
  {% include "vv/images/verification-translation.svg" %}
</figure>

Grey is the exact profile, red is ICE, and the labels mark the elapsed time in units
of $\tau$. The cloud keeps its shape exactly — only the top hat smears, and only
because a discontinuity must. What sets its position is the velocity curve on the
right, which is case A's exponential recovered from a cloud that is now moving.

**The uniformity of $u_p$ is the sharper statement.** It is not uniform at finite
$\Delta t$: the drag source is evaluated on the density that the flux divergence is
about to change, which leaves a spread wherever the density gradient is steep. That
spread is mesh-independent, survives a first-order reconstruction, and falls at
second order in $\Delta t$ —

| `dt-max` | $2\times10^{-5}$ | $1\times10^{-5}$ | $5\times10^{-6}$ | Observed order |
|---|---|---|---|---|
| $(\max u_p - \min u_p)/U_g$ | $9.22\times10^{-6}$ | $2.32\times10^{-6}$ | $5.81\times10^{-7}$ | 1.99 – 2.00 |

— which is the statement that MK's momentum flux and its drag source are consistent.
The mass-weighted $u_p$ matches the exponential to $3\times10^{-8}$ m/s throughout.

A top hat is run alongside the Gaussian for the invariants that survive any amount of
smearing: mass is conserved to $2\times10^{-16}$ across the periodic interface, and the
centroid sits within 0.007 cells of $X(t)$.

### H. A cloud in a straining gas

Case G moves a cloud through a gas of one velocity, so the particle velocity stays
uniform and the cloud only translates. Here the frozen gas carries a uniform strain,
$u_g(x) = U_0 + ax$ with $U_0 = 2$ m/s and $a = 6\ \mathrm{s^{-1}}$, and the cloud can
no longer move rigidly: the particles lag the gas acceleration, which is the defining
behaviour of a dispersed phase. It is the first case with a non-trivial particle
velocity *field*, so the first that can detect a momentum flux which is wrong in a way
a uniform velocity hides.

With $\xi = x + U_0/a$ the gas velocity is $a\xi$, and Stokes drag with a uniform
diameter gives one linear equation along a trajectory:

<figure>
  {% include "vv/images/verification-strain.svg" %}
</figure>

Grey is exact, red is ICE. The cloud is stretched and diluted by the same factor, the
particles run slower than the gas everywhere — the shaded band is the slip — and that
slip approaches the textbook asymptote from below as the Stokes number falls.

$$
\tau\ddot\xi + \dot\xi - a\xi = 0, \qquad
\lambda_\pm = \frac{-1 \pm \sqrt{1 + 4a\tau}}{2\tau}.
$$

Seeding the cloud at the local gas velocity, $\dot\xi(0) = a\xi(0)$, makes the solution
separable — every particle is displaced by the same factor $\xi(t) = C(t)\,\xi_0$, with
$C = Ae^{\lambda_+t} + Be^{\lambda_-t}$, $A + B = 1$ and $A\lambda_+ + B\lambda_- = a$.
The map $x_0 \mapsto x(t)$ is affine, which gives the whole Eulerian solution in closed
form:

$$
\rho_p(x,t) = \frac{1}{C}\,\rho_{p0}\!\left(\frac{\xi}{C}\right), \qquad
u_p(x,t) = \frac{\dot C}{C}\,\xi ,
$$

exact at every Stokes number, with no crossing at any time since $C > 0$. At $St = a\tau = 0.185$ and $t = 0.1$ s the cloud is stretched by $C = 1.7096$.

The centroid lands within 0.004 cells of $C\xi_0$ on the finest mesh, and the
mass-weighted slope of the velocity field matches $\dot C/C$ to $2\times10^{-6}$
relative. The refinement is in [Convergence](#convergence-in-one-picture).

**The asymptotic check.** For small Stokes number the standard result is
$u_p \approx u_g - \tau\,Du_g/Dt$, which here is a slip of $a^2\tau\xi$. Expanding the
exact slip coefficient $a - \dot C/C$ in $St$ gives $a^2\tau(1 - 2St + O(St^2))$, so the
departure from the asymptote must be proportional to $St$ with coefficient 2:

| $St = a\tau$ | 0.200 | 0.100 | 0.050 | 0.025 |
|---|---|---|---|---|
| measured slip / $a^2\tau$ | 0.7167 | 0.8385 | 0.9109 | 0.9529 |
| $(1 - \text{ratio})/St$ | 1.42 | 1.62 | 1.78 | 1.89 |

The ratio approaches 1 monotonically and the last column climbs towards 2, which is
what pins the solver to the right asymptotic behaviour rather than merely to a
plausible one. Against the exact $\dot C/C$ — which holds asymptotically or not — the
velocity field agrees to $10^{-5}$ relative at every one of those Stokes numbers.

### I. A cloud in a prescribed vortex

The frozen gas is in solid-body rotation, $u_g = -\Omega y$, $v_g = \Omega x$, and a
compact cloud is released into it seeded at the local gas velocity. This is the first
two-dimensional verification case: G and H leave the $j$-direction fluxes idle, and a
momentum flux wrong across directions, or a metric wrong on the $j$ faces, cannot show
up in either.

In the complex plane $z = x + iy$ the carrier field is $u_g = i\Omega z$, so the
trajectory equation is case H over the complex field:

$$
\tau\ddot z + \dot z = i\Omega z, \qquad
\lambda_\pm = \frac{-1 \pm \sqrt{1 + 4i\Omega\tau}}{2\tau}.
$$

Seeding $\dot z(0) = i\Omega z(0)$ makes it separable again — $z(t) = C(t)\,z_0$ with
$C = Ae^{\lambda_+t} + Be^{\lambda_-t}$, $A + B = 1$, $A\lambda_+ + B\lambda_- = i\Omega$.
The map is a rotation composed with a dilation: conformal, and non-singular for all
time. So the Eulerian solution is closed form,

$$
\rho_p(z,t) = \frac{1}{|C|^2}\,\rho_{p0}\!\left(\frac{z}{C}\right), \qquad
u_p(z,t) = \frac{\dot C}{C}\,z ,
$$

with **no trajectory crossing at any Stokes number** and therefore no $\delta$-shock:
MK is used inside its domain of validity throughout, which is what makes the case a
clean measurement rather than a test of how the closure fails.

Two numbers carry the physics. $\arg C$ lags $\Omega t$, and $|C|$ exceeds 1 because
inertia throws the particles outwards.

<figure>
  {% include "vv/images/verification-vortex.svg" %}
</figure>

The cloud is released on the dotted gas orbit, and the thick outline is where the
closed form puts it a quarter turn later. At $St = 0.01$ it has ridden round with the
gas and is still the same size; by $St = 10$ it has covered only two thirds of the turn
and inertia has thrown it 1.83 times further from the axis, spreading it by the same
factor. ICE sits inside the outline in every panel.

<figure>
  {% include "vv/images/verification-vortex-transition.svg" %}
</figure>

The same thing as a path and as a curve; lines are the closed form and markers are ICE
throughout. On the left every cloud leaves the same point and spirals outwards, the
more so the heavier it is. On the right both quantities are drawn against a continuum
of Stokes numbers rather than only the four that were run, which is what makes those
four measured points a test rather than a demonstration.

The sweep over four decades of Stokes number, on $96^2$ over a quarter turn (the gas
turns 1.5708 rad):

| $St = \Omega\tau$ | 0.01 | 0.1 | 1 | 10 |
|---|---|---|---|---|
| dilation $\lvert C\rvert$, ICE | 1.015428 | 1.151011 | 1.616801 | 1.828395 |
| dilation $\lvert C\rvert$, exact | 1.015725 | 1.151573 | 1.617528 | 1.828992 |
| phase lag [rad], ICE | 0.002476 | 0.027484 | 0.345841 | 0.538749 |
| phase lag [rad], exact | 0.000310 | 0.025868 | 0.345502 | 0.538606 |

This is the transition the case exists to show: at $St \ll 1$ the cloud is a tracer, at
$St \sim 1$ it lags appreciably, and at $St \gg 1$ it barely responds to the vortex at
all while being flung outwards. The velocity field is checked directly rather than only
through where the cloud ended up — the mass-weighted complex slope of $u_p + iv_p$
against $z$ matches $\dot C/C$ to $10^{-4}$ relative at every Stokes number — and the
total mass stays within $4\times10^{-7}$ of the exact solution, which conserves it
identically because the $1/|C|^2$ factor cancels the area growth.

!!! note "Why the low-$St$ lag is the loosest number in the table"
    At $St = 0.01$ the cloud barely dilates, so it stays as compact as it started and is
    the least resolved case on a fixed mesh; the measured lag of $2.5\times10^{-3}$ rad
    against an exact $3.1\times10^{-4}$ is numerical diffusion, not physics. It is still
    two orders of magnitude below the $St = 1$ lag, which is the statement being made.

## Convergence, in one picture

Four of the cases refine the mesh against an exact solution. Collecting them shows what
the space discretisation actually delivers.

<figure>
  {% include "vv/images/verification-orders.svg" %}
</figure>

| Study | Meshes | $L_1$ density error | Observed order |
|---|---|---|---|
| C. sine, MUSCL / `vanleer` | 25 – 200 | $1.63\times10^{-2}$ → $2.61\times10^{-4}$ | 1.94 – 2.04 |
| C. sine, `first-order` | 25 – 200 | $1.12\times10^{-1}$ → $1.55\times10^{-2}$ | 0.92 – 0.98 |
| G. translation | 50 – 200 | $1.46\times10^{-1}$ → $1.43\times10^{-2}$ | 1.64 – 1.71 |
| H. strain | 100 – 400 | $6.08\times10^{-2}$ → $5.31\times10^{-3}$ | 1.69 – 1.83 |
| I. vortex | $48^2$ – $192^2$ | $3.55\times10^{-1}$ → $2.72\times10^{-2}$ | 1.77 – 1.94 |

Only C reaches the design order cleanly, because only C is smooth everywhere and
everywhere well resolved. The cloud cases carry a compact packet a few cells across at
the coarse end, so the limiter is active over much of the profile and the measured
order sits between 1.6 and 1.9, rising towards 2 as the packet becomes better resolved.
The vortex *centroid*, which does not depend on how well the peak is resolved,
converges faster, at 1.84 – 2.28.

## Running them

```bash
ctest --test-dir build -L verification --output-on-failure
```

A, B, C and K also carry the `fast` label, so the pre-push hook runs them. The rest are in
the default tier only: D, E and F because each correlation is a separate run of the
solver, and G, H and I because each refinement level is. I is the longest, being two
dimensional. `ICE_KEEP_WORK=1` keeps the generated case directories instead of deleting
them on success, which is how to look at what a failing case actually produced.

The figures on this page come from

```bash
cd test/verification && python3 plot_vv.py [figure ...]
```

which takes a few minutes, most of it in the vortex. It does not read the output of a
ctest run — these cases build their inputs on the fly and remove them again, so there
is nothing left behind to plot — and runs what each figure needs itself. It also prints
the numbers quoted above. Named figures are `relaxation`, `correlations`,
`evaporation`, `translation`, `strain`, `vortex` and `orders`; with no argument it
draws all of them.

