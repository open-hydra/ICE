# Code Verification

The cases on the other V&V pages compare ICE against a stored solution of its own, so
they detect change. The cases here compare it against solutions that do not come from
ICE at all — a closed form, or an accurate integration of the equation the solver is
supposed to be integrating — so they can be wrong on the first run.

They live in `test/verification/`. Each one writes its own mesh, initial condition,
boundary conditions and `input.ini` into a scratch directory, runs the solver and
compares: there is no case data in the repository and no reference to regenerate.

All of them use the MK closure on a uniform cloud, so the exact answer is the same in
every cell and the mesh only has to be large enough to exercise the flux loops.

## A. Relaxation under Stokes drag

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

## B. Thermal relaxation

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

## C. Sinusoidal advection on a periodic mesh

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

| Cells | 25 | 50 | 100 | 200 | Observed order |
|---|---|---|---|---|---|
| MUSCL / `vanleer` | $1.63\times10^{-2}$ | $4.25\times10^{-3}$ | $1.07\times10^{-3}$ | $2.61\times10^{-4}$ | 1.94 – 2.04 |
| `first-order` | $1.12\times10^{-1}$ | $5.94\times10^{-2}$ | $3.06\times10^{-2}$ | $1.55\times10^{-2}$ | 0.92 – 0.98 |

Errors are the $L_1$ density error relative to the wave amplitude, at $t = 0.25$ s.
MUSCL reaches its design order on this smooth profile, and the first-order
reconstruction reaches its own. The case also checks that mass is conserved to
round-off across the periodic interface and that the velocity stays uniform — a
pressureless cloud moving at one speed has nothing that could change it.

## D. Every drag correlation

ICE offers twelve drag laws. Case A verifies one of them against a closed form; this
case covers the rest, in three ways.

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

## E. Every Nusselt correlation

The counterpart of D for convective heat exchange, over ICE's seven laws. The cloud
starts hot *and* slipping, so velocity and temperature relax together: $Re$ and $Ma$
change while the particles cool and the Nusselt number follows them. The reference
integrates the coupled pair with RK4.

At vanishing slip the four laws that tend to $Nu = 2$ must reproduce case B's exact
relaxation, and do, to between $5\times10^{-7}$ (Stokes) and $6\times10^{-3}$
(Ranz-Marshall, whose $Re^{1/2}$ correction is the largest of the four at this slip).
JAXA1 tends to zero rather than 2, and the two laws with a Mach correction keep a
finite $Ma/Re$ ratio there, so they are excluded from that check and covered only by
the finite-slip comparison. At $Re = 67$ all seven agree with the reference to
$1.7\times10^{-6}$ of the initial temperature gap.

## F. Evaporation

A uniform cloud of droplets at rest in a uniform, hotter gas at rest. There is no
slip, so $Re = 0$ and every correlation collapses onto its stagnant-film value.
The droplets are water, 100 µm across, at 300 K in air at 800 K and one atmosphere.

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

!!! note "What this check can and cannot catch"
    As with D and E, a correlation written wrongly in the same way in both ICE and the
    reference would pass. What is not vulnerable to that is the $\dot m \propto d$
    scaling, the `ASM`–`CEM` identity at $Re = 0$ and the isothermal identity above:
    those follow from the published forms and would break under any transcription
    error.

## Running them

```bash
ctest --test-dir build -L verification --output-on-failure
```

A, B and C also carry the `fast` label, so the pre-push hook runs them. D, E and F take
longer, because each correlation is a separate run of the solver, and are in the
default tier only.
