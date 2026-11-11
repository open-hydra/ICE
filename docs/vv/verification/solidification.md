# Solidification

Two cases on drops that freeze and melt: N in a closed cell, where the heat source acts
alone, and O in a steady stream, where the freezing is carried by the flow. Most legs run on
MK; some repeat on IG and AG.

| Case | ctest | What it pins down | Reference |
|---|---|---|---|
| [N](#n-a-closed-cell) | `SolidCooling` | Supercooling, nucleation, the plateau, the solid, melting | Closed form |
| [O](#o-a-steady-stream) | `SolidAdvection`, `SolidAdvectionMPI` | The nucleation front and the plateau in a steady stream | Closed form (a recursion at first order) |

## The material

Both cases use IGLOO's solid-box material (`tests/solidification/solid-box` there):

| Property | Value |
|---|---|
| Drop diameter | 30 µm |
| Density $\rho$ | 2950 kg/m³ |
| Specific heats $c_l$, $c_s$ | 1250 and 600 J/(kg K) |
| Latent heat of fusion $h_{fus}$ | 1.07 MJ/kg |
| Melting temperature $T_m$ | 2327 K |
| Nucleation temperature $T_{nuc}$ | 0.8 $T_m$ = 1861.6 K |
| Gas | $k_g$ = 0.026 W/(m K), $Nu = 2$ (`heat-transfer = Stokes`), `drag = NoDrag`, emissivity 0 |

With $m = \rho\pi d^3/6$ and $A = \pi d k_g Nu$, every regime is a closed form:

- the liquid and the solid relax on $\tau = m c/A$ (10.637 ms and 5.106 ms);
- on the plateau, the solid fraction $f$ moves at $A(T_g - T_m)/(m h_{fus})$;
- nucleation from $T_{nuc}$ lands on $f_0 = c_l(T_m - T_{nuc})/h_{fus}$ = 0.543692 at $T_m$
  (or, for $f_0 \ge 1$, on the solid at the temperature that keeps the energy).

## N. A closed cell

**Setup.** Eight uniform cells moving with the gas, so each integrates the heat source alone;
RK2 with $\Delta t = 10^{-5}$ s. The closed form is taken at the solution's own time.

| Leg | Setup | What it checks | Measured |
|---|---|---|---|
| N1 | 2400 K in a 600 K gas | The liquid at $t_n/2$; the plateau at three times, $T_p = T_m$ exactly; the solid at 10 and 15 ms | liquid $4\times10^{-5}$ K; plateau $f$ $2.3\times10^{-4}$ (bound $5.1\times10^{-4}$, the stage-end nucleation); solid 0.20 and 0.074 K (bounds 0.54 and 0.21 K) |
| N2 | `h-fus = 4e5` ($f_0$ = 1.4544) | The jump lands on the solid at 2024.083 K, then relaxes on $\tau_s$ | 0.13 K (bound 0.35 K) |
| N3 | A solid at 1500 K in a 3000 K gas | It reaches $T_m$ at 4.092 ms and melts there: $T_p = T_m$, $f$ = 0.170 at $3\tau_s$ | $f$ $7\times10^{-8}$; MK, IG and AG |
| N4 | A mush ($T_m$, $f$ = 0.8, $\chi$ = 1) in a 3000 K gas | $f$ falls to 0 at 10.82 ms, then the liquid with $\chi = 0$ | $f$ $5\times10^{-14}$; liquid $7\times10^{-5}$ K |
| N5 | N1 and N3 with IG and AG | The same numbers | the same |
| N6 | N1 stopped on the plateau, restarted in place to 15 ms | The straight run's state | identical |
| N7 | `A 1 solidification=on …` / `B 1` ($c_s$ = 1000) | Each material its own history; B's columns those of a run with B alone | identical |
| N8 | N1 on a `Cp` column $c_l = 50 + 0.5\,T$ | The liquid on its exact $t(T)$, $t_n$ = 3.349 ms; $f_0$ = 0.477209; plateau and solid as N1 | liquid $1.2\times10^{-5}$ K; plateau $f$ $2.2\times10^{-4}$ (bound $5.1\times10^{-4}$); solid 0.18 and 0.068 K (bounds 0.53 and 0.21 K) |

N8's liquid follows
$t(T) = \frac{m}{A}\big[b(T_0 - T) + (a + bT_g)\ln\frac{T_0 - T_g}{T - T_g}\big]$, and its
nucleation lands on $f_0 = (e_\ell(T_m) - e_\ell(T_{nuc}))/h_{fus}$.

In every leg, $\chi$ is exactly 0 on the liquid and 1 on the nucleated states, the eight cells
agree to $10^{-10}$, and the velocity, density and number density do not move.

## O. A steady stream

**Setup.** The same material enters a 0.15 m box at 10 m/s through a code-402 inlet, as a
steady run.

**The reference.** At first order the steady state is a recursion from the inlet, with
$a = (A/m)\Delta x/u$:

- while the liquid stays above $T_{nuc}$, $e = (e_{up} + a T_g)/(1 + a/c_l)$;
- otherwise the cell is at the equilibrium state (solid, plateau or liquid), with $\chi$ that of
  the inflow between $T_{nuc}$ and $T_m$.

**A fixed point.** Every run is compared with a twin run at $40\,n_x + 1$ iterations: an odd
gap, and not a multiple of 3, so that a cycle of period 2 or 3 shows. The twins must agree in
$T_p$, $f$ and $\chi$ cell by cell to $10^{-12}$.

| Leg | Setup | What it checks |
|---|---|---|
| O0 | `T-melt = 200` (every cell liquid), MUSCL; MK/Saurel and IG/Rusanov | The plain run's blocks byte for byte; $f$ and $\chi$ all zero |
| O1 | First order, $n_x$ = 120, gas 600 K | Every cell on the recursion, and the energy balance |
| O2 | MUSCL, $n_x$ = 120, 240, 480 | The front in the cell of $x_n$ = 37.805 mm; $T_p = T_m$ exactly and $\chi = 1$ on the plateau; the liquid error over $x < x_n/2$ falling at second order |
| O3 | IG/Rusanov, MUSCL, $n_x$ = 240 | O2's front; $\chi \le 10^{-4}$ upstream |
| O3b | IG/Rusanov, first order, O1's stream at 0.1 m/s, where the inlet's $P$ gives $c$ = 2.5 $u$ and 4 $u$ | Every cell on an independent Rusanov model of $e$ and $\chi$; the front one and two cells above the upwind one |
| O4 | A solid at 1500 K into a 3000 K gas, `h-fus = 4e5`, from an empty domain; MK and IG, $n_x$ = 120, 240 | The melting recursion: solid to $T_m$, the plateau, the liquid; no cell with $f > 0$ above $T_m$ |
| O4u | O4's injection without solidification into a near-empty domain (below) | Both runs fill every cell and agree per cell in $\rho_p$ and $T_p$ |
| O5 | A steady run with the gas at 300 K, then restarted in place with the gas at 1000 K | The front moves back to the 1000 K run's cell: no lock |
| O6 | O2's stream at $n_x$ = 120 on two blocks joined by a connection on the plateau; first order (RK2) and MUSCL (RK3); in an MPI build also on two ranks | Fixed points; two blocks against one; two ranks against the serial two-block run |

| Leg | Measured |
|---|---|
| O0 | identical |
| O1 | front in cell 31, exit at 914.2832 K; $T$ $1.1\times10^{-14}$, $f$ $2.9\times10^{-15}$; balance $3.2\times10^{-15}$ |
| O2 | fronts 31, 61, 122; exits 906.8681, 906.8465, 907.3422 K; liquid errors 0.0290, 0.0075, 0.0019 K; every $n_x$ a fixed point |
| O3 | front 61; $\chi$ $2.5\times10^{-5}$ upstream |
| O3b | fronts 30 and 29; $T$ $8.4\times10^{-14}$, $1.2\times10^{-13}$ |
| O4 | first mush/liquid cells 34/74 and 66/147; exits 2608.6368 and 2610.1672 K; MK $4.3\times10^{-14}$, IG $5.0\times10^{-5}$ in $\chi$ (tolerances $10^{-9}$ and $10^{-4}$) |
| O4u | $2.4\times10^{-14}$, $9.5\times10^{-15}$, $8.3\times10^{-14}$ (tolerance $2\times10^{-12}$) |
| O5 | fronts 26, then 42; per cell $4.5\times10^{-14}$ from a fresh run |
| O6 | first order value for value; MUSCL $T$ $2.2\times10^{-15}$, $f$ $3.3\times10^{-16}$ (tolerances $10^{-13}$ and $10^{-14}$); two ranks value for value |

Notes on the legs:

- **O2.** The liquid error is measured over $x < x_n/2$ because the injected state sits at the
  inlet face. Up to $2\Delta x$ before $x_n$, where the limiter acts at the front's kink, the
  errors are 0.054, 0.022 and 0.0019 K.
- **O3.** The $\chi$ upstream is the Rusanov leak with van Leer ($5.0\times10^{-5}$ at first
  order).
- **O4u.** The domain starts at $10^{-20}$ kg/m³; IG runs at $n_x$ = 120 and 240, AG at 120,
  all at first order. Each runs twice: with the initial condition's pseudo-pressure scaled with
  its content, and with the inlet's $10^{-6}$ ($P/\rho_p = 10^{14}$ m²/s², a signal speed of
  $1.7\times10^{7}$ m/s). The tolerance is ten times the largest twin floor of the scaled runs
  ($5.2\times10^{-15}$, $3.9\times10^{-16}$, $1.8\times10^{-13}$). What makes the two agree is
  the cap of a near-empty cell's dispersion at 1 m²/s². Without it the unscaled IG run ends at
  $\rho_p = 6.4\times10^{247}$ ($n_x$ = 120), and the other two stop on an invalid state.
- **O6.** Under MUSCL the two-block run differs in the last bits because the connection face's
  lengths differ from the interior's.

**Why $T_\ell$ is reconstructed.** The fixed-point test is what a reconstruction of $T_p$
fails: the limiter then sees the recalescence jump, and the front cell alternates between
liquid and nucleated, or its $\chi$ and the $f$ of the cells after it cycle (O2 at
$n_x$ = 480, period 2). With $T_\ell$ reconstructed, O2's stream is a fixed point for RK2 and
RK3, CFL 0.4 and 0.8, and $n_x$ = 240 and 480. The front sits in the cell of $x_n$, and the
exit temperature is the same for every scheme and CFL number: 905.4689 K at $n_x$ = 240 and
906.1498 K at 480. A steady nucleation front does not need a first-order reconstruction.
