# Evaporation and Boiling

Two cases on the mass transfer between the drop and the gas, both on a uniform cloud at rest
in a gas at rest. F runs every evaporation model against the $d^2$ law and against independent
integrations. M checks what happens inside the drop: its temperature profile, the boiling
branch, and the surface temperature the probe writes.

| Case | ctest | What it pins down | Reference |
|---|---|---|---|
| [F](#f-evaporation) | `Evaporation` | All five evaporation models, the latent sink, two exact identities, the boiling clamp, the table's `Psat`, `Lv`, `Mv` and `Cpv` | Closed form and RK4 |
| [M](#m-drop-heating-and-boiling) | `DropHeating` | The parabolic temperature profile inside the drop (`P2T`), the boiling branch (`ZGR`), the surface temperature written by the probe | ICE's time scheme on the same ODEs, and IGLOO |

## F. Evaporation

**Setup.** A uniform cloud of droplets at rest in a uniform, hotter gas at rest. There is no
slip, so $Re = 0$ and every correlation collapses onto its stagnant-film value. The droplets
are water, 100 µm across, at 300 K, in air at 800 K and one atmosphere.

<figure>
  --8<-- "docs/vv/images/verification-evaporation.svg"
</figure>

Lines are the reference and markers are ICE: the closed form on the left, the RK4 integration
in the centre.

- **Left:** every model shrinks the droplet along a straight line in $d^2$, at its own rate.
  `ASM` lies exactly on `CEM`, and is dashed so that both stay visible.
- **Centre:** released, the droplet climbs towards a wet-bulb temperature that differs by
  model, because each resolves the gas-side heat differently.
- **Right:** one combination is rigorously isothermal, and the case uses it as its sharpest
  check.

The case has seven parts:

| Part | What it checks | Reference |
|---|---|---|
| [The $d^2$ law](#the-d2-law) | Each model's rate, the temperature frozen | Published correlations, in Python |
| [The non-equilibrium interface](#the-non-equilibrium-interface) | `evaporation-interface = LK` | Integrated ODE |
| [The energy equation](#the-energy-equation) | The latent sink and the gas-side heat | RK4 of the coupled $(\rho_p, T_p)$ system |
| [An exact check with nothing frozen](#an-exact-check-with-nothing-frozen) | $d^2$ law, $Nu = 2$ and the blowing factor together | An exact identity |
| [Above the boiling point](#above-the-boiling-point) | The boiling clamp | $d^2$ slope; integrated ODE |
| [`Psat` from the property table](#psat-from-the-property-table) | The table's saturation pressure | Integrated `CEM` ODE |
| [`Lv`, `Mv` and `Cpv` from the property table](#lv-mv-and-cpv-from-the-property-table) | The table's latent heat and vapour properties | RK4 |

### The $d^2$ law { #the-d2-law }

ICE transports $\rho_p$ and $n$, and $n$ has no source term, so the diameter comes back out
of the state as $d = 2\,[0.75\rho_p/(n\pi\rho_l)]^{1/3}$. Every model here gives
$\dot m \propto d$, so

$$
\frac{\mathrm d (d^2)}{\mathrm d t} = \frac{4\dot m}{\pi d \rho_l} = -K,
\qquad K \ \text{constant},
$$

and $d^2(t) = d_0^2 - Kt$ exactly, provided the droplet temperature holds. This part forces
that with an enormous specific heat, the idealisation the law is derived under, and compares
$K$ against the published correlations evaluated independently in Python:

| `evaporation` | $K$ [m²/s] | Relative error | $d^2/d_0^2$ left after 0.04 s |
|---|---|---|---|
| `d2-law` | $4.15539205\times10^{-8}$ | $7\times10^{-10}$ | 0.834 |
| `CEM`    | $5.40199730\times10^{-9}$ | $1.6\times10^{-8}$ | 0.978 |
| `CEM-B`  | $3.70421443\times10^{-9}$ | $1.6\times10^{-8}$ | 0.985 |
| `ASM`    | $5.40199730\times10^{-9}$ | $1.5\times10^{-8}$ | 0.978 |
| `TC`     | $7.81868682\times10^{-9}$ | $1.4\times10^{-8}$ | 0.969 |

`ASM` and `CEM` agree to round-off, and must: the Frössling $Sh_0$ is 2 at $Re = 0$, so the
Abramzon-Sirignano film correction $Sh^\star = 2 + (Sh_0-2)/F(B_M)$ has nothing to correct.
The case asserts that identity rather than treating the coincidence as luck. It also checks
that:

- `evaporation = none` leaves the bulk density untouched to round-off;
- the field stays uniform;
- the number density does not move: droplets shrink, they are never destroyed.

### The non-equilibrium interface

`evaporation-interface = LK` depresses the surface mole fraction by a Knudsen-layer term, and
the initial slope drops 0.313 % below the equilibrium one. It is *not* a $d^2$ law: the layer
thickness goes as $1/d$ while the rate driving it does not, so the depression deepens as the
droplet shrinks. The reference is the integrated ODE instead, and ICE matches it to
$4\times10^{-13}$ of the initial bulk density.

### The energy equation

Freezing the temperature says nothing about the latent sink. This part releases it (a real
specific heat, $\tau_T = 0.134$ s) and compares against an RK4 integration of the coupled
$(\rho_p, T_p)$ system. It is the part that reaches the gas-side heat that `ASM` and `TC`
compute for themselves in place of the Nusselt one. After 0.05 s:

| `evaporation` | $T_p$ [K] (ICE) | $T_p$ [K] (RK4) | $\rho_p$ error |
|---|---|---|---|
| `d2-law` | 319.519825 | 319.519825 | $2.4\times10^{-11}$ |
| `ASM`    | 337.352129 | 337.352129 | $6.5\times10^{-12}$ |
| `TC`     | 331.380563 | 331.380563 | $3.1\times10^{-12}$ |

### An exact check with nothing frozen

The $d^2$ law together with $Nu = 2$ and the Miller-Harstad-Bellan blowing factor makes the
Stefan number collapse: $b = -\tfrac32 Pr\,\tau_p\dot m/m_p$ reduces to $\ln(1+B_T)$
identically, every property cancelling. So $f_2 = \ln(1+B_T)/B_T$, and the blown convective
heat equals the latent sink term for term. The combination is rigorously isothermal, for any
gas, any material and any droplet temperature.

| Check | Measured |
|---|---|
| $T_p$ held, over ten thousand steps with a physical specific heat | to $8\times10^{-12}$ K |
| The $d^2$ slope | to $5\times10^{-11}$ |
| The same run without the blowing factor | heats by 19.5 K |

### Above the boiling point

At 380 K the saturation pressure is 1.27 atm, above the gas pressure, so the surface sits on
the boiling clamp $X_s = 1 - 10^{-12}$ and $B_M = 6.2\times10^{11}$. The temperature is frozen
again, and the run lasts until `CEM` has lost 10 % of its mass ($1.2\times10^{-4}$ s).

- Every model reproduces its $d^2$ slope to $4\times10^{-10}$. The tolerance on the slope is
  $10^{-5}$: one ulp of $Y_s$ at the clamp moves the rate by $2.5\times10^{-6}$.
- `CEM` with the `LK` interface, which takes $X_s$ off the clamp to 0.98, matches the
  integrated ODE to $10^{-13}$.

### `Psat` from the property table { #psat-from-the-property-table }

The droplet material is written as a
[property table](../../user/initial-conditions.md#property-table) on 250 to 450 K with a
`Psat` column, which replaces Clausius-Clapeyron. A specific heat of $10^{16}$ holds $T_p$ on
the 300 K row to $6\times10^{-11}$ K, so the interpolation between rows never enters. The
reference is the integrated `CEM` ODE with $p_{sat}$ read from the table as written.

| Column | Effect | ICE against its reference |
|---|---|---|
| 1.2 times the curve | rate × 1.207 | $5\times10^{-13}$ of the initial bulk density |
| Equal to the curve (control) | what the curve gives | $3\times10^{-13}$ |

The tolerance is $10^{-8}$, that of the other integrated references. Reading the curve
instead of the 1.2 column misses it by $6.6\times10^{-3}$.

### `Lv`, `Mv` and `Cpv` from the property table { #lv-mv-and-cpv-from-the-property-table }

This part releases the temperature as the energy-equation part does, with `Lv`, `Mv` and
`Cpv` columns in a table on 250 to 450 K:

- $L_v(T) = 2.26\times10^6 - 2500\,(T-300)$ J/kg;
- $c_{p,v}(T) = 1800 + 2\,(T-300)$ J/(kg K);
- $M_v = 18.015$;
- a `Psat` column that is exact for that $L_v$ (Kirchhoff's law integrated through one
  atmosphere at the boiling temperature);
- `latent-heat` and `vapour-molar-mass` at 0 in the INI.

The reference is the RK4 integration of the energy-equation part with the three columns read
as written, so ICE must take $L_v(T_p)$ in the latent sink and in $B_T$, and $c_{p,v}(T_p)$ in
the gas-side heat of `ASM` and `TC`. After 0.05 s:

| `evaporation` | $T_p$ [K] (ICE) | $T_p$ [K] (RK4) | $T_p$ error | $\rho_p$ error |
|---|---|---|---|---|
| `d2-law` | 319.719179 | 319.719179 | $2.1\times10^{-11}$ | $1.8\times10^{-11}$ |
| `ASM`    | 334.667315 | 334.667315 | $1.6\times10^{-11}$ | $2.8\times10^{-12}$ |
| `TC`     | 328.147575 | 328.147575 | $1.1\times10^{-11}$ | $4.8\times10^{-12}$ |

The $T_p$ error is a fraction of the initial gap. The tolerances, $2.2\times10^{-10}$ on $T_p$
and $1.8\times10^{-10}$ on $\rho_p$, are ten times the largest of these. What fails it:

- a latent heat taken at 298.15 K in the sink instead of at $T_p$ moves $T_p$ by up to
  $4\times10^{-3}$ of the gap;
- a vapour heat capacity taken at the boiling temperature, by up to $3\times10^{-4}$.

The same `d2-law` run with `latent-heat` and `vapour-molar-mass` left out of the INI gives the
same solution file bit for bit.

!!! note "What this check can and cannot catch"
    As with [D and E](exchange.md#d-every-drag-correlation), a correlation written wrongly in
    the same way in ICE and in the reference would pass. Three things are not vulnerable to
    that: the $\dot m \propto d$ scaling, the `ASM`–`CEM` identity at $Re = 0$ and the
    isothermal identity. They follow from the published forms and would break under any
    transcription error.

## M. Drop heating and boiling

**Setup.** A uniform cloud at rest in a gas at rest, with `drag = Stokes` and
`heat-transfer = Stokes` ($Re = 0$, $Nu = 2$). The gas and the vapour are IGLOO's
($L_v$ = 2.26 MJ/kg, $M_v$ = 18), $\rho_p$ = 1 kg/m³ and $n = 1/m_0$, stepped at
`dt-max` = $h$. Every few steps the probe writes the mean temperature, $\rho_p$, $n$ and the
surface temperature `Ts_p1` of a cell. Five legs:

| Leg | Drop and gas | Models | $h$, steps | End (ICE) |
|---|---|---|---|---|
| M1 | $\rho$ 2950, $c$ 1250, $d$ 11.89 µm, $T_0$ 300 K, $T_g$ 600 K, no mass transfer | `P2T`, $k_l$ = 0.026 ($\zeta = 1$) | 1e-5 s, 1501 | $\bar T$ 599.8318 K, $T_s$ 599.8598 K |
| M2 | water, $d$ 50 µm, $T_0$ 640.002 K, $T_g = T_B$ = 600 K | `ZGR`, the adiabatic flash | 1e-6 s, 15001 | $\bar T$ 604.5757 K, $m/m_0$ 0.933822 |
| M3 | water, $d$ 50 µm, $T_0$ 490.0002 K, $T_g$ 600 K, $T_B$ 450 K | `ZGR`, fed by the gas | 1e-6 s, 15001 | $\bar T$ 454.8413 K, $m/m_0$ 0.829101 |
| M4 | water, $d$ 50 µm, $T_0$ 300 K, $T_g$ 800 K | `P2T`, $k_l$ = 0.6, CEM | 1e-5 s, 3001 | $\bar T$ 363.1647 K, $m/m_0$ 0.450683 |
| M5 | water, $d$ 20 µm, $T_0$ 640.002 K, $T_g = T_B$ = 600 K | `ZGR`, then CEM below $T_B$ | 1e-6 s, 75000 | $\bar T$ 463.8254 K (the wet bulb), $m/m_0$ 0.226316 |

### Against the same equations

The reference is ICE's own time scheme on the same equations, written in the test driver: RK2
in $(\rho_p, \rho_p e)$ with $n$ fixed, the sources from the
[models](../../theory/physics.md#conduction-inside-the-drop) at $Re = 0$, each root found by
bisection, and $T_s$ at each row's state.

| Check | Tolerance | Measured |
|---|---|---|
| $\bar T$ and $T_s$ on every row, as a fraction of the leg's gap $T_0 - T_g$ (or $T_0 - T_B$); $m/m_0$ | $10^{-9}$; $10^{-10}$ on $m/m_0$ | at most $2.3\times10^{-14}$ of the gap |
| M5 after the step that crosses $T_B$, which multiplies a difference by about $10^4$ | $10^{-5}$ | $7.9\times10^{-12}$ |
| M2's flash identity, $m\,(L_v - c(\bar T - T_B))$ constant | $10^{-10}$ | $3.7\times10^{-11}$, the vapour fed by the gas's heat |

Each leg also checks the behaviour its model must show:

- **M1** follows the Heun closed form $\bar T - T_g = (T_0 - T_g)(1 - z + z^2/2)^k$, with
  $z = h/\tau$ and $\tau = 1.2\,\tau_{ITC}$, and $T_s = (5\bar T + T_g)/6$.
- **M2 and M3** hold $T_s = T_B$ and $\bar T > T_B$.
- **M5** writes $T_s = T_B$ above $T_B$, and $T_s = \bar T$ below it, where $\bar T$ and $m$
  decrease.
- On every leg, every row's time is a multiple of $h$ (the premise of the reference) and the
  field is uniform.

### Against IGLOO

The same reference shares the equations with ICE, so a material constant mistyped in the
driver would pass. The legs are therefore also checked against IGLOO's trajectory rows of the
same drops, at $t = x/u_g$ (M4 has no IGLOO case). ICE's value at that time (the reference
carries ICE's nearest row) must lie within the sum of:

- IGLOO's own row tolerance against its oracle, from its F12.6 and E13.6 prints and its
  integrator: $\bar T$ within $1.5\times10^{-6} + 5\times10^{-7}|d\bar T/dx|$ K (on M1,
  $10^{-6}(1.01 + |d\bar T/dx|)$ K, and $T_s$ within $10^{-6}$ K), $m$ within its half-ULP
  $+\,5\times10^{-7}|dm/dx| + 10^{-8}m$;
- ICE's time-step error, bounded by $2|\bar T_h - \bar T_{h/2}|$ of the reference, and by
  $h\,|d\bar T/dt|$ times the steps of $\alpha_s$ at 5 and 25 K (the Richardson ratio is 0.67
  on M2, which crosses a step);
- ICE's difference from the reference at the row.

| Leg | IGLOO case | Compared | Measured, as a fraction of the tolerance |
|---|---|---|---|
| M1 | `temp-relax-p2t` | $\bar T$ and $T_s$ at 15 ms | 0.31, 0.58 |
| M2 | `boiling-zgr-flash` | $\bar T$ and $m$ at the outlet, 15 ms | 0.14, 0.94 (the mass is limited by the print: ICE's $5.8551549\times10^{-11}$ kg against IGLOO's 0.585515E-10) |
| M3 | `boiling-zgr-hot` | $\bar T$ and $m$ at the outlet, 15 ms | 0.51, 0.05 (IGLOO's own offset after the step of $\alpha_s$) |
| M5 | `boiling-zgr-switch` | $\bar T$ and $m$ at 60 ms (before the switch) and at 75 ms (the wet bulb) | 0.01, 0.61, 0.22, 0.64 |

A mistyped constant passes every other check and fails these:

- $c$ = 1251 for 1250 on M1: 57 and 112 times the tolerance;
- the water's $c$ = 4187 for 4186: 9 to 48 times on M2, M3 and M5.

### An unresolved switch (M5c)

M5c runs M5 at `dt-max` = 1e-4 s, which does not resolve the switch at $T_B$: the step that
crosses it removes 0.977 of the drop mass (at 61.3 ms, the reference's value to four digits).
ICE must write the
[warning](../../theory/physics.md#boiling-above-the-boiling-temperature) of an unresolved
step exactly once. M2, M3 and M5, whose fractions are 1.3e-4, 1.3e-4 and 0.012, must write
none. Without the warning, M5c fails.

### What fails it

Measured as fractions of the gap:

| Fault | Error |
|---|---|
| `P2T` ignored | 6.7e-2 on M1 (20.1 K), 3.1e-3 on M4 |
| The probe writing $\bar T$ for $T_s$ | 0.16 on M1, 0.51 on M2 |
| The boiling regime bypassed (the boiling clamp in place of `ZGR`) | 3.4 on M2 and 1.8 on M3 (468.07 K and 381.46 K); no row above $T_B$ on M5 |
| The vapour leaving at the mean temperature | 1.2e-2 on M2 and M3, 1.1e-4 on M4 |
| $\alpha_s$'s 2.33 written 2.3 | 9.9e-3 on M2 |
| The film at $k_g$ | 9.2e-4 on M3, 1.4e-2 on its mass |
| The 0.5 K floor of $\Delta T_f$ dropped | no row below $T_B$ on M5 |
| `dt-max` = 1e-4 s on M5 | a mean temperature that rises again below $T_B$ after the step's overshoot (M5c's run, which the warning flags) |
