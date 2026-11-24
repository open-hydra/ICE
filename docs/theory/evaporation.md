# Evaporation

A drop that evaporates loses mass to the gas as vapour, and its surface cools by the latent heat. The symbols are those
of [Drag & Heat Transfer](drag-heat.md#notation).

!!! info "Switching it on"
    `evaporation = <model>` in `[ICE-Physics]` or on a material's line: see
    [Materials & Models](../user/models.md#evaporation).

## What evaporation changes in a cell

Every model returns a rate $\dot m$ **per drop**, negative while the drop loses mass. Per unit volume, the cell's drops
then lose:

- the mass $n\,|\dot m|$;
- the momentum that mass carries;
- the energy that mass carries, plus the latent heat $L_v$ needed to vaporise it.

The number of drops does not change: drops shrink, they never disappear. Their radius always follows from $\rho_p$ and
$n$.

## The surface temperature

The vapour leaves at the surface temperature $T_s$, which every gas-side model takes. Usually $T_s$ is the mean
temperature $\bar T$. It differs from it in two cases: inside a drop with a temperature profile
([conduction inside the drop](#conduction-inside-the-drop)), and in a boiling drop ([boiling](#boiling)). Then the
energy of the cell reads

$$
\rho_p\frac{d\bar e}{dt} = n\left[\dot Q(T_s) + \dot Q_r(T_s) + \dot m(T_s)\,\big(L_v(T_s) + e(T_s) - e(\bar T)\big)\right],
\qquad \frac{d\rho_p}{dt} = n\,\dot m(T_s),
$$

with $e(T)$ the material's energy ($c\,T$, or the property table's). The vapour takes its enthalpy at $T_s$; the mean
temperature receives $e(\bar T) - e(T_s)$ per unit of evaporated mass.

## Saturation at the surface

All five models share one surface condition. The saturation pressure is the `Psat` column of the
[property table](../user/properties.md#the-property-table) when there is one (linear between nodes, constant beyond
the ends). Otherwise it is the Clausius-Clapeyron curve through the boiling point:

$$
p_{sat}(T_p) = p_{B}\,\exp\left[-\frac{L_v M_v}{\mathcal{R}}
  \left(\frac{1}{T_p} - \frac{1}{T_{boil}}\right)\right],
$$

with $M_v$ the vapour molar mass, $T_{boil}$ the boiling temperature and $p_B$ the pressure it is given at (one
atmosphere by default).

From the saturation pressure:

1. The surface mole fraction is $X_s = p_{sat}/p$. It is capped at $1 - 10^{-12}$ once $p_{sat}$ reaches the gas
   pressure, so that every model still returns a finite rate while the drop boils.
2. $X_s$ becomes a mass fraction $Y_s$, through the gas molar mass.
3. The Spalding mass-transfer number is

    $$
    B_M = \frac{Y_s - Y_\infty}{1 - Y_s}.
    $$

Evaporation stops when $Y_s$ falls to $Y_\infty$: a gas already saturated at the drop's temperature takes no more
vapour. The vapour diffusivity is $D_v = k_g/(\rho_g c_{p,g} Le)$, and $Sc = Pr\,Le$.

With the `Lv`, `Mv` and `Cpv` columns in the property table, $L_v$ and the vapour specific heat $c_{p,v}$ follow the
drop temperature wherever a model reads them, and $M_v$ is the column's value.

## Evaporation models

|  Model | $\dot m$ | Gas-side heat |
|---|---|---|
| **d2-law** | $-2\pi d\,\dfrac{k_g}{c_{p,g}}\ln(1+B_T)$, with $B_T = c_{p,g}(T_g-T_p)/L_v$ | Nusselt |
| **CEM** | $-\pi d\,\rho_g D_v\,Sh\,\ln(1+B_M)$, $Sh = 2+0.6\,Re^{1/2}Sc^{1/3}$ | Nusselt |
| **CEM-B** | CEM with $\rho_g$, $\mu_g$, $k_g$ re-evaluated at the 1/3-rule film temperature | Nusselt |
| **ASM** | CEM with $Sh^\star = 2 + (Sh_0-2)/F(B_M)$, $Sh_0$ Frössling, $F(B)=(1+B)^{0.7}\ln(1+B)/B$ | its own |
| **TC** | $-\pi d\,\rho_g D_v\,Sh\,\hat m$, $\hat m$ from a monotone transcendental solved by Newton | its own |

- `d2-law` is the classical Godsave-Spalding stagnant-film rate, the only model driven by the thermal Spalding number.
- `CEM` adds the Ranz-Marshall convective correction; `CEM-B` adds the Hubbard film rule on top.
- `ASM` (Abramzon-Sirignano) and `TC` (Tonini-Cossali) compute the gas-side heat inside their own film, Stefan flow
  included. That heat replaces the Nusselt one.

## Non-equilibrium interface

The Langmuir-Knudsen interface (Miller-Harstad-Bellan model M2) replaces the equilibrium surface condition. A Knudsen
layer lowers the surface mole fraction in proportion to the evaporation rate itself:

$$
X_s = X_{s,eq} - \frac{2L_K}{d}\beta, \qquad
\beta = -\frac{\dot m\,Pr}{2\pi\mu_g d}.
$$

$\alpha_e$ in $L_K$ is the accommodation coefficient. The relation is implicit in $\dot m$; a damped fixed-point
iteration solves it. Since $L_K/d$ grows as the drop shrinks, the effect grows over a drop's life.

## Stefan blowing

The outgoing vapour thickens the thermal film and reduces the convective heat by the Miller-Harstad-Bellan factor

$$
f_2 = \frac{b}{e^{b}-1}, \qquad b = -\frac{3}{2}\,Pr\,\tau_p\,\frac{\dot m}{m_p}.
$$

It applies only where the Nusselt correlation is in use, so not under `ASM` and `TC`. $f_2 \to 1$ as the rate
vanishes.

!!! note "`d2-law` with $Nu = 2$ and blowing keeps the drop temperature fixed"
    Substituting the $d^2$-law rate into $b$ cancels every property and leaves $b = \ln(1+B_T)$, so
    $f_2 = \ln(1+B_T)/B_T$. The blown convective heat $2\pi d k_g (T_g-T_p) f_2$ then equals the latent sink
    $-\dot m L_v$ exactly, for any gas, material and drop temperature. The drop temperature stays fixed and the
    textbook $d^2$ law is reproduced exactly:
    [case F](../vv/verification/phase-change.md#an-exact-check-with-nothing-frozen) uses this as its sharpest check.

## Conduction inside the drop

By default (`ITC`) a drop has one uniform temperature. Under `P2T` its temperature is parabolic in the radius
(Dombrovsky & Sazhin 2003):

$$
T(r) = T_c + (T_s - T_c)\,(r/R)^2, \qquad \bar T = \frac{2T_c + 3T_s}{5}.
$$

**Surface temperature.** The heat the liquid conducts inward at the surface is the heat the surface absorbs. That
balance fixes $T_s$:

$$
T_s - \bar T = \frac{\dot Q_{net}(T_s)}{10\pi\,k_{eff}\,d}, \qquad
\dot Q_{net}(T) = \dot Q_g(T) + \dot Q_r(T) + \dot m(T)\,L_v(T),
$$

where $\dot Q_g = \pi d\,k_g\,Nu\,(T_g - T)$ is the convective heat (times $f_2$ with Stefan blowing, or the model's
own heat under ASM and TC), $\dot Q_r$ the radiative exchange, and $\dot m$ the evaporation rate with the surface
state at $T$.

**Circulation inside the drop.** A drop in a slip flow circulates inside, which raises its effective conductivity to
$k_{eff} = \chi\,k_l$ (Abramzon & Sirignano 1989, eqs. 30, 32, 39):

$$
U_s = \frac{12.69}{32}\,|\Delta u|\,\frac{\mu_g}{\mu_l}\,\frac{Re^{1/3}}{1 + B_M}, \qquad
Pe_L = \frac{U_s\,d\,\rho_l\,c_l}{k_l}, \qquad
\chi = 1.86 + 0.86\tanh\left[2.245\log_{10}\frac{Pe_L}{30}\right],
$$

with $\chi = 1$ without slip. $k_l$ and $\mu_l$ are the liquid conductivity and viscosity (`k-liq`, `mu-liq`, or the
`Conductivity` and `Viscosity` columns). $\rho_l$, $c_l$ and $k_l$ are taken at $\bar T$, and $\mu_l$ at the surface
temperature, as Abramzon & Sirignano do.

**How $T_s$ is found.** Brent's method solves the balance, to $4\varepsilon$ in at most 100 evaluations. The
bracket starts between $\bar T$ and $\bar T + \dot Q_{net}(\bar T)/(10\pi k_{eff} d)$ and widens to $T_g$ for a heated
drop, or to the low end for a cooled one. Every gas-side quantity then takes $T_s$: $p_{sat}$, $Y_s$, $B_M$, $L_v$ and
$c_{p,v}$, the convective and radiative heats, and the temperature ratio of the drag laws. The density, the diameter
and the energy stay at $\bar T$.

**A closed form.** At $Re = 0$ without mass transfer, $T_s = (5\bar T + \zeta T_g)/(5 + \zeta)$ with
$\zeta = Nu\,k_g/(2k_l)$, and $\bar T$ relaxes with the time constant $\tau_{ITC}(1 + \zeta/5)$.

**Edge cases.**

- A cell below the empty-cell density keeps $T_s = \bar T$.
- A surface balance with no root gives $T_s$ = NaN. The run then stops on the invalid state after the update and
  names the cell: ICE has no step rejection to retry it with.

## Boiling

By default (`clamp`) the surface vapour fraction stays below 1 at any temperature. Under `ZGR` (Zuo, Gomes & Rutland
2000) a drop **boils** once its mean temperature reaches the boiling temperature at the gas pressure,
$\bar T \ge T_B(p)$ with $p = \rho_g R_g T_g$. $T_B$ inverts the Clausius-Clapeyron curve,

$$
T_B = \left(\frac{1}{T_{boil,1}} - \frac{\mathcal{R}}{L_v M_v}\ln\frac{p}{p_{atm}}\right)^{-1},
$$

with $T_{boil,1}$ the boiling temperature at one atmosphere, on the same curve:
$1/T_{boil,1} = 1/T_{boil} + \ln(p_B/p_{atm})\,\mathcal{R}/(L_v M_v)$. With a `Psat` column, $T_B$ inverts the column
instead: a column below $p$ at its last node has no boiling regime, and one at or above $p$ at its first node gives
$T_{min}$.

Below $T_B$ the material's evaporation model runs as usual.

While the drop boils, its surface sits at $T_B$ and two mechanisms vaporise it.

**The flash, fed by the superheat.** With $\Delta T = \bar T - T_B$ and $\Delta T_f = \max(\Delta T, 0.5\ \mathrm{K})$,

$$
G_f = \frac{\alpha_s\,\Delta T_f\,\pi d^2}{L_v}, \qquad
\alpha_s = 760\,\Delta T_f^{0.26}\ (\Delta T_f < 5),\quad 27\,\Delta T_f^{2.33}\ (< 25),\quad
13800\,\Delta T_f^{0.39}\ \mathrm{W\,m^{-2}\,K^{-1}}.
$$

A `P2T` drop takes its parabolic conduction instead, $G_f = 10\pi k_l d\,\Delta T_f/L_v$.

**The heat from the gas**, through the film the vapour blows. With $A = c_{p,g}(T_g - T_B)/L_v$, the film
conductivity $k_f = k_g (T_f/T_g)^{0.7}$ at $T_f = T_B + (T_g - T_B)/3$, and $B = \pi k_f d\,Nu/c_{p,g}$, the
heat-fed rate $G$ solves

$$
G\,(1 + G_f/G) = B\ln\big(1 + A(1 + G_f/G)\big)
$$

($G = 0$ for $A \le 0$). With $S = G + G_f$ the drop then loses mass and receives heat as

$$
\dot m = -S, \qquad \dot Q_s = \frac{S\,c_{p,g}(T_g - T_B)}{e^{S/B} - 1},
$$

$\dot Q_s$ replacing the convective heat. It is negative for a gas colder than $T_B$ and 0 at $T_g = T_B$, the
adiabatic flash. $L_v$ and the radiative exchange are taken at $T_B$; the interface and blowing models do not act while
the drop boils.

!!! warning "The time step near the boiling temperature"
    Just below $T_B$ the evaporation rate is far larger than just above it. On a 20 µm water drop in a gas at $T_B$,
    $d\bar T/dt = -7.7\times10^6$ K/s at $10^{-6}$ K below, against $-24$ K/s at 0.1 K above. A step that carries a
    drop across $T_B$ overshoots, the more the longer the step, and can even leave a negative density.

    Under `ZGR`, ICE watches the fraction $\Delta t\,|S_1|/\rho_p$ of a cell's drop mass that one step removes. The
    first time it exceeds 0.5 for a family, ICE warns once:

        [WARNING] boiling = ZGR, family <p> (material <m>, rank <r>): at iteration <n>, time <t>, a step removes <f> of a
        cell's drop mass (more than 0.5): the step does not resolve the evaporation at the boiling temperature and the
        result is not trusted; lower dt-max

    The run goes on. The warning has two causes: a drop crossing $T_B$ with too long a step, or a cell whose drops
    evaporate completely within one step. In the first case the solution is not to be trusted: lower `dt-max`. On the drop above, the fraction is 0.012 at `dt-max` = 1e-6 s and 0.37 at 4e-5 s, where the result is
    still right; at 5e-5 s it is 0.57 and the mean temperature undershoots to 382 K.
