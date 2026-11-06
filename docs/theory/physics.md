# Particle Physics

Source terms relate the particles with the gas phase.

## The dimensionless groups

From the local slip $\Delta\mathbf u = \mathbf u_g - \mathbf u_p$ and the particle radius $R_p$ derived from $\rho_p$ and $n$:

$$
Re = \frac{2\rho_g R_p |\Delta\mathbf u|}{\mu_g}, \qquad
Ma = \frac{|\Delta\mathbf u|}{\sqrt{\gamma R T_g}}, \qquad
Pr = \frac{\gamma}{\gamma-1}\frac{R\,\mu_g}{k_g}, \qquad
T_r = \frac{T_p}{T_g}.
$$

$Re$ is built on the **diameter**. $T_r$ and $\gamma$ are used only by the
compressible drag correlations.


## Drag

The drag coefficient sets a relaxation time,

$$
\tau_p = \frac{8\,\rho_{al}(T_p)\,R_p}{3\,\rho_g\,C_d\,|\Delta\mathbf u|},
$$

and the momentum source is that relaxation applied to the bulk density,

$$
\mathbf S_{\rho_p \mathbf u} = \rho_p \frac{\mathbf u_g - \mathbf u_p}{\tau_p},
$$

with the matching work term $\rho_p\,\mathbf u_p\cdot(\mathbf u_g-\mathbf u_p)/\tau_p$
in the energy equation. For $C_d = 24/Re$ the relaxation time reduces exactly to the
Stokes value $\rho_{al}d_p^2/(18\mu_g)$, which is the identity
[case A](../vv/verification.md) verifies.

| Name | $C_d$ | Notes |
|---|---|---|
| `Newton` | $0.45$ | Constant |
| `Stokes` | $24/Re$ | Creeping flow |
| `Schlichting` | $\frac{24}{Re}\left(1 + \frac{3}{16}Re\right)$ | Oseen correction |
| `Schiller-Naumann` | $\frac{24}{Re}\left(1+0.15\,Re^{0.687}\right)$ | Standard sub-critical fit |
| `Wen-Yu` | Schiller-Naumann for $Re \le 1000$, else $0.43$ | |
| `Putnam` | $\frac{24}{Re}\left(1+\frac{1}{6}Re^{2/3}\right)$ for $Re < 1000$, else $0.4392$ | |
| `Clift-Gauvin` | $\frac{24}{Re}\left(1+0.15Re^{0.687}+\frac{0.0175\,Re}{1+4.25\times10^{4}Re^{-1.16}}\right)$ | Valid through the drag crisis |
| `Morsi-Alexander` | $a_1 + a_2/Re + a_3/Re^2$ | Piecewise over eight $Re$ ranges; the first is exactly Stokes |
| `Carlson-Hoglund` | Wen-Yu times a rarefaction and compressibility factor in $Ma$ and $Re$ | |
| `Henderson` | Separate subsonic and supersonic fits, linearly bridged over $1 < Ma < 1.75$ | The bridge is continuous at both ends |
| `Crowe` | Wen-Yu blended towards $C_d = 2$ by $Ma$, with a $\tanh(\log_{10} Re)$ function | |
| `Hermsen` | Same structure as Crowe with a rational $Re$ function | |
| `NoDrag` | $0$ | No momentum exchange; the relaxation time is infinite |

Each correlation is a pure function of $(Re, Ma, \gamma, T_r)$; the choice travels as an
integer, so nothing mutable is shared between threads. Every one of them is checked
against an independent integration in
[case D](../vv/verification.md#d-every-drag-correlation), including at transonic and
supersonic slip for the four that have a compressible branch.

## Convective heat transfer

$$
\dot q_{\text{conv}} = 2\,Nu\;k_g\,\pi R_p\,n\,(T_g - T_p)
$$

per unit volume, which is the familiar $h A \Delta T$ with $h = Nu\,k_g/d_p$ summed over the $n$ particles in the cell.

| Name | $Nu$ | Limit as $Re \to 0$ |
|---|---|---|
| `Stokes` | $2$ | 2 |
| `JAXA1` | $2.5\,Re^{0.15} + 0.04\,Re$ | 0 |
| `JAXA2` | $2 + 0.37\,Re^{0.6}Pr^{1/3}$ | 2 |
| `JAXA3` | $2 + 0.459\,Re^{0.55}Pr^{1/3}$ (Chang) | 2 |
| `JAXA4` | $\left[\left(2+0.654\,Re^{1/2}Pr^{1/2}\right)^{-1} + \dfrac{3.42\,Ma}{Re\,Pr}\right]^{-1}$ | rarefaction-dependent |
| `Ranz-Marshall` | $2 + 0.6\,Re^{1/2}Pr^{1/3}$ | 2 |
| `Kavanau-Drake` | $\dfrac{N}{1 + 3.42\,Ma\,N/(Re\,Pr)}$, $N = 2+0.459\,Re^{0.55}Pr^{0.33}$ | rarefaction-dependent |
| `NoHeat` | $0$ | 0 (no convective exchange) |

The `JAXA` names and their constants are IGLOO's, so a case that runs both solvers can name one law for both. In a
coupled run the word `Chang` is refused with a pointer to `JAXA3`; the constant 0.654 and the exponent $Pr^{1/2}$ of
`JAXA4` are NASA SP-8039 Table II's last row (Shimada 2006 eq. 50 prints $Pr^{1/3}$).

## Radiation

$$
\dot q_{\text{rad}} = \varepsilon\,\sigma\,4\pi R_p^2\,n\,\big(T_g^4 - T_p^4\big),
$$

with $\sigma = 5.670374419\times10^{-8}$ W m⁻² K⁻⁴ (CODATA). The area factor is $4\pi R_p^2 = \pi d_p^2$, the whole sphere surface. Radiation is exchanged with the local gas temperature, not with a wall or a far-field temperature. The emissivity $\varepsilon$ is `emissivity` in `[ICE-Physics]`, one value per material: 0, the default, switches radiation off; a grey body cannot exceed 1, so a value above 1 is cut to 1 with a warning, and a negative value stops the run. A material whose property table carries an `Emissivity` column takes $\varepsilon(T_p)$ from it instead, at the particle temperature, linear between the nodes and the end value outside.

## Evaporation

Particles evaporation provides a mass source term. Every model returns a rate $\dot m$ **per particle**, negative while the droplet loses mass. The source routine multiplies it by the number density, so the bulk density loses $n\dot m$ per unit volume, the energy equation loses both the
enthalpy that mass carries away and the latent heat $L_v$ needed to vaporise it, and the momentum equation loses the momentum it carries. The number density has no source term at all: droplets shrink, they never disappear, and the radius follows from $\rho_p$ and $n$ as it always does.

The mass leaves at the surface temperature $T_s$, which every gas-side model takes. Without the two models below
$T_s$ is the mean temperature $\bar T$ the energy equation carries. Under `liquid-conduction = P2T`, and in the
boiling branch of `boiling = ZGR`, $T_s$ differs from it, and the energy equation of the cell reads

$$
\rho_p\frac{d\bar e}{dt} = n\left[\dot Q(T_s) + \dot Q_r(T_s) + \dot m(T_s)\,\big(L_v(T_s) + e(T_s) - e(\bar T)\big)\right],
\qquad \frac{d\rho_p}{dt} = n\,\dot m(T_s),
$$

with $e(T)$ the material's energy ($c\,T$, or the table's): the vapour takes its enthalpy at $T_s$, and the mean
temperature receives $e(\bar T) - e(T_s)$ per unit of evaporated mass.

When the [property table](../user/initial-conditions.md#property-table) carries the `Lv`, `Mv` and `Cpv` columns, $L_v$ and the vapour heat capacity $c_{p,v}$ are functions of the particle temperature, $L_v(T_p)$ and $c_{p,v}(T_p)$, wherever a model reads them: the latent sink, $B_T$ of the d2-law, $\varphi$ and the gas-side heat of ASM, and the vapour Lewis number and the gas-side heat of TC. $M_v$ is then the column's constant. Without the columns the `[ICE-Physics]` values hold, and a vapour heat capacity of 0 stands for the gas one.

The model, the interface and the accommodation coefficient belong to the material: each takes its phase-file tokens, or the `[ICE-Physics]` default without them (see [Materials](../user/input.md#materials)), together with its own vapour properties.

### The surface state

All five models share one surface condition. The saturation pressure is the `Psat` column
of the [property table](../user/initial-conditions.md#property-table) when it has a non-zero one, linear
between the nodes and constant beyond its ends, and otherwise Clausius-Clapeyron, anchored at
the boiling point:

$$
p_{sat}(T_p) = p_{B}\,\exp\left[-\frac{L_v M_v}{\mathcal{R}}
  \left(\frac{1}{T_p} - \frac{1}{T_{boil}}\right)\right],
$$

with $M_v$ as the vapour molar mass, $T_{boil}$ the boiling-temperature, $p_B$ the boiling-pressure at which it is
given (one atmosphere by default), $L_v$ the latent-heat. ICE evaluates it anchored at one atmosphere, through the
boiling point of the same line there, $1/T_{boil,1} = 1/T_{boil} + \ln(p_B/p_{atm})\,\mathcal{R}/(L_v M_v)$. The surface mole fraction is
$X_s = p_{sat}/p$, clamped to $1 - 10^{-12}$ once $p_{sat}$ reaches the local gas pressure — the
boiling regime. The cap keeps $B_M$ finite, about $10^{12}\,M_v/M_g$, so every model returns a finite
rate while the droplet boils (IGLOO clamps at the same value). Converting to a mass fraction $Y_s$ against the gas
molar mass gives the Spalding mass-transfer number

$$
B_M = \frac{Y_s - Y_\infty}{1 - Y_s},
$$

Evaporation stops when $Y_s$ falls to
$Y_\infty$: a gas already saturated at the droplet's own temperature takes no more vapour. The vapour diffusivity comes from the Lewis number as
$D_v = k_g/(\rho_g c_{p,g} Le)$, and $Sc = Pr\,Le$.

### The models

|  Model | $\dot m$ | Gas-side heat |
|---|---|---|
| **d2-law** | $-2\pi d\,\dfrac{k_g}{c_{p,g}}\ln(1+B_T)$, with $B_T = c_{p,g}(T_g-T_p)/L_v$ | Nusselt |
| **CEM** | $-\pi d\,\rho_g D_v\,Sh\,\ln(1+B_M)$, $Sh = 2+0.6\,Re^{1/2}Sc^{1/3}$ | Nusselt |
| **CEM-B** | CEM with $\rho_g$, $\mu_g$, $k_g$ re-evaluated at the 1/3-rule film temperature | Nusselt |
| **ASM** | CEM with $Sh^\star = 2 + (Sh_0-2)/F(B_M)$, $Sh_0$ Frössling, $F(B)=(1+B)^{0.7}\ln(1+B)/B$ | its own |
| **TC** | $-\pi d\,\rho_g D_v\,Sh\,\hat m$, $\hat m$ from a monotone transcendental solved by Newton | its own |

`d2-law` is the only one driven by the thermal Spalding number rather than the mass one; it is the classical Godsave-Spalding stagnant-film rate. `CEM` adds theRanz-Marshall convective correction, `CEM-B` the Hubbard film rule on top of it.

`ASM` (Abramzon-Sirignano) and `TC` (Tonini-Cossali) also resolve the **gas-side heat** inside their own film, Stefan flow included, and that value replaces the Nusselt one
for those two models.

### Non-equilibrium interface

This model replaces the equilibrium surface condition with the
Langmuir-Knudsen one (Miller-Harstad-Bellan model M2): the surface mole fraction is depressed below its equilibrium value by a Knudsen-layer term proportional to the evaporation rate itself,

$$
X_s = X_{s,eq} - \frac{2L_K}{d}\beta, \qquad
\beta = -\frac{\dot m\,Pr}{2\pi\mu_g d},
$$

which is implicit in $\dot m$ and is closed by a damped Picard iteration. $\alpha_e$ in $L_K$ is the accommodation coefficient. Because $L_K/d$ grows as the droplet shrinks while
$\beta$ does not, the depression deepens over a droplet's life.

### Stefan blowing

It multiplies the convective heat by the Miller-Harstad-Bellan reduction factor

$$
f_2 = \frac{b}{e^{b}-1}, \qquad b = -\frac{3}{2}\,Pr\,\tau_p\,\frac{\dot m}{m_p},
$$

which accounts for the outgoing vapour thickening the thermal film. It applies only
where the Nusselt correlation is in use, so it is ignored under `ASM` and `TC`, and
$f_2 \to 1$ as the rate vanishes.

!!! note "`d2-law` with $Nu = 2$ and blowing is exactly isothermal"
    Substituting the $d^2$-law rate into $b$ makes every property cancel and leaves
    $b = \ln(1+B_T)$ identically, so $f_2 = \ln(1+B_T)/B_T$ and the blown convective
    heat $2\pi d k_g (T_g-T_p) f_2$ equals the latent sink $-\dot m L_v$ term for
    term — for any gas, any material and any droplet temperature. That combination
    holds $T_p$ fixed and reproduces the textbook $d^2$ law exactly, which is what
    [case F](../vv/verification.md#f-evaporation) uses as its sharpest check.

### Conduction inside the drop

`liquid-conduction = ITC` (the default) gives the drop a uniform temperature. Under `P2T` the liquid temperature
is parabolic in the radius (Dombrovsky & Sazhin 2003), $T(r) = T_c + (T_s - T_c)(r/R)^2$, so that the mean
temperature is $\bar T = (2T_c + 3T_s)/5$. The heat the liquid conducts inward at the surface is the heat the
surface absorbs, which fixes $T_s$:

$$
T_s - \bar T = \frac{\dot Q_{net}(T_s)}{10\pi\,k_{eff}\,d}, \qquad
\dot Q_{net}(T) = \dot Q_g(T) + \dot Q_r(T) + \dot m(T)\,L_v(T),
$$

$\dot Q_g = \pi d\,k_g\,Nu\,(T_g - T)$ the convective heat (times $f_2$ under `evaporation-blowing = LK`, or the
model's own heat under ASM and TC), $\dot Q_r$ the radiative exchange and $\dot m$ the material's evaporation model
with the surface state at $T$. The effective conductivity $k_{eff} = \chi\,k_l$ accounts for the circulation
inside a drop in a slip flow (Abramzon & Sirignano 1989, eqs. 30, 32, 39):

$$
U_s = \frac{12.69}{32}\,|\Delta u|\,\frac{\mu_g}{\mu_l}\,\frac{Re^{1/3}}{1 + B_M}, \qquad
Pe_L = \frac{U_s\,d\,\rho_l\,c_l}{k_l}, \qquad
\chi = 1.86 + 0.86\tanh\left[2.245\log_{10}\frac{Pe_L}{30}\right],
$$

with $\chi = 1$ without slip, $B_M$ the gas-side model's Spalding number at $\bar T$, $\rho_l$ and $c_l$ at
$\bar T$ ($c_l$ the chord $e(\bar T + 1/2) - e(\bar T - 1/2)$ when the specific heat varies), $k_l$ the material's
`k-liq` or its `Conductivity` column at $\bar T$, and $\mu_l$ its `mu-liq` or its `Viscosity` column at the
surface temperature, as Abramzon & Sirignano evaluate it: with the column, $U_s$, $Pe_L$ and $k_{eff}$ follow each
trial $T_s$ of the solve. The drop's diameter is the family's one size in the cell. $T_s$ is found by
Brent's method on the bracket between $\bar T$ and $\bar T + \dot Q_{net}(\bar T)/(10\pi k_{eff} d)$ (clipped to
$[\max(1, T_{min} - 1), \max(T_g, T_{max})]$ K, the table's range), widened to $T_g$ for a heated drop or to the
low end for a cooled one, to $4\varepsilon$ in at most 100 evaluations. Every gas-side quantity then takes $T_s$:
$p_{sat}$, $Y_s$ and $B_M$, $L_v$ and $c_{p,v}$ of the table, the convective and radiative heats, and the
temperature ratio of the drag laws; the density, the diameter and the energy stay at $\bar T$.

At $Re = 0$ without mass transfer the solution is closed: $T_s = (5\bar T + \zeta T_g)/(5 + \zeta)$ with
$\zeta = Nu\,k_g/(2k_l)$, and $\bar T$ relaxes with the time constant $\tau_{ITC}(1 + \zeta/5)$. A cell below the
empty-cell density keeps $T_s = \bar T$. A surface balance without a root gives $T_s$ = NaN, and the run stops on
the invalid state after the update, naming the cell: ICE has no step rejection to retry it with.

### Boiling above the boiling temperature

`boiling = clamp` (the default) keeps the surface vapour fraction below 1 at any temperature. Under `ZGR` a drop
boils when its mean temperature reaches the boiling temperature at the gas pressure, $\bar T \ge T_B(p)$,
$p = \rho_g R_g T_g$ (Zuo, Gomes & Rutland 2000), with $T_B$ the inverse of Clausius-Clapeyron,

$$
T_B = \left(\frac{1}{T_{boil,1}} - \frac{\mathcal{R}}{L_v M_v}\ln\frac{p}{p_{atm}}\right)^{-1},
$$

or of the `Psat` column (on its linear interpolant: a column below $p$ at its last node has no boiling regime, one
at or above $p$ at its first node gives $T_{min}$). Below $T_B$ the material's model runs as under the clamp. In the
regime the surface sits at $T_B$. With $\Delta T = \bar T - T_B$ and $\Delta T_f = \max(\Delta T, 0.5\ \mathrm{K})$,
the flash fed by the superheat is

$$
G_f = \frac{\alpha_s\,\Delta T_f\,\pi d^2}{L_v}, \qquad
\alpha_s = 760\,\Delta T_f^{0.26}\ (\Delta T_f < 5),\quad 27\,\Delta T_f^{2.33}\ (< 25),\quad
13800\,\Delta T_f^{0.39}\ \mathrm{W\,m^{-2}\,K^{-1}},
$$

and for a `P2T` drop the parabolic conduction $G_f = 10\pi k_l d\,\Delta T_f/L_v$ ($k_l$ at $\bar T$). The gas heat reaches the surface
through the film blown by the vapour: with $A = c_{p,g}(T_g - T_B)/L_v$, the film conductivity
$k_f = k_g (T_f/T_g)^{0.7}$ at $T_f = T_B + (T_g - T_B)/3$ and $B = \pi k_f d\,Nu/c_{p,g}$, the heat-fed rate $G$
solves $G(1 + G_f/G) = B\ln(1 + A(1 + G_f/G))$ ($G = 0$ for $A \le 0$), by Brent's method in $\ln(G_f/G)$ on
$[-60, 60]$. With $S = G + G_f$,

$$
\dot Q_s = \frac{S\,c_{p,g}(T_g - T_B)}{e^{S/B} - 1}, \qquad \dot m = -S,
$$

$\dot Q_s$ replacing the convective heat; it is negative for a gas colder than $T_B$ and 0 at $T_g = T_B$, the
adiabatic flash, where $m\,(L_v - c\,\Delta T)$ is conserved. $L_v$ and the radiative exchange are taken at $T_B$;
the interface and blowing models do not act in the regime.

Just below $T_B$ the clamp's rate is far larger than the boiling branch's just above it: on a 20 µm water drop in a
gas at $T_B$, $d\bar T/dt = -7.7\times10^6$ K/s at $10^{-6}$ K below against $-24$ K/s at 0.1 K above. The explicit
step that carries a drop across $T_B$ overshoots the continuous transient, the more the longer the step, and a step
close to $m/|\dot m|$ at the clamp's rate ($7\times10^{-5}$ s for that drop) can leave a negative density and stop
the run. ICE adds no limiter of its own for this: with `boiling = ZGR` it watches, at every stage, the fraction
$\Delta t\,|S_1|/\rho_p$ of a cell's drop mass that the step removes, and the first time a family's exceeds 0.5 it
writes once

    [WARNING] boiling = ZGR, family <p> (material <m>, rank <r>): at iteration <n>, time <t>, a step removes <f> of a
    cell's drop mass (more than 0.5): the step does not resolve the evaporation at the boiling temperature and the
    result is not trusted; lower dt-max

and the run goes on unchanged. A run that writes it has either crossed $T_B$ with a step the clamp does not resolve,
or emptied a cell of the family within a step (a drop evaporating to its end); in the first case its solution is not
to be trusted: take a smaller `dt-max`. On that drop, integrated by ICE's scheme, the fraction is 0.012 at `dt-max` =
1e-6 s and 0.37 at 4e-5 s, where the mean temperature still ends on the wet bulb from above (463.8254 K, $m/m_0$
within 0.7 %); at 5e-5 s it is 0.57 and the mean temperature undershoots to 382 K, at 1e-4 s 0.98 and 84 K.

## Solidification

A material with `solidification = on` on its line of the phase file freezes and melts, as IGLOO's model 6 does: a
molten droplet supercools, recalesces, freezes on a plateau at the melting point and cools as a solid, and a solid or
a partly frozen droplet heated past the melting point melts there. The families of that material carry two more
variables after $n_p$, the frozen fraction $f$ and the nucleated fraction $\chi$, both carried with the mass.

### The energy

The latent heat is part of the energy per unit mass,

$$
e(T, f) = c_l\,T - f\,L(T), \qquad L(T) = h_{fus} - (c_l - c_s)(T_m - T),
$$

with $c_l$ the material's constant specific heat (the liquid's), $c_s$ = `cp-solid`, $h_{fus}$ = `h-fus` and
$T_m$ = `T-melt`: a fraction $f$ of solid and $1 - f$ of liquid at one temperature. On the liquid ($f = 0$), the
plateau ($T = T_m$) and the solid ($f = 1$) it is IGLOO's enthalpy without its datum. With a `Cp` column that varies,
$c_l T$ becomes the table's liquid energy $e_l(T)$ and $L(T) = h_{fus} - \big(e_l(T_m) - e_l(T)\big) + c_s(T_m - T)$;
$T_\ell = e_l^{-1}(e)$, and below every $c_l(T_m - T_\ell)$ reads $e_l(T_m) - e_l(T_\ell)$.

### The phase of a cell

After every update, $T_p$ and $f$ follow from the energy and $\chi$. With $T_\ell = e/c_l$, the temperature the content
would have if it were all liquid:

| Energy | $\chi$ | State |
|---|---|---|
| $T_\ell \ge T_m$ | any | liquid, $T_p = T_\ell$, $f = 0$; $\chi$ is cleared |
| $T_\ell \le T_{nuc}$ | any | nucleated; $\chi$ is set to 1 |
| $T_{nuc} < T_\ell < T_m$ | $\ge \tfrac12$ | nucleated; $\chi$ kept |
| $T_{nuc} < T_\ell < T_m$ | $< \tfrac12$ | undercooled liquid, $T_p = T_\ell$, $f = 0$; $\chi$ kept |

A nucleated content takes the lever rule at $T_m$: the plateau, $T_p = T_m$ with $f = c_l(T_m - T_\ell)/h_{fus}$,
while that is below 1, and otherwise the solid at $T_p = T_m - \big(c_l(T_m - T_\ell) - h_{fus}\big)/c_s$ with $f = 1$.

So a liquid supercools down to $T_{nuc}$ = `T-nuc` and then recalesces at constant energy, to the plateau with
$f_0 = c_l(T_m - T_{nuc})/h_{fus}$ or, when $f_0 \ge 1$, straight to the solid; it freezes on the plateau at the rate the
gas takes its heat, $h_{fus}\,\mathrm{d}f/\mathrm{d}t = -\dot q/\rho_p$, and cools as a solid. Heated, a nucleated content
warms as a solid to $T_m$, melts on the plateau at the same rate, and is liquid again at $f = 0$: nothing is solid above
$T_m$, and a melted content supercools again before it refreezes. The convective and radiative exchanges use the cell's
temperature; the recalescence exchanges nothing with the gas, and the latent heat leaves or enters it through the
plateau's exchange.

### The nucleated fraction

$\chi$ is the fraction of a cell's mass that has nucleated and not melted completely since. Nucleation sets it,
complete melting clears it, and otherwise it is only carried with the mass, so a cell follows the history of what
flows through it: liquid flowing into a cell that once nucleated replaces its content and returns it to the liquid,
and the nucleation front moves wherever the gas moves it. The frozen fraction itself could not hold that memory: it
is derived from the energy, so a nucleated cell would rebuild it every step and nucleate all the liquid that enters.
A cell where particles of different histories meet ($0 < \chi < 1$) holds one state, decided by its energy and the
majority.

### Accuracy

Along a single history, a closed cell, the formulation reproduces the particle model exactly. Its one time-step error
is the nucleation switch, which happens at the end of a stage: the plateau is shifted by at most
$\pi d_p k_g Nu\,(T_m - T_{nuc})\,\Delta t/(m h_{fus})$ in $f$. Melting has none, the heat rate being continuous through
both of its switches.

In a steady stream the nucleation front and both melting points are located within one cell. The reconstruction works
on $T_\ell$ in place of $T_p$ ([Reconstruction](numerics.md#reconstruction)): $T_\ell$ is continuous where a content
nucleates, since the recalescence keeps its energy, so $T_\ell$ has no jump at the nucleation front ($\chi$ still has
one, and its limited value decides a face in the band), and in every configuration measured a steady run converges
there with MUSCL as at first order: `space-reconstruction = first-order` is not needed for that. Three approximations
remain. A cell where contents of different histories meet ($0 < \chi < 1$) holds one state, decided by the majority;
carrying the nucleated and the liquid parts as two populations would lift it. At the nucleation front the liquid and
the nucleated state of the front cell can both be steady when its inflow lies in a narrow window, which puts the front
one cell earlier or later. In a time-accurate run nucleation happens at the end of a stage (above).

The Rusanov and HLLE fluxes leak $\chi$ upstream of a sharp front. At first order Rusanov puts $c/(2u + c)$ into the
first cell and HLLE nothing while $c \le u$ and $(c - u)/(c + u)$ above, where $u$ is the flow speed and $c$ the
closure's signal speed; with MUSCL and the van Leer limiter the Rusanov leak is about half as large
($2.5\times10^{-5}$ against $5.0\times10^{-5}$ at first order on the stream of the verification case O3). No cell
flips while $c < 2u$ (Rusanov) or $c < 3u$ (HLLE); above that one or two cells upstream of the front are nucleated by
the leak, and the tail decays, so the front moves by that much and no further. An IG or AG stream injected with the
inlet's pseudo-pressure ($P = 10^{-6}$) is far below that bound unless its bulk density is tiny.

### Inputs

All on the material's line of the phase file, as IGLOO reads them:

| Token | Meaning | Default |
|---|---|---|
| `solidification` | `on` or `off` | `off` |
| `T-melt` | melting temperature $T_m$ [K] | 2327 |
| `h-fus` | heat of fusion [J/kg] | required, > 0 |
| `T-nuc` | nucleation temperature [K], below `T-melt` | 0.8 `T-melt` |
| `cp-solid` | specific heat of the solid [J/(kg K)] | required, > 0 |

The liquid's specific heat and the density are the material's; the specific heat may vary with the temperature (a
`Cp` column, the liquid's energy then being the table's), the density must be constant: a `Density` column that varies
stops the run, as does a material that also evaporates. The setup prints each solidifying material's values.
For alumina the NIST-JANAF tables give $T_m$ = 2327 K, $h_{fus}$ = 1.09 MJ/kg, a liquid specific heat of
1.89 kJ/(kg K) and a solid one of 1.36 kJ/(kg K) at the melting point; the default `T-nuc` is a nominal supercooling.
