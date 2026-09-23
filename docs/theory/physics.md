# Particle Physics

Every source term below is evaluated only when a gas field is present. Without
`INPUT/gas.tec` the run is 0-way coupled and the cloud is transported with no exchange.
The gas supplies, per cell, its density, velocity, temperature, $R$, $\gamma$, thermal
conductivity and viscosity; ICE never modifies it.

## The dimensionless groups

From the local slip $\Delta\mathbf u = \mathbf u_g - \mathbf u_p$ and the particle
radius $R_p$ derived from $\rho_p$ and $n$:

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

`drag` in `[ICE-Physics]` selects the correlation, globally for all families:

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

Each correlation is a pure function of $(Re, Ma, \gamma, T_r)$; the choice travels as an
integer, so nothing mutable is shared between threads. Every one of them is checked
against an independent integration in
[case D](../vv/verification.md#d-every-drag-correlation), including at transonic and
supersonic slip for the four that have a compressible branch.

!!! note "`Chang` is not a separate model"
    The Chang correlation, $\frac{24}{Re}(1+0.15Re^{0.687}) + \frac{0.42}{1+4.25\times10^4 Re^{-1.16}}$,
    is Clift-Gauvin rewritten: $24 \times 0.0175 = 0.42$. Asking for it stops the run
    with a message pointing at `Clift-Gauvin`.

## Convective heat transfer

$$
\dot q_{\text{conv}} = 2\,Nu\;k_g\,\pi R_p\,n\,(T_g - T_p)
$$

per unit volume, which is the familiar $h A \Delta T$ with $h = Nu\,k_g/d_p$ summed
over the $n$ particles in the cell. `heat-transfer` in `[ICE-Physics]` selects $Nu$:

| Name | $Nu$ | Limit as $Re \to 0$ |
|---|---|---|
| `Stokes` | $2$ | 2 |
| `JAXA1` | $2.5\,Re^{0.15} + 0.04\,Re$ | 0 |
| `JAXA2` | $2 + 0.37\,Re^{0.6}Pr^{1/3}$ | 2 |
| `JAXA3` | $\left[\left(2+0.645\,Re^{1/2}Pr^{1/3}\right)^{-1} + \dfrac{3.42\,Ma}{Re\,Pr}\right]^{-1}$ | rarefaction-dependent |
| `Chang` | $2 + 0.459\,Re^{0.55}Pr^{1/3}$ | 2 |
| `Ranz-Marshall` | $2 + 0.6\,Re^{1/2}Pr^{1/3}$ | 2 |
| `Kavanau-Drake` | $\dfrac{N}{1 + 3.42\,Ma\,N/(Re\,Pr)}$, $N = 2+0.459\,Re^{0.55}Pr^{0.33}$ | rarefaction-dependent |

With `heat-transfer = Stokes` the exchange is a linear relaxation of $T_p$ towards $T_g$ over
$\tau_T = \rho_{al}c_s d_p^2/(12 k_g)$, which is what
[case B](../vv/verification.md#b-thermal-relaxation) verifies; the other six are
compared against an independent integration in
[case E](../vv/verification.md#e-every-nusselt-correlation).

Unlike the drag models, `Chang` and `Kavanau-Drake` are genuinely different here: they
differ in the Prandtl exponent and in the rarefaction denominator.

## Radiation

$$
\dot q_{\text{rad}} = \varepsilon\,\sigma\,2\pi R_p^2\,n\,\big(T_g^4 - T_p^4\big),
$$

with $\sigma = 5.67\times10^{-8}$ W m⁻² K⁻⁴ and $\varepsilon$ from `emissivity` in
`[ICE-Physics]`. The area factor is $2\pi R_p^2$, half the sphere surface. Radiation
is exchanged with the local gas temperature, not with a wall or a far-field
temperature, and `emissivity = 0` switches it off.

## Material properties

$\rho_{al}$ and $c_s$ are needed at every conversion between conservative and primitive
variables, so their temperature dependence matters even without any source term. They
come either from `[ICE-Physics]` `density` and `specific-heat` as constants, or from a table in
`INPUT/part-properties.dat` indexed by integer temperature. ICE prints which of the two
it is using at startup — an absent table is a legitimate configuration, but a
*silently* absent one is indistinguishable from a mis-named file.

## Evaporation

`evaporation` in `[ICE-Physics]` selects a mass-transfer model, globally for all
families. Every model returns a rate $\dot m$ **per particle**, negative while the
droplet loses mass. The source routine multiplies it by the number density, so the
bulk density loses $n\dot m$ per unit volume, the energy equation loses both the
enthalpy that mass carries away and the latent heat $L_v$ needed to vaporise it, and
the momentum equation loses the momentum it carries. The number density has no source
term at all: droplets shrink, they never disappear, and the radius follows from
$\rho_p$ and $n$ as it always does.

### The surface state

All five models share one surface condition. The saturation pressure is
Clausius-Clapeyron, anchored at the boiling point rather than at a tabulated curve:

$$
p_{sat}(T_p) = p_{atm}\,\exp\left[-\frac{L_v M_v}{\mathcal{R}}
  \left(\frac{1}{T_p} - \frac{1}{T_{boil}}\right)\right],
$$

with $M_v$ from `vapour-molar-mass`, $T_{boil}$ from `boiling-temperature`, $L_v$ from
`latent-heat` and $\mathcal R = 8314.46$ J kmol⁻¹ K⁻¹. The surface mole fraction is
$X_s = p_{sat}/p$, clamped to 1 once $p_{sat}$ reaches the local gas pressure — the
boiling regime. Converting to a mass fraction $Y_s$ against the gas molar mass gives
the Spalding mass-transfer number

$$
B_M = \frac{Y_s - Y_\infty}{1 - Y_s},
$$

with $Y_\infty$ from `vapour-mass-fraction`. Evaporation stops when $Y_s$ falls to
$Y_\infty$: a gas already saturated at the droplet's own temperature takes no more
vapour. The vapour diffusivity comes from `lewis-number` as
$D_v = k_g/(\rho_g c_{p,g} Le)$, and $Sc = Pr\,Le$.

### The models

| `evaporation` | $\dot m$ | Gas-side heat |
|---|---|---|
| `none` | $0$ | Nusselt |
| `d2-law` | $-2\pi d\,\dfrac{k_g}{c_{p,g}}\ln(1+B_T)$, with $B_T = c_{p,g}(T_g-T_p)/L_v$ | Nusselt |
| `CEM` | $-\pi d\,\rho_g D_v\,Sh\,\ln(1+B_M)$, $Sh = 2+0.6\,Re^{1/2}Sc^{1/3}$ | Nusselt |
| `CEM-B` | CEM with $\rho_g$, $\mu_g$, $k_g$ re-evaluated at the 1/3-rule film temperature | Nusselt |
| `ASM` | CEM with $Sh^\star = 2 + (Sh_0-2)/F(B_M)$, $Sh_0$ Frössling, $F(B)=(1+B)^{0.7}\ln(1+B)/B$ | its own |
| `TC` | $-\pi d\,\rho_g D_v\,Sh\,\hat m$, $\hat m$ from a monotone transcendental solved by Newton | its own |

`d2-law` is the only one driven by the thermal Spalding number rather than the mass
one; it is the classical Godsave-Spalding stagnant-film rate. `CEM` adds the
Ranz-Marshall convective correction, `CEM-B` the Hubbard film rule on top of it.

`ASM` (Abramzon-Sirignano) and `TC` (Tonini-Cossali) also resolve the **gas-side heat**
inside their own film, Stefan flow included, and that value replaces the Nusselt one
for those two models rather than correcting it. `heat-transfer` therefore has no effect
under `ASM` or `TC`. Neither includes the latent heat: that is added once, by the
closure, so it cannot be double-counted.

### Non-equilibrium interface

`evaporation-interface = LK` replaces the equilibrium surface condition with the
Langmuir-Knudsen one (Miller-Harstad-Bellan model M2): the surface mole fraction is
depressed below its equilibrium value by a Knudsen-layer term proportional to the
evaporation rate itself,

$$
X_s = X_{s,eq} - \frac{2L_K}{d}\beta, \qquad
\beta = -\frac{\dot m\,Pr}{2\pi\mu_g d},
$$

which is implicit in $\dot m$ and is closed by a damped Picard iteration.
`evaporation-coefficient` is the accommodation coefficient $\alpha_e$ in $L_K$; the
default 1 is full accommodation. Because $L_K/d$ grows as the droplet shrinks while
$\beta$ does not, the depression deepens over a droplet's life, and a run with `LK`
no longer follows an exact $d^2$ law. The default `VLE` is the equilibrium condition
and costs nothing.

### Stefan blowing

`evaporation-blowing = LK` multiplies the convective heat by the
Miller-Harstad-Bellan reduction factor

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
