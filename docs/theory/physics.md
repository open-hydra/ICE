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
| `JAXA4` | $\left[\left(2+0.654\,Re^{1/2}Pr^{1/3}\right)^{-1} + \dfrac{3.42\,Ma}{Re\,Pr}\right]^{-1}$ | rarefaction-dependent |
| `Ranz-Marshall` | $2 + 0.6\,Re^{1/2}Pr^{1/3}$ | 2 |
| `Kavanau-Drake` | $\dfrac{N}{1 + 3.42\,Ma\,N/(Re\,Pr)}$, $N = 2+0.459\,Re^{0.55}Pr^{0.33}$ | rarefaction-dependent |
| `NoHeat` | $0$ | 0 (no convective exchange) |

The `JAXA` names and their constants are IGLOO's, so a case that runs both solvers can name one law for both. In a
coupled run the word `Chang` is refused with a pointer to `JAXA3`; the constant 0.654 of `JAXA4` is Shimada's (2006,
eq. 50, after NASA SP-8039).

## Radiation

$$
\dot q_{\text{rad}} = \varepsilon\,\sigma\,2\pi R_p^2\,n\,\big(T_g^4 - T_p^4\big),
$$

with $\sigma = 5.67\times10^{-8}$ W m⁻² K⁻⁴. The area factor is $2\pi R_p^2$, half the sphere surface. Radiation is exchanged with the local gas temperature, not with a wall or a far-field temperature.

## Evaporation

Particles evaporation provides a mass source term. Every model returns a rate $\dot m$ **per particle**, negative while the droplet loses mass. The source routine multiplies it by the number density, so the bulk density loses $n\dot m$ per unit volume, the energy equation loses both the
enthalpy that mass carries away and the latent heat $L_v$ needed to vaporise it, and the momentum equation loses the momentum it carries. The number density has no source term at all: droplets shrink, they never disappear, and the radius follows from $\rho_p$ and $n$ as it always does.

The model, the interface and the accommodation coefficient belong to the material: each takes its phase-file tokens, or the `[ICE-Physics]` default without them (see [Materials](../user/input.md#materials)), together with its own vapour properties.

### The surface state

All five models share one surface condition. The saturation pressure is the `Psat` column
of the [property table](../user/initial-conditions.md#property-table) when it has a non-zero one, linear
between the nodes and constant beyond its ends, and otherwise Clausius-Clapeyron, anchored at
the boiling point:

$$
p_{sat}(T_p) = p_{atm}\,\exp\left[-\frac{L_v M_v}{\mathcal{R}}
  \left(\frac{1}{T_p} - \frac{1}{T_{boil}}\right)\right],
$$

with $M_v$ as the vapour molar mass, $T_{boil}$ the boiling-temperature, $L_v$ the latent-heat. The surface mole fraction is
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
plateau ($T = T_m$) and the solid ($f = 1$) it is IGLOO's enthalpy without its datum.

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

In a steady stream the nucleation front and both melting points are located within one cell. The reconstruction
works on $T_\ell$ in place of $T_p$ ([Reconstruction](numerics.md#reconstruction)): $T_\ell$ is continuous where a
content nucleates, since the recalescence keeps its energy, so the limiter sees no jump at the nucleation front, and in
every configuration measured a steady run converges there with MUSCL as at first order: `space-reconstruction =
first-order` is not needed for that. Three approximations remain. A cell where contents of
different histories meet ($0 < \chi < 1$) holds one state, decided by the majority; carrying the nucleated and the
liquid parts as two populations would lift it. At the nucleation front the liquid and the nucleated state of the
front cell can both be steady when its inflow lies in a narrow window, which puts the front one cell earlier or later.
In a time-accurate run nucleation happens at the end of a stage (above).

The Rusanov and HLLE fluxes leak $\chi$ upstream of a sharp front: $c/(2u + c)$ into the first cell with Rusanov,
nothing with HLLE while $c \le u$ and $(c - u)/(c + u)$ above, where $u$ is the flow speed and $c$ the closure's signal
speed. No cell flips while $c < 2u$ (Rusanov) or $c < 3u$ (HLLE); above that one or two cells upstream of the front
are nucleated by the leak, and the tail decays, so the front moves by that much and no further. An IG or AG stream
injected with the inlet's pseudo-pressure ($P = 10^{-6}$) is far below that bound unless its bulk density is tiny.

### Inputs

All on the material's line of the phase file, as IGLOO reads them:

| Token | Meaning | Default |
|---|---|---|
| `solidification` | `on` or `off` | `off` |
| `T-melt` | melting temperature $T_m$ [K] | 2327 |
| `h-fus` | heat of fusion [J/kg] | required, > 0 |
| `T-nuc` | nucleation temperature [K], below `T-melt` | 0.8 `T-melt` |
| `cp-solid` | specific heat of the solid [J/(kg K)] | required, > 0 |

The liquid's specific heat and the density are the material's, and must be constant: a `Cp` or `Density` column that
varies stops the run, as does a material that also evaporates. The setup prints each solidifying material's values.
For alumina the NIST-JANAF tables give $T_m$ = 2327 K, $h_{fus}$ = 1.09 MJ/kg, a liquid specific heat of
1.89 kJ/(kg K) and a solid one of 1.36 kJ/(kg K) at the melting point; the default `T-nuc` is a nominal supercooling.
