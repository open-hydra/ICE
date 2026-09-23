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

`drag` in `[ICE-Scheme]` selects the correlation, globally for all families:

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
over the $n$ particles in the cell. `heat` in `[ICE-Scheme]` selects $Nu$:

| Name | $Nu$ | Limit as $Re \to 0$ |
|---|---|---|
| `Stokes` | $2$ | 2 |
| `JAXA1` | $2.5\,Re^{0.15} + 0.04\,Re$ | 0 |
| `JAXA2` | $2 + 0.37\,Re^{0.6}Pr^{1/3}$ | 2 |
| `JAXA3` | $\left[\left(2+0.645\,Re^{1/2}Pr^{1/3}\right)^{-1} + \dfrac{3.42\,Ma}{Re\,Pr}\right]^{-1}$ | rarefaction-dependent |
| `Chang` | $2 + 0.459\,Re^{0.55}Pr^{1/3}$ | 2 |
| `Ranz-Marshall` | $2 + 0.6\,Re^{1/2}Pr^{1/3}$ | 2 |
| `Kavanau-Drake` | $\dfrac{N}{1 + 3.42\,Ma\,N/(Re\,Pr)}$, $N = 2+0.459\,Re^{0.55}Pr^{0.33}$ | rarefaction-dependent |

With `heat = Stokes` the exchange is a linear relaxation of $T_p$ towards $T_g$ over
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

with $\sigma = 5.67\times10^{-8}$ W m⁻² K⁻⁴ and $\varepsilon$ from `emiss` in
`[ICE-Physics]`. The area factor is $2\pi R_p^2$, half the sphere surface. Radiation
is exchanged with the local gas temperature, not with a wall or a far-field
temperature, and `emiss = 0` switches it off.

## Material properties

$\rho_{al}$ and $c_s$ are needed at every conversion between conservative and primitive
variables, so their temperature dependence matters even without any source term. They
come either from `[ICE-Physics]` `rho` and `cs` as constants, or from a table in
`INPUT/part-properties.dat` indexed by integer temperature. ICE prints which of the two
it is using at startup — an absent table is a legitimate configuration, but a
*silently* absent one is indistinguishable from a mis-named file.

## Mass transfer

The closures carry a mass-transfer rate, and the latent heat `lv` multiplies it in the
energy equation, but the source routine sets the rate to zero: vaporisation and
combustion are not evaluated. `lv` and `q` in `[ICE-Physics]` consequently have no
effect on a run.
