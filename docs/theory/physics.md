# Particle Physics

## Drag

The momentum exchange between particles and the carrier gas is modelled via a drag force per unit volume:

$$
\mathbf{f}_{drag} = \frac{3}{4} \frac{\alpha_p \rho_g}{d_p} C_D(Re_p, Ma) |\mathbf{u}_g - \mathbf{u}_p| (\mathbf{u}_g - \mathbf{u}_p)
$$

where $d_p$ is the particle diameter, $Re_p = \rho_g |\mathbf{u}_g - \mathbf{u}_p| d_p / \mu_g$ the particle Reynolds number, and $C_D$ the drag coefficient.

Available drag models are selected via the `drag` keyword in the `[part-N]` section of `input.ini`.

## Heat transfer

Convective heat exchange between particle and gas is modelled via:

$$
\dot{q}_{conv} = \frac{6 \alpha_p \lambda_g}{d_p^2}\, Nu(Re_p, Pr, Ma) \,(T_g - T_p)
$$

### Nusselt correlations

| Model | Expression | Range |
|-------|-----------|-------|
| `Stokes` | $Nu = 2$ | $Re_p \to 0$ |
| `Ranz-Marshall` | $Nu = 2 + 0.6\,Re_p^{1/2}\,Pr^{1/3}$ | Moderate $Re_p$ |

The model is selected via the `heat` keyword in `input.ini`.

## Radiation

Thermal radiation is modelled via the Stefan–Boltzmann law:

$$
\dot{q}_{rad} = \varepsilon \sigma (T_{ref}^4 - T_p^4)
$$

where $\varepsilon$ is the particle emissivity (`emiss` in `Parameters_m`) and $\sigma = 5.67 \times 10^{-8}$ W m$^{-2}$ K$^{-4}$ the Stefan–Boltzmann constant.

## Phase change

Particle vaporisation and combustion are modelled as mass source terms. The latent heat $L_v$ and heat of combustion $q_p$ for the particle material are specified in the `[part-N]` section.

---
