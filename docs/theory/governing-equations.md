# Governing Equations

## Eulerian description of the condensed phase

ICE models the dispersed (particulate) phase using an **Eulerian** framework: instead of tracking individual particles, the condensed phase is treated as a continuum characterised by local volume fraction, velocity, and energy fields.

The **monokinetic (MK) closure** assumes that all particles within a group share a single velocity at each point in space — i.e., the velocity distribution function is a Dirac delta. This yields a closed set of hyperbolic conservation laws.

## Conservative form

For each particle group $g$, the governing equations read:

$$
\frac{\partial \mathbf{U}_g}{\partial t} + \nabla \cdot \mathbf{F}_g = \mathbf{S}_g
$$

with conservative variables

$$
\mathbf{U}_g = \begin{pmatrix}
  \alpha_g \rho_g \\
  \alpha_g \rho_g \mathbf{u}_g \\
  \alpha_g \rho_g e_g \\
  \alpha_g
\end{pmatrix}
$$

where $\alpha_g$ is the volume fraction, $\rho_g$ the particle material density (constant), $\mathbf{u}_g$ the particle velocity, and $e_g$ the specific total energy.

## Convective flux

The convective flux tensor is

$$
\mathbf{F}_g = \mathbf{u}_g \otimes \mathbf{U}_g
$$

Under the MK closure there is no pressure contribution from the dispersed phase (pressureless gas model).

## Source terms

The source vector $\mathbf{S}_g$ collects:

- **Drag** — momentum exchange with the carrier gas.
- **Heat transfer** — energy exchange with the gas phase and radiative losses.
- **Phase change** — mass transfer between particle groups and the gas due to vaporisation or combustion.

See [Particle Physics](physics.md) for model details.

## Multi-group system

For $N_\text{grp}$ groups, $N_\text{grp}$ independent systems of the form above are solved simultaneously, coupled only through the shared gas-phase fields provided by MOSE.

---
