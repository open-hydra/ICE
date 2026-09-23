# Governing Equations

## Eulerian description of the condensed phase

ICE carries the dispersed phase as a continuum. Instead of tracking parcels, each cell
holds field values describing the particles inside it, and those fields obey
conservation laws of the usual form

$$
\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{U}) = \mathbf{S}.
$$

The state is a **bulk density** $\rho_p$ — the mass of condensed material per unit
volume of mixture, in kg/m³ — together with a **number density** $n$, in m⁻³. Volume
fraction is never used, and the particle radius is not carried either: it follows from
the two,

$$
R_p = \left(\frac{3}{4\pi}\,\frac{\rho_p}{n\,\rho_{al}(T_p)}\right)^{1/3},
$$

where $\rho_{al}$ is the density of the material itself, either the constant `density` from
`[ICE-Physics]` or a tabulated $\rho_{al}(T_p)$.

Both $\rho_p$ and $n$ are transported, and neither has a flux the other does not, so
the local radius is free to vary in space. Nothing in ICE enforces a size distribution;
separate sizes are represented by separate families.

## The closure problem

At the kinetic level the cloud is a distribution $f(\mathbf{x},\mathbf{v},t)$ over
particle velocity. The moments of $f$ give the fields above, and the flux of each
moment involves the next one, so the hierarchy has to be closed. ICE offers three
closures.

### MK — monokinetic

Every particle in a cell shares one velocity: $f$ is a Dirac delta in $\mathbf{v}$.
Six variables,

$$
\mathbf{P} = (\rho_p,\ u,\ v,\ w,\ T_p,\ n), \qquad
\mathbf{U} = \Big(\rho_p,\ \rho_p u,\ \rho_p v,\ \rho_p w,\
                  \rho_p\big[c_s T_p + \tfrac12 |\mathbf{u}|^2\big],\ n\Big).
$$

The flux carries each quantity at the local velocity and nothing else:

$$
\mathbf{F}\cdot\hat{\mathbf n} = u_n\,\mathbf{U}, \qquad u_n = \mathbf{u}\cdot\hat{\mathbf n}.
$$

There is no particle pressure, so the sound speed is zero and the system is only
weakly hyperbolic. Two streams that should cross cannot: the monokinetic assumption
has no way to hold two velocities in one cell, and the crossing becomes a
$\delta$-shock. This is the closure's defining limitation, not an artefact of the
discretisation, and it is what the [crossing-jets](../vv/crossing-jets.md) case
measures.

### IG — isotropic Gaussian

$f$ is a Maxwellian of uniform width: one extra scalar $P$, the velocity dispersion,
which acts exactly like a pressure. Seven variables,

$$
\mathbf{P} = (\rho_p,\ u,\ v,\ w,\ P,\ T_p,\ n),
$$

with the energy variable $E = \rho_p|\mathbf{u}|^2 + 3P$ and the flux

$$
\mathbf{F}\cdot\hat{\mathbf n} =
\Big(\rho_p u_n,\ \ \rho_p u_n\mathbf{u} + P\hat{\mathbf n},\ \
     \rho_p u_n|\mathbf{u}|^2 + 5Pu_n,\ \ \dots,\ n u_n\Big).
$$

The system is strictly hyperbolic, with sound speed

$$
a = \sqrt{3P/\rho_p}.
$$

### AG — anisotropic Gaussian

$f$ is a Maxwellian whose width differs by direction: the dispersion becomes a
symmetric tensor $P_{ij}$. Twelve variables,

$$
\mathbf{P} = (\rho_p,\ u,\ v,\ w,\ P_{11},\ P_{12},\ P_{13},\ P_{22},\ P_{23},\ P_{33},\ T_p,\ n),
$$

whose conservative form is $\rho_p u_i u_j + P_{ij}$ for the six tensor components. The
off-diagonal terms are what let the closure carry a shear: two streams crossing at an
angle show up as a $P_{12}$ that the transport equation then carries along, which is
the quantity the [Berthon Riemann problems](../vv/berthon.md) check. The sound speed
uses the mean normal stress, $a = \sqrt{3\bar P/\rho_p}$ with
$\bar P = (P_{11}+P_{22}+P_{33})/3$.

## Energy and temperature

In all three closures the particle temperature enters through the internal energy
$c_s(T_p)\,T_p$, with $c_s$ the specific heat of the condensed material — a constant or
a table lookup. Inverting the energy variable for $T_p$ is therefore implicit when the
table is used, and ICE does it in two passes: a first estimate with the constant $c_s$,
then one correction with $c_s$ evaluated at that estimate.

## Source terms

$\mathbf{S}$ collects, for every closure,

- **drag**, as a relaxation of the particle velocity towards the gas over a time
  $\tau_p$, together with the work it does on the energy;
- **convective heat exchange** with the gas, through a Nusselt number;
- **radiative exchange** with the gas, as a grey body;
- **mass transfer** between the phases, from the selected evaporation model. The
  leaving mass carries its own enthalpy and momentum out of the condensed phase, and
  the latent heat `latent-heat` is applied to it in the energy equation. The number
  density is not a source of anything, so the droplets shrink rather than vanish.
  Combustion is not implemented.

The expressions are in [Particle Physics](physics.md). All of them vanish when no gas
field is present, which is the 0-way coupled mode.

## Several families

Each `[ICE-FamilyN]` section adds one family with its own closure. The families are
advanced one after another inside each time step, over the same mesh and against the
same gas field. They do not exchange anything with each other, so $N$ families are $N$
independent systems sharing a grid — the way a polydisperse cloud is represented is by
giving each size its own family.
