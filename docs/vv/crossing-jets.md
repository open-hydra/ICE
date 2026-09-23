# Crossing Jets

Two jets of particles enter a square domain at ±45° and cross at right angles. There is no
carrier gas (0-way coupling, no drag), so every particle moves in a straight line at
constant velocity. The exact solution is therefore **free streaming**: the two jets pass
through each other unchanged, and their densities add up where they overlap.

Moment methods cannot always reproduce this, because each closure assumes a shape for the
local velocity distribution. How close each closure gets to free streaming is the point of
the test:

- **MK** (monokinetic): one velocity per point. Crossing particles cannot coexist, so the
  jets collide and merge.
- **IG** (isotropic Gaussian): a velocity dispersion that is the same in every direction.
- **AG** (anisotropic Gaussian): a dispersion that can differ between directions, which
  allows two crossing streams to be represented.

## Problem setup

| Parameter | Value |
|---|---|
| Domain | $[0, 1] \times [-0.5, 0.5]$ m, one cell in $z$ |
| Grid | $100 \times 100$ cells, uniform |
| Lower inlet (face 1) | $y \in [-0.3, -0.2]$, direction $+45°$ |
| Upper inlet (face 1) | $y \in [0.2, 0.3]$, direction $-45°$ |
| Inlet condition | ATLAS `402`: mass flux $g = 0.7071$ kg/(m² s), speed $10$ m/s, $T_p = 300$ K, $r_p = 1$ µm |
| Inlet density | $\rho = g / \lvert \mathbf{u}\cdot\mathbf{n} \rvert = 0.1$ kg/m³ |
| Rest of face 1 | Symmetry (`300`) |
| Faces 2–4 | Extrapolation (`400`) |
| Coupling | 0-way (no gas) |
| Time | Local time stepping, 5000 iterations, CFL 0.8 |
| Space reconstruction | MUSCL with shock detector, Van Leer limiter |

In the exact solution each jet is a 45° band of density 0.1 kg/m³. The bands overlap in a
diamond centred at $(0.25, 0)$, where the density is 0.2 kg/m³.

## Results

<figure>
  {% include "vv/images/doisneau-fields.svg" %}
</figure>

<figure>
  {% include "vv/images/doisneau-profiles.svg" %}
</figure>

| Model | Peak density [kg/m³] | Behaviour at the crossing |
|---|---|---|
| Exact | 0.200 | Jets pass through each other |
| MK | 2.00 | Jets merge into a single jet along $y = 0$; mass concentrates in the axis cells (delta shock) |
| IG | 0.131 | Jets merge into a single spreading jet |
| AG | 0.100 | Jets cross and continue as two separate, diffused jets |

- **MK** carries exactly the nominal inlet mass flux (0.1414 kg/(m s) through every
  vertical line before the jets reach the outer boundaries). Its merged jet, a delta shock,
  is the expected monokinetic answer to crossing streams.
- **IG** cannot hold two velocities at the same point either, but its dispersion spreads
  the merged jet instead of concentrating it.
- **AG** is the only closure that lets the jets cross. Downstream they are much wider and
  weaker than the exact bands.

The cases are also part of the regression suite (`Doisneau/MK`, `Doisneau/IG`,
`Doisneau/AG`), which checks the density field against a stored reference.
