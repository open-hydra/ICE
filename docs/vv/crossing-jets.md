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
| Upper inlet (face 1) | $y \in [0.2, 0.3]$ (MK: $[0.21, 0.31]$), direction $-45°$ |
| Inlet condition | ATLAS `402`: mass flux $g = 0.7071$ kg/(m² s), speed $10$ m/s, $T_p = 300$ K, $r_p = 1$ µm |
| Inlet density | $\rho = g / \lvert \mathbf{u}\cdot\mathbf{n} \rvert = 0.1$ kg/m³ |
| Rest of face 1 | Symmetry (`300`) |
| Faces 2–4 | Extrapolation (`400`) |
| Coupling | 0-way (no gas) |
| Time | Local time stepping, 5000 iterations, CFL 0.8 |
| Space reconstruction | MUSCL, Van Leer limiter, no shock detector |

In the exact solution each jet is a 45° band of density 0.1 kg/m³. The bands overlap in a
diamond centred at $(0.25, 0)$, where the density is 0.2 kg/m³.

MK's upper inlet sits one cell higher: with both inlets symmetric about $y = 0$, its two
pressureless jets collide on the face between rows 50 and 51, and which row takes the merged
delta is decided at round-off (it differed between compilers). One cell higher, the jets are
mirror images about the centre of row 51 ($y = 0.005$), which takes the delta.

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
| MK | 2.00 | Jets merge into a single jet along $y = 0.005$; mass concentrates in the cells of row 51 (delta shock) |
| IG | 0.346 | Jets merge where they meet, then spread as a single jet |
| AG | 0.234 | Jets cross and continue as two separate, diffused jets; the peak lies beside the crossing |

- **MK** merges the jets into one along $y = 0.005$, midway between its inlets, and piles
  their mass into the cells of row 51. This delta shock is the expected monokinetic answer
  to crossing streams.
- **IG** cannot hold two velocities at the same point either. The jets merge where they
  meet, with a peak above the exact diamond's, and its dispersion then spreads the merged
  jet instead of concentrating it.
- **AG** is the only closure that lets the jets cross. Downstream they are much wider and
  weaker than the exact bands.

The three closures run without the [shock detector](../theory/numerics.md#shock-detector).
For MK and AG this is a trade. Without it MK stays nearer its exact solution (an RMS distance
of $5.6 \times 10^{-3}$ kg/m³ against $1.4 \times 10^{-2}$ with the detector; at $x = 0.15$,
upstream of where the jets meet, row 51 holds $8 \times 10^{-4}$ kg/m³ against 0.069), and so
does AG ($2.6 \times 10^{-2}$ against $2.8 \times 10^{-2}$; the crossing, 0.2 kg/m³ in the
exact solution, holds 0.06 on average against 0.04). But neither converges: MK's residual
stalls at $1.6 \times 10^{-5}$ and its jets are $1.2 \times 10^{-2}$ from their mirror image at
the outer edges of the inlets, and AG's residual cycles between $4 \times 10^{-5}$ and
$9 \times 10^{-5}$. AG also gathers mass on either side of the crossing, up to 0.234 kg/m³
where the exact solution has at most 0.15 (the field's maximum is 0.115 with the detector),
and clips its pressure tensor at the floor in 22 nearly empty cells ($1$ to
$3 \times 10^{-3}$ kg/m³, where the exact solution is empty): 16 between the jets next to the
inlet face, 6 outside the jets beside the crossing. Both tests therefore run their 5000
iterations; with the detector, MK reaches the default residual threshold at iteration 2493
and AG at iteration 6907. IG converges without it (a residual of $1.4 \times 10^{-6}$ at 5000
iterations), and its chimera and multi-block versions then match it to round-off
([Chimera Overset](chimera.md), [Multi-block and MPI](multiblock-mpi.md)).

The cases are also part of the regression suite (`Doisneau/MK`, `Doisneau/IG`,
`Doisneau/AG`), which checks the density and number-density fields against a stored reference.
