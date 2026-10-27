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
| Time | Local time stepping, 5000 iterations, CFL 0.8, `dt-max` $10^{-4}$ s (the default) |
| Space reconstruction | MUSCL, Van Leer limiter, no shock detector |

The CFL condition sets the local step only during the first part of a run (never for MK, whose
cells step at the `dt-max` ceiling from the first iteration); after it, every cell steps at the
ceiling.

In the exact solution each jet is a 45° band of density 0.1 kg/m³. The bands overlap in a
diamond centred at $(0.25, 0)$, where the density is 0.2 kg/m³.

MK's upper inlet sits one cell higher: with both inlets symmetric about $y = 0$, its two
pressureless jets collide on the face between rows 50 and 51, and which row takes the merged
delta is decided at round-off (it differed between compilers). One cell higher, the jets are
mirror images about the centre of row 51 ($y = 0.005$), which takes the delta.

## Results

IG and AG are set beside free streaming, the kinetic solution they approximate. MK, a
pressureless gas, is set beside the exact solution of its own equations, the **sticky
merger**: each jet streams freely until it reaches the mirror line $y = 0.005$, where the two
join a single jet moving along the line at their common $x$ velocity, 7.07 m/s. The merged
jet moves in $x$ as fast as the jets did, so the column mass per unit $x$ is conserved: once
both jets have reached the line, from $x = 0.305$ on, the merged jet carries 0.02 kg/m²,
2.00 kg/m³ in the cells of row 51.

<figure>
  {% include "vv/images/doisneau-fields.svg" %}
</figure>

<figure>
  {% include "vv/images/doisneau-profiles.svg" %}
</figure>

| Model | Peak density [kg/m³] | Behaviour at the crossing |
|---|---|---|
| Free streaming (exact) | 0.200 | Jets pass through each other |
| MK | 2.00 | Jets merge into a single jet along $y = 0.005$; mass concentrates in the cells of row 51 (delta shock), as in the exact merger |
| IG | 0.346 | Jets merge where they meet, then spread as a single jet |
| AG | 0.234 | Jets cross and continue as two separate, diffused jets; the peak lies beside the crossing |

- **MK** merges the jets into one along $y = 0.005$, midway between its inlets, and piles
  their mass into the cells of row 51: 2.00 kg/m³ downstream of the meeting point, the value
  of the exact merger. This delta shock is the expected monokinetic answer to crossing
  streams.
- **IG** cannot hold two velocities at the same point either. The jets merge where they
  meet, with a peak above the exact diamond's, and its dispersion then spreads the merged
  jet instead of concentrating it.
- **AG** is the only closure that lets the jets cross. Downstream they are much wider and
  weaker than the exact bands.

The three closures run without the [shock detector](../theory/numerics.md#shock-detector).
For MK and AG this is a trade. Without it MK stays nearer its exact solution (an RMS distance
of $5.5 \times 10^{-3}$ kg/m³ against $1.4 \times 10^{-2}$ with the detector; at $x = 0.145$,
upstream of where the jets meet, row 51 holds $8.4 \times 10^{-4}$ kg/m³ against 0.069), and
AG nearer free streaming ($2.6 \times 10^{-2}$ against $2.8 \times 10^{-2}$; the crossing,
0.2 kg/m³ in the exact solution, holds 0.06 on average against 0.04). But neither converges: MK's residual
stalls at $1.6 \times 10^{-5}$ and its jets are $1.2 \times 10^{-2}$ from their mirror image at
the outer edges of the inlets (rows 21 and 81), and AG's residual cycles between
$4 \times 10^{-5}$ and $9 \times 10^{-5}$. AG also gathers mass in two bands just above and
below the crossing, at $0.035 \le |y| \le 0.075$ from $x = 0.175$ to $0.275$: up to
0.234 kg/m³ inside a single jet, where the exact solution is 0.1, and up to 0.16 just beyond
the top and bottom corners of the crossing, outside both jets, where it is empty. With the
detector the field's maximum is 0.115, on each jet just upstream of the crossing (exact 0.1).
Without it AG also clips its pressure tensor at the floor in 22 nearly empty cells ($1$ to
$3 \times 10^{-3}$ kg/m³, where the exact solution is empty): 16 between the jets next to the
inlet face, 6 outside the jets above and below the crossing, at $|y| \approx 0.12$. Both tests
therefore run their 5000 iterations; with the detector, MK reaches the default residual
threshold at iteration 2493 and AG at iteration 6907. IG's residual falls without it, to
$1.7 \times 10^{-6}$ at 5000 iterations apart from brief spikes to $3 \times 10^{-5}$ every 80
iterations or so, and its chimera and multi-block versions match it to round-off
([Chimera Overset](chimera.md), [Multi-block and MPI](multiblock-mpi.md)).

The cases are also part of the regression suite (`Doisneau/MK`, `Doisneau/IG`,
`Doisneau/AG`), which checks the density and number-density fields against a stored reference.
