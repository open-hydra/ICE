# Gaussian Closures

Three properties of the IG and AG closures that the crossing jets and the Berthon problems
cannot see, because they compare density, velocity and pressure only: that the temperature
is carried with the mass, that the fastest wave in a direction is $\sqrt{3P_{nn}/\rho}$ and
not an isotropic average, and that a symmetry plane reflects the pressure tensor. Each is a
case directory under `test/` with its own `make_case.py` (the inputs are regenerated from it)
and a `verify.py` that compares against a solution that does not come from ICE.

| Case | What it pins down | Reference | Tolerance | Measured |
|---|---|---|---|---|
| `Thermal/IG` | Pressure work $2.5\,P u_n$ in the IG energy flux | $T$ uniform | $\max\lvert T-1\rvert < 10^{-9}$ | $1.0\times10^{-15}$ |
| `Thermal/AG` | Pressure work $\mathbf{u}\cdot(\mathsf{P}\mathbf{n})$ in the AG energy flux | $T$ uniform | $\max\lvert T-1\rvert < 10^{-9}$ | $2.0\times10^{-15}$ |
| `Riemann/AG` | The wave speed $\sqrt{3P_{nn}/\rho}$ in Rusanov and in the time step | Exact γ = 3 Riemann solution | L1(ρ) ≤ 0.05, TV ≤ 1.05 × exact | L1 0.0200, TV 1.001 / 1.009 |
| `Reflect/AG` | The reflected tensor $\mathsf{H}\mathsf{P}\mathsf{H}$ at a code-200 plane | $M_x$ conserved | $\lvert\Delta M_x\rvert/M_x < 10^{-9}$ | $2.7\times10^{-14}$ |

### Thermal/IG and Thermal/AG: the temperature of an expanding blob

A Gaussian blob, $\rho = 0.01 + e^{-r^2/(2\sigma^2)}$ with $\sigma = 0.05$, at rest in the centre
of a closed $1\times1$ box (100 × 100 cells, walls 300 on the four sides, a slab 0.01 thick),
with $P = \rho\theta$, $\theta = 1$ (for AG the isotropic tensor $P_{11} = P_{22} = P_{33} = \rho$),
and a uniform $T = 1$. It expands for $t = 0.1$, RK3 with MUSCL and the van Leer limiter.

The energy flux carries the pressure work, $2.5\,P u_n$ for IG and $\mathbf{u}\cdot(\mathsf{P}\mathbf{n})$
for AG. With it, the flux of the thermal energy $E - \tfrac12\operatorname{tr}\mathsf{M}$ is
$c_s T$ times the mass flux wherever $T$ is uniform: the central parts of the two fluxes differ
by exactly the pressure work, and the Rusanov dissipation acts on both alike. A uniform $T$
therefore stays uniform to round-off, however violently the blob expands. The case checks
$\max|T-1| < 10^{-9}$ over every cell with $\rho > 10^{-3}$ (all of them); it measures
$1.0\times10^{-15}$ (IG) and $2.0\times10^{-15}$ (AG). Without the pressure work the thermal
energy is not advected with the mass and $T$ moves by $O(\theta/(c_s T)) = O(1)$: 1.0 for IG,
where the expansion drives $T$ to zero and the floor clips it, and 10.9 for AG.

Two settings make that signature $O(1)$ and must not be changed. `specific-heat = 1` puts the
thermal energy $\rho c_s T = \rho$ on the scale of the pressure energy $\tfrac32\rho$; and
`INPUT/` holds no `part-properties.dat`, because a table overrides the INI value (a copied
Doisneau table, $c_p = 1500$, would shrink the signature to $6.7\times10^{-4}$). `verify.py`
refuses to judge when a table is present. The slab thickness is 0.01, not the writer's
default $10^{-3}$: the time step takes the smallest of the three directional limits, and a thin
slab would make $z$ bind. IG runs at CFL 0.5, AG at 0.8. The step is the minimum of the
per-direction limits, so on the diagonals of a radial flow the two Courant numbers add; for the
odd-even mode of a first-order Rusanov flux that puts RK3 at $z = -4\,\mathrm{CFL}$, outside its
real stability interval $[-2.51, 0]$ at CFL 0.8. The IG blob, the faster of the two, grows that
mode and floors 24 cells by $t = 0.08$, which spoils its mass balance by 0.5 %; at CFL 0.5 it
stays clean.

The total energy $\sum E V$ is printed but is not the gate: the scheme is conservative whatever
the flux, so the sum moves only through floors and boundaries. It is conserved to
$5\times10^{-15}$ (IG) and $6\times10^{-13}$ (AG).

### Riemann/AG: an anisotropic Sod tube against the γ = 3 solution

400 × 4 cells over $[0, 1]$ (the slab 0.01 thick, faces 3 and 4 walls, faces 1 and 2
extrapolation), left $\rho = 1$, $P_{11} = 10^5$, right $\rho = 0.125$, $P_{11} = 10^4$, at rest,
$P_{22} = P_{33} = 0.01\rho$, no shear, $T = 300$. Explicit Euler at first order, CFL 0.8, to
$t = 2.7\times10^{-4}$ s. The pressures are $10^5$ times those of a unit problem so that the CFL
step (3.7 µs) stays under `dt-max`, which the INI sets to $10^{-3}$ s.

In one dimension with no shear, $(\rho, \rho u, \rho u^2 + P_{11})$ obey the Euler equations
with γ = 3, and $P_{22}/\rho$, $P_{33}/\rho$ and $T$ are carried with the mass, so the exact
solution is the γ = 3 Riemann solution: a rarefaction to the left, a contact, a shock to the
right, $p^* = 27291$, $u^* = 192.45$. `verify.py` computes it (Toro, *Riemann Solvers*, ch. 4:
Newton on the pressure function, then the standard sampling) after checking its own solver
against Toro's Sod values. The checks, on the row average:

* every value finite and $\min P_{11} > 10^3$: no floor hits (the floor is $10^{-25}$);
* $\sum|\rho - \rho_{ex}| / \sum\rho_{ex} \le 0.05$: measured 0.0200, the smearing of the contact
  and the shock by first-order Rusanov over 400 cells;
* $\max|P_{22}/\rho - 0.01| < 10^{-9}$: measured $6\times10^{-17}$;
* the total variation of $\rho$ and of $P_{11}$ at most 1.05 times the exact one: measured 1.001
  and 1.009.

The last check is what sees the wave speed. The fastest wave in $x$ moves at $\sqrt{3P_{11}/\rho}$,
$\sqrt3$ faster than the isotropic $\sqrt{\operatorname{tr}\mathsf{P}/\rho}$ on these states. With
the isotropic speed the Rusanov flux is under-dissipated and the step too long, but the run does
not blow up: the global step follows the fastest cell, and it shrinks from 6.3 µs to 2.6 µs as
the star region speeds up, which brings the left state back to the edge of stability. The field
oscillates instead: the total variation is 1.56 ($\rho$) and 3.59 ($P_{11}$) times the exact one
and $u$ overshoots $u^*$ by 61 %, while the L1 error of $\rho$, 0.0175, is as small as the
correct scheme's. An L1 norm alone cannot tell the two apart.

### Reflect/AG: no shear through a symmetry plane

The upper half of a sheared blob moving along a symmetry plane: 200 × 100 cells over $2\times1$
(slab 0.01 thick), face 3 ($y = 0$) code 200, faces 1, 2 and 4 extrapolation. The blob is
$\rho = 0.01 + b$, $b = e^{-((x-1)^2+y^2)/(2\sigma^2)}$, $\sigma = 0.1$, moving at $u = 1$, with
$P_{11} = P_{22} = P_{33} = \rho$, $P_{12} = 0.3\,b\,y/\sigma$ and $T = 1$, run to $t = 0.1$ with RK3,
MUSCL and the van Leer limiter. $P_{12}$ is odd in $y$ and zero on the plane: the field is the
upper half of a symmetric full-domain one.

Code 200 mirrors both ghost layers, so the reconstructed states at the plane are an exact
mirror pair. The reflected tensor $\mathsf{H}\mathsf{P}\mathsf{H}$, $\mathsf{H} = \mathsf{I} - 2\mathbf{n}\mathbf{n}^T$,
flips $P_{12}$, the plane's $x$-momentum flux $\tfrac12(P_{12,l} + P_{12,r})$ vanishes, and the
total $x$-momentum $M_x = \sum\rho u V$ can change only through faces 1, 2 and 4, where the
background's fluxes cancel exactly and the blob's tails are $e^{-40}$. The case checks
$|M_x(t) - M_x(0)|/M_x(0) < 10^{-9}$ and measures $2.7\times10^{-14}$. With the tensor copied
into the ghost instead of reflected, the plane carries the shear $P_{12}$ of the first cell and
$M_x$ drifts by $2.0\times10^{-2}$.

The shear is carried by the blob alone. On the background, $0.3\cdot0.01\,y/\sigma$ would leave
through face 4 as the $x$-momentum flux $P_{12}$ and move $M_x$ by about 0.1 whatever the plane
does. The plane is a code-200 face and not a 300 wall, because a 300 face fills its second ghost
by extrapolation: the pair at the plane is then not a mirror pair, and a $10^{-9}$ threshold
would not be defensible on correct code.
