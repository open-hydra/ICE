# Axisymmetric Wedge

An axisymmetric case runs on a wedge one cell thick in the azimuthal direction.

- The two **side faces** are not boundaries of the flow: they are the planes the solution is rotated across. Their
  ghost cells hold the mirror of the neighbouring interior cell (code `200`).
- The **axis face** holds a plain symmetry condition (code `300`).

The three cases share one mesh: a wedge of one degree about the $x$ axis, $0.1 \times 0.05$ m, with
$40 \times 20 \times 1$ cells. The axis row of nodes sits at $r = 10^{-8}$ m, where ATLAS places it, so the axis face
is a regular face of tiny area rather than a degenerate one.

| Case | Closure | State | What it pins down | Tolerance |
|---|---|---|---|---|
| `Axis/MK` | MK | Cloud moving along $x$, density a function of $r$ only | The side-face and axis ghosts are written; the boundary census | $10^{-9}$ on the axis-row density |
| `Axis/IG` | IG | Uniform cloud at rest, $P = 1$ Pa | The side faces carry the pressure that holds the cloud in equilibrium | $10^{-10}$ on $\lvert\mathbf{u}\rvert/c$, $\rho$ and $P$ |
| `Axis/AG` | AG | Uniform cloud at rest, $\mathbf{P} = \mathbf{I}$ Pa | The same, through the tensor closure | $10^{-10}$ on $\lvert\mathbf{u}\rvert/c$, $\rho$ and $P_{11}, P_{22}, P_{33}$ |

## A pressureless cloud along the axis (`Axis/MK`)

**Setup.** The cloud moves along $x$ at 20 m/s with a density that grows away from the axis,
$\rho_p = 1 + 4\,(r/R)^2$, fed through a mass-flux inlet (code `402`) with the same profile. The run makes 200
iterations with RK3, MUSCL and the shock detector.

**Expected.** Nothing depends on $x$ and nothing moves radially, so the field is stationary at any order of the scheme:
every $i$-face carries the same flux in and out, and every $j$-face and axis face carries none.

**Checks.**

- The axis row ($j = 1$) keeps its initial density to $10^{-9}$. ICE holds it to $2\times10^{-15}$, the round-off of
  the sweep, and the whole field to $8\times10^{-15}$.
- The boundary census printed at start-up shows 80 symmetry faces (the axis and the outer radius) and 1600
  axisymmetry faces (the two side faces), so a side face read as anything else is caught.
- The output is finite everywhere and comes from the full run.

**Notes.**

- The inlet values are written at full precision: the inlet ghost density is the mass flux divided by the speed, and
  seven significant digits would shift the rows by up to $2\times10^{-7}$ and fail the check.
- `res-threshold = 0` makes the run end on `iter-threshold`: the density residual is at round-off from the first
  iteration, so it never falls by the relative threshold.

## A cloud at rest in a closed wedge (`Axis/IG`, `Axis/AG`)

**Setup.** A uniform cloud at rest, $\rho_p = 1$ kg/m³ and $P = 1$ Pa ($\mathbf{P} = \mathbf{I}$ for AG), fills the
wedge. Every face is closed: the four meridional faces are symmetry (`300`), the side faces `200`. The run makes 500
iterations with RK2 and MUSCL, without the shock detector.

**Expected.** A cloud at rest stays at rest, and discretely this holds to round-off:

- every cell is a closed polyhedron, so $\sum_f A_f \mathbf{n}_f = 0$;
- a uniform state reconstructs to itself, so every face carries exactly $P \mathbf{n} A$;
- the radial faces of a cell differ in area by $2\,\Delta x\,\Delta r \sin(\delta/2)$, and the pressure on the two
  side faces, the hoop term $P/r$ of the axisymmetric equations, closes the balance.

**Checks.** Over every cell, $\max(\lvert u\rvert, \lvert v\rvert, \lvert w\rvert)/c < 10^{-10}$ with
$c = \sqrt{3P/\rho}$ (for AG, $P$ is the trace over three), and $\rho$ and the diagonal pressures within $10^{-10}$
of 1. ICE holds the velocity to $3\times10^{-16}$ of the sound speed, and the density and pressures exactly.

**Why the tolerance catches a fault.** Without the side-face pressure each cell feels a net force
$P\,\delta\,\Delta x\,\Delta r$ toward the axis, $2P/(\rho\,\Delta r) = 800$ m/s² on the axis row, and the cloud
collapses onto the axis: for IG, $\lvert v\rvert/c = 0.23$ after 50 iterations and 1.1 after 300. The tolerance sits
ten decades below that signal and five above what ICE achieves.

As in `Axis/MK`, `res-threshold = 0` makes the run end on `iter-threshold`: at equilibrium the density residual is
exactly zero.

## Running them

```bash
ctest --test-dir build -R Axis --output-on-failure
```

The three cases carry the `fast` label, so the pre-push hook runs them.
