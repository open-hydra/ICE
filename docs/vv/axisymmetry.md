# Axisymmetric Wedge

An axisymmetric case is run on a wedge one cell thick in the azimuthal direction. Its two
side faces are not boundaries of the flow: they are the planes the solution is rotated
across, so their ghost cells hold the mirror of the neighbouring interior cell (ATLAS code
200), and the axis face holds a plain symmetry condition (code 300). The three cases below
share one mesh, a wedge of one degree about the $x$ axis, $0.1 \times 0.05$ m, with
$40 \times 20 \times 1$ cells. The axis row of nodes sits at $r = 10^{-8}$ m, where ATLAS
places it, so the axis face is a regular face of tiny area rather than a degenerate one.

| Case | Closure | State | What it pins down | Tolerance |
|---|---|---|---|---|
| `Axis/MK` | MK | Cloud moving along $x$, density a function of $r$ only | The side-face and axis ghosts are written; the boundary census | $10^{-9}$ on the axis-row density |
| `Axis/IG` | IG | Uniform cloud at rest, $P = 1$ Pa | The side faces carry the pressure that holds the cloud in equilibrium | $10^{-10}$ on $\lvert\mathbf{u}\rvert/c$, $\rho$ and $P$ |
| `Axis/AG` | AG | Uniform cloud at rest, $\mathbf{P} = \mathbf{I}$ Pa | The same, through the tensor closure | $10^{-10}$ on $\lvert\mathbf{u}\rvert/c$, $\rho$ and $P_{11}, P_{22}, P_{33}$ |

### A pressureless cloud along the axis (`Axis/MK`)

The cloud moves along $x$ at 20 m/s with a density that grows away from the axis,
$\rho_p = 1 + 4\,(r/R)^2$, fed through a mass-flux inlet (code 402) that carries the same
profile. Nothing depends on $x$ and nothing moves radially, so the field is stationary at
any order of the scheme: every $i$-face carries the same flux in and out, and every
$j$-face and axis face carries none. The run makes 200 iterations with RK3, MUSCL and the
Jameson detector, and must keep the axis row ($j = 1$) at its initial density to
$10^{-9}$; it holds it to $2\times10^{-15}$, the roundoff of the sweep, and the whole
field to $8\times10^{-15}$. The inlet payload is written at full precision because the
inlet ghost density is the mass flux divided by the speed: seven significant digits would
shift the rows by up to $2\times10^{-7}$ and fail the check.

The case also checks the boundary census the solver prints at start-up: 80 symmetry faces
(the axis and the outer radius) and 1600 axisymmetry faces (the two side faces), so a side
face read as anything else is caught before the field is looked at. The output must be
finite everywhere and come from the full run (`Time of operation` in the log, and a
solution time of 200 iterations).

`res-threshold = 0` is required: the density residual of this field is at roundoff from
the first iteration, and the default threshold of $10^{-10}$ would end the run there.

### A cloud at rest in a closed wedge (`Axis/IG`, `Axis/AG`)

A uniform cloud at rest, $\rho_p = 1$ kg/m³ and $P = 1$ Pa ($\mathbf{P} = \mathbf{I}$ for
AG), fills the wedge; every face is closed (the four meridional faces are symmetry, code
300, the side faces 200). A gas at rest stays at rest, and discretely this holds to
roundoff: every cell is a closed polyhedron, so $\sum_f A_f \mathbf{n}_f = 0$, and a uniform
state reconstructs to itself, so every face carries exactly $P \mathbf{n} A$. The radial
faces of a cell differ in area by $2\,\Delta x\,\Delta r \sin(\delta/2)$, and it is the
pressure on the two side faces — the hoop term $P/r$ of the axisymmetric equations — that
closes the balance.

The run makes 500 iterations with RK2 and MUSCL, without the detector, and requires, over
every cell, $\max(\lvert u\rvert, \lvert v\rvert, \lvert w\rvert)/c < 10^{-10}$ with
$c = \sqrt{3P/\rho}$ (for AG $P$ is the trace over three), and $\rho$ and the diagonal
pressures within $10^{-10}$ of 1. ICE holds the velocity to $3\times10^{-16}$ of the sound
speed and the density and pressures exactly. The tolerance leaves more than five decades of
margin over that and sits ten decades below the signal it guards against: without the side-face
pressure each cell feels a net force $P\,\delta\,\Delta x\,\Delta r$ towards the axis,
$2P/(\rho\,\Delta r) = 800$ m/s² on the axis row, and the cloud collapses onto the axis —
for IG, $\lvert v\rvert/c = 0.23$ after 50 iterations and 1.1 after 300.

`res-threshold = 0` is required here too: at equilibrium the density residual is exactly
zero, and the default threshold would stop the run after one iteration.
