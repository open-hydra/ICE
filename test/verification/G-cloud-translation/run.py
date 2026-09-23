"""G. A particle cloud released into a uniform frozen gas - exact translation.

Case A relaxes a spatially uniform cloud, so its transport operator never fires; case
C transports a cloud at a fixed velocity, so its source terms never fire. This case is
the product of the two, and it is the smallest problem in which getting either one
wrong, or coupling them wrongly, shows up in the answer.

A frozen gas moves at a uniform `UG` and the cloud starts at rest in a localized packet.
The Stokes relaxation time

    tau = rho_al dp^2 / (18 mu)

depends only on the particle diameter, which stays uniform because `n` is seeded
proportional to rho_p and the two obey the same transport equation. Every particle in
the domain therefore feels the same drag, the particle velocity stays uniform in space
for all time, and it follows the single ODE

    du_p/dt = (UG - u_p) / tau    ->    u_p(t) = UG (1 - e^{-t/tau}),

so the cloud translates rigidly:

    X(t) = UG t - UG tau (1 - e^{-t/tau}),   rho_p(x, t) = rho_p(x - X(t), 0).

Three things are measured, on a periodic mesh:

  * refining the mesh at fixed dt gives the order of the space discretisation under a
    time-varying advection velocity (smooth profile: a wrapped Gaussian);
  * refining dt at a fixed mesh gives the order at which u_p becomes uniform. It is not
    uniform at finite dt: the drag source is evaluated on the density the flux
    divergence is about to change, which leaves an O(dt^2) spread wherever the density
    gradient is steep. That spread vanishing at second order is the statement that
    MK's momentum flux and its drag source are consistent;
  * a top hat carries the invariants that survive any amount of smearing - mass, and
    the position of the centroid.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, observed_order                # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)

L, UG, RHO0 = 1.0, 10.0, 10.0
XC, SIG, HALF = 0.25, 0.06, 0.10        # cloud centre, Gaussian width, top-hat half width
FLOOR = 1.0e-6                          # background haze, keeps the vacuum clamp out of it
DT_REF = 1.0e-5                         # dt ceiling for the mesh-refinement runs


def wrapped_gaussian(x):
    """Gaussian on a periodic line: the image sum makes it exactly L-periodic."""
    return sum(math.exp(-(x - XC - m * L) ** 2 / (2.0 * SIG ** 2)) for m in (-1, 0, 1))


def top_hat(x):
    d = (x - XC + 0.5 * L) % L - 0.5 * L
    return 1.0 if abs(d) <= HALF else 0.0


def profile(shape, x):
    return RHO0 * (shape(x) + FLOOR)


def cell_average(shape, xl, xr, shift, m=32):
    """Midpoint average over one cell of the exact profile displaced by `shift`."""
    h = (xr - xl) / m
    return sum(profile(shape, (xl + (i + 0.5) * h - shift) % L) for i in range(m)) / m


def displacement(t, tau):
    return UG * t - UG * tau * (1.0 - math.exp(-t / tau))


def cloud_case(ph, nx, shape, t_end, dt_max, name):
    case = Case(WORK / name, nx=nx, Lx=L)
    volume = ph.rho_al * (4.0 / 3.0) * math.pi * ph.rp ** 3

    def rho(x):
        return profile(shape, x)

    # u_p = 0 at t = 0; n proportional to rho_p keeps the diameter, hence tau, uniform
    case.particles(rho=rho, u=0.0, v=0.0, w=0.0, T=ph.Tg, n=lambda x: rho(x) / volume)
    case.gas(rho=ph.rho_g, u=UG, v=0.0, w=0.0, T=ph.Tg,
             R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('periodic-x')
    case.ini(t_end=t_end, cfl=0.4, dt_max=dt_max, drag='Stokes', heat='Stokes',
             rho_al=ph.rho_al, cs=ph.cs)
    return case.run()


def centroid(rho, dx):
    """Centre of mass on the periodic line, through the circular mean."""
    ang = sum(r * complex(math.cos(2 * math.pi * (i + 0.5) * dx / L),
                          math.sin(2 * math.pi * (i + 0.5) * dx / L))
              for i, r in enumerate(rho))
    return (math.atan2(ang.imag, ang.real) * L / (2 * math.pi)) % L


def check_common(rep, sol, shape, nx, tau, tag):
    """The checks every run gets: exact u_p, exact mass, exact centroid."""
    t, dx = sol['time'], L / nx
    rho, up = sol['var'][RHO], sol['var'][U]

    exact_u = UG * (1.0 - math.exp(-t / tau))
    mean_u = sum(a * b for a, b in zip(rho, up)) / sum(rho)
    rep.check(abs(mean_u - exact_u) <= 1e-6 * UG,
              '%s: u_p = %.8f vs exact %.8f' % (tag, mean_u, exact_u))

    mass = sum(rho) * dx
    mass0 = sum(cell_average(shape, i * dx, (i + 1) * dx, 0.0) for i in range(nx)) * dx
    rep.check(abs(mass / mass0 - 1.0) <= 1e-12,
              '%s: mass conserved to %.1e' % (tag, abs(mass / mass0 - 1.0)))

    X = displacement(t, tau)
    drift = abs((centroid(rho, dx) - (XC + X)) % L)
    drift = min(drift, L - drift)
    rep.check(drift <= 0.05 * dx,
              '%s: centroid off by %.3f cells' % (tag, drift / dx))
    return X


def main():
    ph = Physics(ug=UG)
    tau = ph.tau_stokes
    t_end = 3.0 * tau
    rep = Report('G. Cloud translation, uniform gas')
    print('   tau = %.4e s, t_end = %.4e s (3 tau), u_p(t_end)/UG = %.6f'
          % (tau, t_end, 1.0 - math.exp(-3.0)))

    # --- 1. space refinement on the smooth profile ---------------------------------
    grids, errors = (50, 100, 200), []
    for nx in grids:
        sol = cloud_case(ph, nx, wrapped_gaussian, t_end, DT_REF, 'gauss-n%d' % nx)
        X = check_common(rep, sol, wrapped_gaussian, nx, tau, 'Gaussian n=%3d' % nx)
        dx = L / nx
        ref = [cell_average(wrapped_gaussian, i * dx, (i + 1) * dx, X) for i in range(nx)]
        errors.append(sum(abs(a - b) for a, b in zip(sol['var'][RHO], ref)) / sum(ref))

    orders = observed_order(errors, grids)
    print('   Gaussian   L1     %s' % ' '.join('%.3e' % e for e in errors))
    print('              orders %s' % ' '.join('%.2f' % o for o in orders))
    rep.check(min(orders) >= 1.5,
              'profile converges in space at order %.2f (expected >= 1.50)' % min(orders))

    # --- 2. time refinement: u_p becomes uniform at second order --------------------
    steps, spreads = (2.0e-5, 1.0e-5, 5.0e-6), []
    for dt_max in steps:
        sol = cloud_case(ph, 100, wrapped_gaussian, t_end, dt_max, 'dt%g' % dt_max)
        up = sol['var'][U]
        spreads.append((max(up) - min(up)) / UG)
    orders = observed_order(spreads, [1.0 / d for d in steps])
    print('   u_p spread        %s' % ' '.join('%.3e' % s for s in spreads))
    print('              orders %s' % ' '.join('%.2f' % o for o in orders))
    rep.check(min(orders) >= 1.8,
              'u_p becomes uniform at order %.2f in dt (expected >= 1.80)' % min(orders))
    rep.check(spreads[-1] <= 1e-6,
              'u_p uniform to %.1e of UG at the finest dt' % spreads[-1])

    # --- 3. the discontinuous profile: invariants only ------------------------------
    nx = 100
    sol = cloud_case(ph, nx, top_hat, t_end, DT_REF, 'tophat-n%d' % nx)
    X = check_common(rep, sol, top_hat, nx, tau, 'top hat  n=%3d' % nx)
    print('   displacement X(t_end) = %.6f m (gas would give %.6f)' % (X, UG * t_end))

    rep.close(WORK)


if __name__ == '__main__':
    main()
