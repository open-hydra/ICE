"""H. A particle cloud in a linearly straining frozen gas - exact affine solution.

Case G moves a cloud through a gas of one velocity, so the particle velocity stays
uniform and the cloud only translates. Here the frozen gas carries a uniform strain,

    u_g(x) = U0 + a x,

and the cloud can no longer move rigidly: the particles lag the gas acceleration, which
is the defining behaviour of a dispersed phase. It is the first case in which the
particle velocity field itself is non-trivial, so it is the first that can detect a
momentum flux that is wrong in a way a uniform velocity hides.

With xi = x + U0/a the gas velocity is a xi, and Stokes drag with a diameter that stays
uniform gives one linear ODE along a trajectory,

    tau xi'' + xi' - a xi = 0,    lambda± = (-1 ± sqrt(1 + 4 a tau)) / (2 tau).

Seeding the cloud at the local gas velocity, xi'(0) = a xi(0), makes the solution
separable - every particle is displaced by the same factor,

    xi(t) = C(t) xi_0,    C = A e^{lambda+ t} + B e^{lambda- t},

with A + B = 1 and A lambda+ + B lambda- = a. The map x_0 -> x(t) is therefore affine,
which gives the whole Eulerian solution in closed form:

    rho_p(x, t) = rho_p0(xi / C) / C,     u_p(x, t) = (C'/C) xi,

exact at every Stokes number, with no crossing at any time since C > 0.

Two things are measured. The profile and the velocity field are compared against that
closed form under mesh refinement. Then, separately, the slip

    u_g - u_p = (a - C'/C) xi

is compared against the small-Stokes asymptote u_p ~ u_g - tau Du_g/Dt, which here is
a^2 tau xi: the departure has to fall linearly in St = a tau, and that is what pins the
solver to the right asymptotic behaviour rather than merely to a plausible one.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, observed_order                # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)

L, U0, A_STRAIN = 1.0, 2.0, 6.0         # domain, u_g = U0 + a x
X0, SIG, RHO0 = 0.20, 0.035, 10.0       # initial cloud
FLOOR = 1.0e-6
T_END, DT_MAX = 0.10, 2.0e-5


def xi(x):
    return x + U0 / A_STRAIN


def stretch(tau, t):
    """C(t) and C'(t) for a cloud seeded at the local gas velocity."""
    root = math.sqrt(1.0 + 4.0 * A_STRAIN * tau)
    lp, lm = (-1.0 + root) / (2.0 * tau), (-1.0 - root) / (2.0 * tau)
    a = (A_STRAIN - lm) / (lp - lm)
    b = (lp - A_STRAIN) / (lp - lm)
    return (a * math.exp(lp * t) + b * math.exp(lm * t),
            a * lp * math.exp(lp * t) + b * lm * math.exp(lm * t))


def profile(x):
    return RHO0 * (math.exp(-(x - X0) ** 2 / (2.0 * SIG ** 2)) + FLOOR)


def strain_case(ph, nx, dp, name):
    case = Case(WORK / name, nx=nx, Lx=L)
    rp = 0.5 * dp
    volume = ph.rho_al * (4.0 / 3.0) * math.pi * rp ** 3
    # Seed u_p at the local gas velocity; n proportional to rho_p keeps dp, hence tau,
    # uniform, which is what makes the trajectory ODE linear and the map affine.
    case.particles(rho=profile, u=lambda x: U0 + A_STRAIN * x, v=0.0, w=0.0, T=ph.Tg,
                   n=lambda x: profile(x) / volume)
    case.gas(rho=ph.rho_g, u=lambda x: U0 + A_STRAIN * x, v=0.0, w=0.0, T=ph.Tg,
             R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=T_END, cfl=0.4, dt_max=DT_MAX, drag='Stokes', heat='Stokes',
             rho_al=ph.rho_al, cs=ph.cs)
    return case.run()


def cell_average(f, xl, xr, m=32):
    h = (xr - xl) / m
    return sum(f(xl + (i + 0.5) * h) for i in range(m)) / m


def main():
    ph = Physics()
    rep = Report('H. Cloud in a linear strain')

    # --- 1. the closed-form solution under mesh refinement --------------------------
    dp = ph.dp
    tau = ph.rho_al * dp ** 2 / (18.0 * ph.mu)
    print('   a = %.1f 1/s, tau = %.4e s, St = a tau = %.4f' % (A_STRAIN, tau, A_STRAIN * tau))

    grids, errors = (100, 200, 400), []
    for nx in grids:
        sol = strain_case(ph, nx, dp, 'exact-n%d' % nx)
        t, dx = sol['time'], L / nx
        C, Cdot = stretch(tau, t)
        rho, up = sol['var'][RHO], sol['var'][U]
        xc = [(i + 0.5) * dx for i in range(nx)]

        ref = [cell_average(lambda x: profile(xi(x) / C - U0 / A_STRAIN) / C,
                            i * dx, (i + 1) * dx) for i in range(nx)]
        errors.append(sum(abs(a - b) for a, b in zip(rho, ref)) / sum(ref))

        # Centroid: the affine map sends the Gaussian centre to C xi_0
        mass = sum(rho)
        xbar = sum(r * x for r, x in zip(rho, xc)) / mass
        exact_x = C * xi(X0) - U0 / A_STRAIN
        rep.check(abs(xbar - exact_x) <= 0.05 * dx,
                  'n=%3d: centroid %.6f vs exact %.6f (%.3f cells)'
                  % (nx, xbar, exact_x, abs(xbar - exact_x) / dx))

        # Velocity field: mass-weighted least squares of u_p against xi, through 0
        slope = (sum(r * u * xi(x) for r, u, x in zip(rho, up, xc))
                 / sum(r * xi(x) ** 2 for r, x in zip(rho, xc)))
        rep.check(abs(slope / (Cdot / C) - 1.0) <= 2e-4,
                  'n=%3d: du_p/dxi = %.6f vs exact %.6f' % (nx, slope, Cdot / C))

        mass_ex = sum(ref)
        rep.check(abs(mass / mass_ex - 1.0) <= 2e-3,
                  'n=%3d: cloud mass within %.1e of the exact stretch'
                  % (nx, abs(mass / mass_ex - 1.0)))

    orders = observed_order(errors, grids)
    print('   C(t_end) = %.6f  (cloud stretched by that factor)' % C)
    print('   L1     %s' % ' '.join('%.3e' % e for e in errors))
    print('   orders %s' % ' '.join('%.2f' % o for o in orders))
    rep.check(min(orders) >= 1.5,
              'converges at order %.2f (expected >= 1.50)' % min(orders))

    # --- 2. the small-Stokes asymptote ----------------------------------------------
    # The exact slip coefficient is a - C'/C, and expanding it in St = a tau gives
    # a^2 tau (1 - 2 St + O(St^2)): the asymptote u_p ~ u_g - tau Du_g/Dt is the
    # leading term, and the departure from it must be proportional to St.
    print('   slip coefficient against a^2 tau (the u_g - tau Du_g/Dt asymptote):')
    ratios, stokes = [], (0.20, 0.10, 0.05, 0.025)
    for St in stokes:
        tau_s = St / A_STRAIN
        dp_s = math.sqrt(18.0 * ph.mu * tau_s / ph.rho_al)
        sol = strain_case(ph, 200, dp_s, 'st%g' % St)
        t, dx = sol['time'], L / 200
        rho, up = sol['var'][RHO], sol['var'][U]
        xc = [(i + 0.5) * dx for i in range(200)]
        slope = (sum(r * u * xi(x) for r, u, x in zip(rho, up, xc))
                 / sum(r * xi(x) ** 2 for r, x in zip(rho, xc)))
        C, Cdot = stretch(tau_s, t)
        measured = A_STRAIN - slope                      # slip coefficient, per unit xi
        ratios.append(measured / (A_STRAIN ** 2 * tau_s))

        # First, the exact solution: this has to hold at every St, asymptotic or not
        rep.check(abs(slope / (Cdot / C) - 1.0) <= 5e-4,
                  'St=%-6.3f du_p/dxi = %.6f vs exact %.6f' % (St, slope, Cdot / C))
        print('     St=%-6.3f slip/xi = %.6f, asymptote %.6f, ratio %.5f'
              % (St, measured, A_STRAIN ** 2 * tau_s, ratios[-1]))

    rep.check(all(a < b for a, b in zip(ratios[:-1], ratios[1:])),
              'the asymptote is approached monotonically as St falls')
    rep.check(all(1.3 * s <= 1.0 - r <= 2.2 * s for r, s in zip(ratios, stokes)),
              'the departure from the asymptote is proportional to St, coefficient ~2')
    rep.check(abs(ratios[-1] - 1.0) <= 0.05,
              'slip is within %.1f%% of the asymptote at St = %.3f'
              % (100 * abs(ratios[-1] - 1.0), stokes[-1]))

    rep.close(WORK)


if __name__ == '__main__':
    main()
