"""I. A particle cloud in a prescribed vortex - exact solution at every Stokes number.

The frozen gas is in solid-body rotation,

    u_g = -Omega y,    v_g = Omega x,

and a compact cloud is released into it seeded at the local gas velocity. This is the
first two-dimensional verification case: cases G and H leave the j-direction fluxes
idle, and a momentum flux that is wrong across directions, or a metric that is wrong on
the j faces, cannot show up in either.

In the complex plane z = x + i y the carrier field is u_g = i Omega z, so Stokes drag
with a diameter that stays uniform gives one linear equation along a trajectory,

    tau z'' + z' = i Omega z,   lambda± = (-1 ± sqrt(1 + 4 i Omega tau)) / (2 tau).

Seeding z'(0) = i Omega z(0) makes it separable, exactly as in case H but over the
complex field: every particle is displaced by the same complex factor,

    z(t) = C(t) z_0,   C = A e^{lambda+ t} + B e^{lambda- t},

with A + B = 1 and A lambda+ + B lambda- = i Omega. The map is a rotation composed with
a dilation - conformal, and non-singular for all time - so the whole Eulerian solution
is closed form,

    rho_p(z, t) = rho_p0(z / C) / |C|^2,    u_p(z, t) = (C'/C) z,

with no trajectory crossing at any Stokes number and therefore no delta-shock: MK is
used inside its domain of validity throughout.

The two numbers that carry the physics are arg C, which lags Omega t, and |C|, which
exceeds 1 because inertia throws the particles outwards. Both are measured against the
closed form over St = Omega tau = 0.01 to 10, which is the transition TEST.md asks for:
at St << 1 the cloud is a tracer, at St ~ 1 it lags appreciably, at St >> 1 it barely
responds to the vortex at all.
"""
import cmath
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, observed_order                # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)

OMEGA, HALF = 1.0, 2.0                  # rotation rate, half width of the square domain
X0, SIG, RHO0 = 0.6, 0.12, 10.0         # cloud centre on the x axis, width, peak
FLOOR = 1.0e-8
QUARTER = 0.5 * math.pi / OMEGA         # a quarter turn of the gas


def stretch(tau, t):
    """C(t) and C'(t) for a cloud seeded at the local gas velocity."""
    root = cmath.sqrt(1.0 + 4.0j * OMEGA * tau)
    lp, lm = (-1.0 + root) / (2.0 * tau), (-1.0 - root) / (2.0 * tau)
    a = (1.0j * OMEGA - lm) / (lp - lm)
    b = (lp - 1.0j * OMEGA) / (lp - lm)
    return (a * cmath.exp(lp * t) + b * cmath.exp(lm * t),
            a * lp * cmath.exp(lp * t) + b * lm * cmath.exp(lm * t))


def profile(x, y):
    return RHO0 * (math.exp(-((x - X0) ** 2 + y ** 2) / (2.0 * SIG ** 2)) + FLOOR)


def exact(x, y, C):
    """rho_p0 pulled back through the conformal map, with the area Jacobian |C|^2."""
    z0 = complex(x, y) / C
    return profile(z0.real, z0.imag) / abs(C) ** 2


def vortex_case(ph, n, tau, t_end, name):
    dp = math.sqrt(18.0 * ph.mu * tau / ph.rho_al)
    volume = ph.rho_al * (4.0 / 3.0) * math.pi * (0.5 * dp) ** 3
    case = Case(WORK / name, nx=n, ny=n, Lx=2.0 * HALF, Ly=2.0 * HALF,
                Lz=2.0 * HALF / n, x0=-HALF, y0=-HALF)
    # Seeded at the local gas velocity; n proportional to rho_p keeps dp, hence tau,
    # uniform, which is what makes the trajectory equation linear.
    case.particles(rho=profile,
                   u=lambda x, y: -OMEGA * y, v=lambda x, y: OMEGA * x, w=0.0,
                   T=ph.Tg, n=lambda x, y: profile(x, y) / volume)
    case.gas(rho=ph.rho_g,
             u=lambda x, y: -OMEGA * y, v=lambda x, y: OMEGA * x, w=0.0,
             T=ph.Tg, R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=t_end, cfl=0.4, dt_max=1.0e-3, drag='Stokes', heat='Stokes',
             rho_al=ph.rho_al, cs=ph.cs)
    return case, case.run()


def measure(case, sol):
    """Cloud mass, complex centroid and the complex slope of the velocity field."""
    rho, up, vp = sol['var'][RHO], sol['var'][U], sol['var'][V]
    pts = [complex(x, y) for x, y in case.centres]
    mass = sum(rho)
    centre = sum(r * z for r, z in zip(rho, pts)) / mass
    slope = (sum(r * complex(a, b) * z.conjugate() for r, a, b, z in zip(rho, up, vp, pts))
             / sum(r * abs(z) ** 2 for r, z in zip(rho, pts)))
    return mass, centre, slope


def l1_error(case, sol, C):
    ref = [exact(x, y, C) for x, y in case.centres]
    return sum(abs(a - b) for a, b in zip(sol['var'][RHO], ref)) / sum(ref)


def main():
    ph = Physics()
    rep = Report('I. Cloud in a prescribed vortex')

    # --- 1. mesh refinement against the closed form, at St = 0.1 -------------------
    tau = 0.1 / OMEGA
    grids, errors, drifts = (48, 96, 192), [], []
    for n in grids:
        case, sol = vortex_case(ph, n, tau, QUARTER, 'refine-n%d' % n)
        C, _ = stretch(tau, sol['time'])
        errors.append(l1_error(case, sol, C))
        mass, centre, _ = measure(case, sol)
        h = 2.0 * HALF / n
        drifts.append(abs(centre - C * X0))
        rep.check(drifts[-1] <= 0.1 * h,
                  'n=%3d: centroid (%+.5f,%+.5f) vs exact (%+.5f,%+.5f), %.3f cells'
                  % (n, centre.real, centre.imag, (C * X0).real, (C * X0).imag,
                     drifts[-1] / h))
    orders = observed_order(errors, grids)
    drift_orders = observed_order(drifts, grids)
    print('   St = 0.1, quarter turn:  L1       %s' % ' '.join('%.3e' % e for e in errors))
    print('                            orders   %s' % ' '.join('%.2f' % o for o in orders))
    print('                            centroid %s' % ' '.join('%.3e' % d for d in drifts))
    print('                            orders   %s' % ' '.join('%.2f' % o for o in drift_orders))
    rep.check(min(orders) >= 1.5,
              'the profile converges at order %.2f (expected >= 1.50)' % min(orders))
    rep.check(min(drift_orders) >= 1.5,
              'the centroid converges at order %.2f (expected >= 1.50)' % min(drift_orders))

    # --- 2. the Stokes transition ---------------------------------------------------
    n = 96
    print('   Stokes sweep at n=%d over a quarter turn (gas turns %.4f rad):' % (n, QUARTER))
    print('     St      |C| ICE   |C| exact  lag [rad]  exact lag   L1')
    lags, radii = [], []
    for St in (0.01, 0.1, 1.0, 10.0):
        tau_s = St / OMEGA
        case, sol = vortex_case(ph, n, tau_s, QUARTER, 'st%g' % St)
        C, Cdot = stretch(tau_s, sol['time'])
        mass, centre, slope = measure(case, sol)
        err = l1_error(case, sol, C)

        scale, exact_scale = abs(centre) / X0, abs(C)
        lag = OMEGA * sol['time'] - cmath.phase(centre)
        exact_lag = OMEGA * sol['time'] - cmath.phase(C)
        lags.append(lag)
        radii.append(scale)
        print('     %-6g  %.6f  %.6f   %.6f   %.6f   %.3e'
              % (St, scale, exact_scale, lag, exact_lag, err))

        rep.check(abs(scale / exact_scale - 1.0) <= 5e-3,
                  'St=%-5g: radial dilation %.6f vs exact %.6f' % (St, scale, exact_scale))
        rep.check(abs(lag - exact_lag) <= 5e-3,
                  'St=%-5g: phase lag %.6f vs exact %.6f rad' % (St, lag, exact_lag))
        # The velocity field itself, not just where the cloud ended up
        rep.check(abs(slope / (Cdot / C) - 1.0) <= 1e-2,
                  'St=%-5g: du_p/dz = %.5f%+.5fj vs exact %.5f%+.5fj'
                  % (St, slope.real, slope.imag, (Cdot / C).real, (Cdot / C).imag))
        rep.check(err <= 0.15,
                  'St=%-5g: L1 error %.3e against the exact map' % (St, err))
        rep.check(abs(mass / sum(exact(x, y, C) for x, y in case.centres) - 1.0) <= 1e-3,
                  'St=%-5g: mass within %.1e of the exact solution'
                  % (St, abs(mass / sum(exact(x, y, C) for x, y in case.centres) - 1.0)))

    # --- 3. the physics the sweep is there to show ----------------------------------
    rep.check(all(a < b for a, b in zip(lags[:-1], lags[1:])),
              'the phase lag grows monotonically with St')
    rep.check(all(a < b for a, b in zip(radii[:-1], radii[1:])),
              'the centrifugal ejection grows monotonically with St')
    rep.check(lags[0] <= 5e-3 and abs(radii[0] - 1.0) <= 0.02,
              'at St=0.01 the cloud is a tracer: lag %.1e rad, dilation %.4f'
              % (lags[0], radii[0]))
    rep.check(lags[-1] >= 0.4,
              'at St=10 the cloud barely responds: lag %.4f rad' % lags[-1])

    rep.close(WORK)


if __name__ == '__main__':
    main()
