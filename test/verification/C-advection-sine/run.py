"""C. Sinusoidal density wave on a periodic mesh - exact solution.

Cases A and B leave the transport operator idle. Here there is no gas at all (0-way
coupling, no drag), the cloud moves at a uniform velocity U and the density carries a
sine. Free transport gives

    rho_p(x, t) = rho_0 + A sin(k (x - U t)),   k = 2 pi / L,

so after any time the profile is the initial one, shifted. The velocity must stay
uniform: a pressureless cloud moving at one speed has nothing to change it.

Refining the mesh at a fixed time step (ICE's 1e-4 s ceiling holds dt constant across
these grids, so the time error is the same for all of them) gives the order of the
space discretisation. Two reconstructions are compared: MUSCL with the Van Leer
limiter, and IORD, which is first order.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, observed_order            # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)

L, SPEED, RHO0, AMP, T_END = 1.0, 1.0, 10.0, 5.0, 0.25
K = 2.0 * math.pi / L


def cell_average(xl, xr, t):
    """Exact cell average of rho_0 + A sin(k(x - U t)) over [xl, xr]."""
    phase = K * (xr - SPEED * t), K * (xl - SPEED * t)
    return RHO0 + AMP * (math.cos(phase[1]) - math.cos(phase[0])) / (K * (xr - xl))


def advect_case(ph, nx, reconstruction, limiter, name):
    case = Case(WORK / name, nx=nx, Lx=L)
    volume = ph.rho_al * (4.0 / 3.0) * math.pi * ph.rp ** 3

    def rho(x):
        return RHO0 + AMP * math.sin(K * x)

    case.particles(rho=rho, u=SPEED, v=0.0, w=0.0, T=ph.Tg,
                   n=lambda x: rho(x) / volume)
    case.boundaries('periodic-x')
    case.ini(t_end=T_END, reconstruction=reconstruction, limiter=limiter,
             rho_al=ph.rho_al, cs=ph.cs)
    return case.run()


def main():
    ph = Physics()
    rep = Report('C. Periodic sine advection')
    grids = (25, 50, 100, 200)
    results = {}

    for reconstruction, limiter, tag in (('MUSCL', 'VANLEER', 'MUSCL/VanLeer'),
                                         ('MUSCL', 'IORD', 'IORD (1st order)')):
        errors = []
        for nx in grids:
            sol = advect_case(ph, nx, reconstruction, limiter,
                              '%s-n%d' % (limiter.lower(), nx))
            dx = L / nx
            rho = sol['var'][RHO]
            ref = [cell_average(i * dx, (i + 1) * dx, sol['time']) for i in range(nx)]
            errors.append(sum(abs(a - b) for a, b in zip(rho, ref)) / nx / AMP)

            mass = sum(rho) * dx
            mass0 = sum(cell_average(i * dx, (i + 1) * dx, 0.0) for i in range(nx)) * dx
            rep.check(abs(mass - mass0) <= 1e-12 * mass0,
                      '%s n=%3d: mass conserved to %.1e' % (tag, nx, abs(mass / mass0 - 1)))
            rep.check(max(abs(u - SPEED) for u in sol['var'][U]) <= 1e-10 * SPEED,
                      '%s n=%3d: velocity stays uniform' % (tag, nx))

        orders = observed_order(errors, grids)
        results[tag] = errors
        print('   %-18s errors %s' % (tag, ' '.join('%.3e' % e for e in errors)))
        print('   %-18s orders %s' % ('', ' '.join('%.2f' % o for o in orders)))
        expected = 1.5 if limiter == 'VANLEER' else 0.85
        rep.check(min(orders) >= expected,
                  '%s converges at order %.2f (expected >= %.2f)' % (tag, min(orders), expected))

    rep.check(all(a < b for a, b in zip(results['MUSCL/VanLeer'], results['IORD (1st order)'])),
              'MUSCL is more accurate than first order on every grid')
    rep.close(WORK)


if __name__ == '__main__':
    main()
