"""A. Particle relaxation under Stokes drag - exact solution.

A uniform cloud at rest in a uniform gas stream. There are no spatial gradients, so
the momentum equation reduces to du_p/dt = (u_g - u_p)/tau_p and every cell must
follow the same exponential,

    u_p(t) = u_g + (u_p0 - u_g) exp(-t/tau_p),   tau_p = rho_al dp^2 / (18 mu),

which is exactly what ICE's tau reduces to with Cd = 24/Re. The density, the number
density and the temperature must not move at all: the drag work in the energy
equation cancels the kinetic energy it produces, so T_p stays put.

The second half refines the time step at fixed mesh. The field stays uniform, so the
spatial operator contributes nothing and the error is purely the time integration of
the source term - the observed order should be that of the Runge-Kutta scheme.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, observed_order            # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)


def exact(ph, t, up0=0.0):
    return ph.ug + (up0 - ph.ug) * math.exp(-t / ph.tau_stokes)


def relax_case(ph, t_end, cfl=0.8, rk='RK2', nx=8, Lx=1.0, Ly=None, dt_max=None, name='relax'):
    """Uniform cloud, one-way coupled to a uniform gas, run to t_end."""
    case = Case(WORK / name, nx=nx, Lx=Lx, Ly=Ly)
    case.particles(rho=ph.rho_p, u=0.0, v=ph.vg, w=0.0, T=ph.Tg, n=ph.n)
    case.gas(rho=ph.rho_g, u=ph.ug, v=ph.vg, w=0.0, T=ph.Tg,
             R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=t_end, cfl=cfl, rk=rk, drag='Stokes', heat='Stokes',
             rho_al=ph.rho_al, cs=ph.cs, dt_max=dt_max)
    return case.run()


def main():
    ph = Physics()
    rep = Report('A. Stokes drag relaxation')
    tau = ph.tau_stokes
    print('   tau_p = rho_al dp^2 / (18 mu) = %.6f s' % tau)

    # --- The relaxation curve, sampled where the reference table is tabulated ------
    print('   t/tau      u_p/u_g (ICE)   u_p/u_g (exact)    rel. error')
    worst = 0.0
    for ratio in (0.1, 0.5, 1.0, 2.0, 3.0, 5.0):
        sol = relax_case(ph, ratio * tau, name='t%s' % ratio)
        u = sol['var'][U]
        t = sol['time']
        ref = exact(ph, t)
        err = abs(u[0] - ref) / ph.ug
        worst = max(worst, err)
        print('   %5.1f      %.9f     %.9f      %.2e' % (ratio, u[0] / ph.ug, ref / ph.ug, err))

        spread = max(u) - min(u)
        rho, n, temp = sol['var'][RHO], sol['var'][N], sol['var'][T]
        rep.check(spread <= 1e-12 * ph.ug, 't = %.2f tau: uniform to %.1e' % (ratio, spread))
        rep.check(abs(max(rho) - ph.rho_p) <= 1e-12 * ph.rho_p and
                  abs(min(rho) - ph.rho_p) <= 1e-12 * ph.rho_p,
                  't = %.2f tau: density unchanged' % ratio)
        rep.check(max(abs(x - ph.n) for x in n) <= 1e-9 * ph.n,
                  't = %.2f tau: number density unchanged' % ratio)
        rep.check(max(abs(x - ph.Tg) for x in temp) <= 1e-9 * ph.Tg,
                  't = %.2f tau: temperature unchanged (no drag heating)' % ratio)

    rep.check(worst <= 1e-5, 'velocity matches the exact solution to %.1e of u_g' % worst)

    # --- dt-max really is the step, not the step before the CFL factor --------------
    # dt_max here is well below the CFL-limited step, so every step should be exactly
    # dt_max and the end time an exact multiple of it. Were the ceiling applied before
    # the CFL factor, the step would be cfl * dt_max and the end time a multiple of
    # that instead, which is what this distinguishes.
    cfl, dt_max, target = 0.8, 2.0e-6, 0.5 * tau
    sol = relax_case(ph, target, cfl=cfl, dt_max=dt_max, name='dtmax')
    over = sol['time'] - target
    steps = sol['time'] / dt_max
    print('   dt-max = %.0e s: end time %.9e s = %.4f steps, overshoot %.2e s'
          % (dt_max, sol['time'], steps, over))
    rep.check(0.0 <= over <= 1.001 * dt_max, 'dt-max bounds the time step')
    rep.check(abs(steps - round(steps)) <= 1e-6,
              'every step is exactly dt-max, so the ceiling is on the step itself')
    rep.check(abs(sol['var'][U][0] - exact(ph, sol['time'])) / ph.ug <= 1e-5,
              'the shorter step still matches the exact solution')

    # --- Time-step refinement ------------------------------------------------------
    # A transverse velocity common to gas and particles keeps the CFL condition, not
    # ICE's 1e-4 s ceiling, in charge of dt, so halving the CFL halves the step. The
    # cell sizes put dt just at that ceiling at CFL 0.8, making the steps as long as
    # the solver allows so that RK3 stays clear of the round-off floor.
    ph2 = Physics(vg=50.0)
    t_end, floor = 0.02, 1e-12
    cfls = (0.8, 0.4, 0.2)
    for rk, expected in (('RK2', 2.0), ('RK3', 3.0)):
        errors = []
        for cfl in cfls:
            sol = relax_case(ph2, t_end, cfl=cfl, rk=rk, nx=10, Lx=0.01, Ly=5e-3,
                             name='%s-cfl%s' % (rk.lower(), cfl))
            errors.append(abs(sol['var'][U][0] - exact(ph2, sol['time'])) / ph2.ug)
        orders = observed_order(errors, [1.0 / c for c in cfls])
        usable = [o for o, e in zip(orders, errors[1:]) if e > floor]
        print('   %-5s errors %s' % (rk, ' '.join('%.2e' % e for e in errors)))
        print('         orders %s' % ' '.join('%.2f' % o for o in orders))
        if not rep.check(bool(usable), '%s stays above the round-off floor' % rk):
            continue
        rep.check(min(usable) >= expected - 0.35,
                  '%s converges at order %.2f (expected %.0f)' % (rk, min(usable), expected))

    rep.close(WORK)


if __name__ == '__main__':
    main()
