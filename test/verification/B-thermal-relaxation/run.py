"""B. Thermal relaxation of a particle cloud - exact solution.

The companion of case A for the energy equation. The particles move with the gas, so
there is no slip, no drag force and no drag work; with heat = Stokes the Nusselt
number is exactly 2 and the convective exchange term reduces to a linear relaxation

    dT_p/dt = (T_g - T_p)/tau_T,   tau_T = rho_al cs dp^2 / (12 kg),

so that T_p(t) = T_g + (T_p0 - T_g) exp(-t/tau_T). Radiation is switched off with
emiss = 0; it is the only other term in the energy source.

The velocity, density and number density must not move.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report                            # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)
TP0 = 400.0                      # particles start hotter than the gas


def exact(ph, t):
    return ph.Tg + (TP0 - ph.Tg) * math.exp(-t / ph.tau_thermal)


def thermal_case(ph, t_end, name):
    """Cloud in equilibrium with the gas in velocity, out of equilibrium in temperature."""
    case = Case(WORK / name, nx=8, Lx=1.0)
    case.particles(rho=ph.rho_p, u=ph.ug, v=0.0, w=0.0, T=TP0, n=ph.n)
    case.gas(rho=ph.rho_g, u=ph.ug, v=0.0, w=0.0, T=ph.Tg,
             R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=t_end, drag='Stokes', heat='Stokes', rho_al=ph.rho_al, cs=ph.cs)
    return case.run()


def main():
    ph = Physics()
    rep = Report('B. Thermal relaxation (Nu = 2)')
    tau = ph.tau_thermal
    print('   tau_T = rho_al cs dp^2 / (12 kg) = %.6f s' % tau)

    print('   t/tau_T    T_p [K] (ICE)   T_p [K] (exact)    rel. error')
    worst = 0.0
    for ratio in (0.1, 0.5, 1.0, 2.0, 3.0, 5.0):
        sol = thermal_case(ph, ratio * tau, 't%s' % ratio)
        temp = sol['var'][T]
        ref = exact(ph, sol['time'])
        err = abs(temp[0] - ref) / (TP0 - ph.Tg)
        worst = max(worst, err)
        print('   %5.1f      %.9f   %.9f      %.2e' % (ratio, temp[0], ref, err))

        rep.check(max(temp) - min(temp) <= 1e-10 * TP0,
                  't = %.2f tau_T: uniform to %.1e K' % (ratio, max(temp) - min(temp)))
        rep.check(max(abs(x - ph.ug) for x in sol['var'][U]) <= 1e-12 * ph.ug,
                  't = %.2f tau_T: velocity unchanged (no slip, no drag)' % ratio)
        rep.check(max(abs(x - ph.rho_p) for x in sol['var'][RHO]) <= 1e-12 * ph.rho_p,
                  't = %.2f tau_T: density unchanged' % ratio)
        rep.check(max(abs(x - ph.n) for x in sol['var'][N]) <= 1e-9 * ph.n,
                  't = %.2f tau_T: number density unchanged' % ratio)

    rep.check(worst <= 1e-5,
              'temperature matches the exact solution to %.1e of the initial gap' % worst)
    rep.close(WORK)


if __name__ == '__main__':
    main()
