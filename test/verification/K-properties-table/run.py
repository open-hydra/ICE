"""K. The condensed-material property table - exact solution at a fixed temperature.

Case A's cloud (at rest in a uniform stream, Stokes drag) with heat-transfer = NoHeat
and emissivity = 0: nothing heats the particles and the drag work cancels the kinetic
energy it produces, so T_p keeps its initial value and the table is read at one fixed
temperature for the whole run. The velocity is then case A's exponential,

    u_p(t) = u_g (1 - exp(-t/tau_p)),   tau_p = rho_mat dp^2 / (18 mu),

with rho_mat the table's density at T_p. ICE does not carry the radius: it recovers
Rp = (3 rho_p / (4 pi n rho_mat))^(1/3) from the state, so at fixed (rho_p, n) tau_p
goes as rho_mat^(1/3). rho_p and n are set so that Rp = dp/2 at the rho_mat the oracle
expects, which it reads from the table as written: linear between the integer-kelvin
rows, saturated past either end. The table owns density and cp, so the INI sets
neither (common.Case.ini).

  K1   interpolation. Rows 1..5000 K, cp constant,
       rho 2000 up to 300 K and 1000 from 301 K, T_p = 300.4 K: linear gives 1600.
       The nearest kelvin reads 2000, tau x (2000/1600)^(1/3) = 1.077 and u(tau) is
       4.3 % low. A fixed T_p avoids the cancellation a temperature sweep would give.
  K1c  non-vacuity: K1's table with rho = 1600 on every row; any reader passes it.
  K6   range: a constant table on 280..400 K (rho 1500, cp 2000), T_p = 300.4 K.
  K7   out of range: K6's rows with rho = 1500 + 5 (T - 280), T_p = 250 K below Tmin.
       The table saturates: rho(Tmin) = 1500. A linear extrapolation (1350) would put
       u(tau) 2.0 % high, the row at Tmax (2100) 6.5 % low.

Tolerances. With dt-max = tau/1000, RK2 leaves e^-1 (dt/tau)^2 / 6 = 6e-8 u_g on u(tau),
1e-7 relative (the same formula gives case A's measured 6.4e-7 u_g at 309 steps per tau):
the 1e-4 tolerance sits three decades above it and 2.6 below K1's signal. T_p drifts by
round-off only (case A holds it to 1e-9 relative); 1e-4 K is far above that and far below
the 0.1 K that would move nint(300.4) to the next row.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, table_lookup, try_run    # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)
TOL_U, TOL_T = 1.0e-4, 1.0e-4
CP = 2000.0


def leg(rep, label, Tp, table):
    """Case A's cloud at T_p = T_g = Tp with its own table, run to one tau: u_p(t)
    against the closed form at the table's linear density, T_p fixed."""
    case = Case(WORK / label.lower(), nx=8, Lx=1.0)
    written = case.properties(**table)
    rho_mat = table_lookup(written, 'rho', Tp)
    near = min(max(int(math.floor(Tp + 0.5)), written['T'][0]), written['T'][-1])
    rho_near = written['rho'][int(near - written['T'][0])]
    ph = Physics(rho_al=rho_mat, Tg=Tp)
    tau = ph.tau_stokes
    case.particles(rho=ph.rho_p, u=0.0, v=0.0, w=0.0, T=Tp, n=ph.n)
    case.gas(rho=ph.rho_g, u=ph.ug, v=0.0, w=0.0, T=Tp,
             R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=tau, drag='Stokes', heat='NoHeat', dt_max=tau / 1000.0)
    print('   %-3s rows %d..%d K, T_p = %.1f K: rho_mat %.3f (linear, saturated), '
          'nearest row %.3f' % (label, written['T'][0], written['T'][-1], Tp, rho_mat, rho_near))

    sol = try_run(rep, case.run)
    if sol is None:
        return
    for line in (case.dir / 'log').read_text().splitlines():
        if 'condensed properties' in line:
            print('       reader: %s' % line.strip())
    t, u, temp = sol['time'], sol['var'][U][0], sol['var'][T]
    exact = ph.ug * (1.0 - math.exp(-t / tau))
    near_u = ph.ug * (1.0 - math.exp(-t / (tau * (rho_near / rho_mat) ** (1.0 / 3.0))))
    err = (u - exact) / exact
    seen = rho_mat * (-t / math.log(1.0 - u / ph.ug) / tau) ** 3 if 0.0 < u < ph.ug else float('nan')
    print('       t = %.4f tau: u_p/u_g ICE %.9f, exact %.9f, rel. error %+.3e '
          '(nearest-row prediction %+.3e)' % (t / tau, u / ph.ug, exact / ph.ug, err,
                                               (near_u - exact) / exact))
    print('       density ICE applied, inverted from u_p(t): %.3f kg/m^3' % seen)
    rep.check(abs(err) <= TOL_U, '%s: u_p(t) matches the closed form at rho_mat = %.1f to %.1e'
              % (label, rho_mat, abs(err)))
    drift = max(abs(x - Tp) for x in temp)
    rep.check(drift <= TOL_T, '%s: T_p stays at %.1f K to %.1e K' % (label, Tp, drift))


def main():
    rep = Report('K. Property table (fixed T)')
    step = lambda T: 2000.0 if T <= 300 else 1000.0          # noqa: E731

    leg(rep, 'K1', 300.4, dict(Tmin=1, Tmax=5000, cp=CP, rho=step))
    leg(rep, 'K1c', 300.4, dict(Tmin=1, Tmax=5000, cp=CP, rho=1600.0))
    leg(rep, 'K6', 300.4, dict(Tmin=280, Tmax=400, cp=CP, rho=1500.0))
    leg(rep, 'K7', 250.0, dict(Tmin=280, Tmax=400, cp=CP, rho=lambda T: 1500.0 + 5.0 * (T - 280)))

    rep.close(WORK)


if __name__ == '__main__':
    main()
