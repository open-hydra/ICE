"""K. The condensed-material property table - exact solutions at a fixed and a varying T.

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

A varying cp. The energy of a tabulated cp is rho_p e(T) with e = h - hOff, linear between
the rows, hOff the enthalpy at 0 K along the first segment; T comes back from e by
inverting that line. cp = 1000 + 2 T on rows 1..5000 K, h its exact integral (integers).

  K2   stationarity. No gas: a uniform cloud moving at 10 m/s at T_p = 500.37 K (between
       two rows, so e(T) and its inverse both interpolate), 200 RK2 steps, nothing acts on
       it. Each step's round trips T -> e -> T must return T_p, |T_p - T_p0| <= 1e-9 K.
       The state rho_p cp(T) T, recovered in two passes seeded by specific-heat (1598, the
       table forbids setting it), walks 500 -> 444.1 -> 409.2 -> 385.3 K per round trip
       towards 299 K.
  K2c  non-vacuity: K2 with cp = 2000 on every row, which any energy state keeps.
  K3   heating. The cloud moves with a gas at 800 K (no slip, drag = NoDrag), T_p0 = 300 K,
       heat = Stokes (Nu = 2): rho_p de/dt = 4 pi kg Rp n (T_g - T), with T(e) the table's
       line. The oracle integrates it with RK4 in e. The state rho_p cp(T) T heats with the
       capacity cp + T dcp/dT instead of cp and loses 56 K per round trip on top.
  K3c  control: K3 with cp = 2000 on every row against the same oracle.
  K4   datum: K3's table as Enthalpy_abs with h - 1.5e7: part-field.tec bit for bit K3's.
       e = h - hOff, and on integer rows both subtractions are exact.
  K5   permuted columns: K3's table written Temperature, Enthalpy, Density, Cp: part-field.tec
       bit for bit K3's.

Tolerances (varying cp). K2: a round trip returns T to a few ulp (1e-13 K at 500 K); 400
of them stay below 1e-10 K, so 1e-9 K is a decade above the worst accumulation and eleven
decades below the RED signal. K3: RK2 at dt-max = tau_T/1000 leaves (dt/tau)^2/6 e^-1 =
6e-8 of the 500 K gap at t = tau_T (measured: see the printed error, the constant-cp
control shares it); the kinks of e at the rows add a first-order term per crossed row,
(dt/tau)^2 * 2/cp ~ 1e-9 each over ~300 rows. 1e-6 of the gap sits a decade above both
and five below the signal (the capacity error alone is 25 %).
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, table_lookup, try_run    # noqa: E402
import bisect                                                       # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)
TOL_U, TOL_T = 1.0e-4, 1.0e-4
TOL_STATIONARY, TOL_HEAT = 1.0e-9, 1.0e-6
CP = 2000.0
VARYING = dict(Tmin=1, Tmax=5000, cp=lambda T: 1000.0 + 2.0 * T, rho=1000.0)


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


class Energy(object):
    """The oracle's e(T) and T(e): e = h - hOff, linear between the rows as written, the first
    segment down to e(0) = 0, the last one beyond Tmax."""

    def __init__(self, table):
        Ts, h, cp = table['T'], table['h'], table['cp']
        if all(c == cp[0] for c in cp):
            off = h[0] - cp[0] * Ts[0]
        else:
            off = h[0] - Ts[0] * (h[1] - h[0])
        self.T, self.e = Ts, [x - off for x in h]

    def of(self, T):
        Ts, e = self.T, self.e
        if T < Ts[0]:
            return e[0] * T / Ts[0]
        i = min(max(bisect.bisect_right(Ts, T) - 1, 0), len(Ts) - 2)
        return e[i] + (e[i + 1] - e[i]) * (T - Ts[i]) / (Ts[i + 1] - Ts[i])

    def T_of(self, x):
        Ts, e = self.T, self.e
        if x < e[0]:
            return Ts[0] * x / e[0]
        i = min(max(bisect.bisect_right(e, x) - 1, 0), len(e) - 2)
        return Ts[i] + (x - e[i]) / (e[i + 1] - e[i]) * (Ts[i + 1] - Ts[i])


def permute_columns(path, order):
    """Rewrite a written table with its value columns in `order` (header names), bytes kept."""
    lines = Path(path).read_text().splitlines()
    names = [n.strip().strip('"') for n in lines[1].split('=', 1)[1].split(',')]
    pick = [names.index(n) for n in ['Temperature'] + order]
    lines[1] = 'VARIABLES = ' + ', '.join('"%s"' % names[k] for k in pick)
    for i in range(4, len(lines)):
        f = lines[i].split()
        lines[i] = ' '.join(f[k] for k in pick)
    Path(path).write_text('\n'.join(lines) + '\n')


def stationary(rep, label, table, Tp=500.37):
    """No gas: a uniform cloud moving at 10 m/s at Tp, 200 steps. T_p must stay at Tp."""
    ph = Physics()
    case = Case(WORK / label.lower(), nx=8, Lx=1.0)
    case.properties(**table)
    case.particles(rho=ph.rho_p, u=ph.ug, v=0.0, w=0.0, T=Tp, n=ph.n)
    case.boundaries('extrapolation')
    case.ini(iters=200, drag='NoDrag', heat='NoHeat')
    sol = try_run(rep, case.run)
    if sol is None:
        return
    drift = max(abs(x - Tp) for x in sol['var'][T])
    print('   %-3s no gas, 200 steps at T_p = %.2f K: max|T_p - T_p0| = %.3e K' % (label, Tp, drift))
    rep.check(drift <= TOL_STATIONARY, '%s: T_p stays at %.2f K to %.1e K' % (label, Tp, drift))


def heating(label, table, **kw):
    """A cloud moving with a gas at 800 K, T_p0 = 300 K, Nu = 2: T_p(t) against the RK4
    integration of rho_p de/dt = 4 pi kg Rp n (T_g - T(e)). Writes the case, runs nothing."""
    ph = Physics(Tg=800.0)
    case = Case(WORK / label.lower(), nx=8, Lx=1.0)
    written = case.properties(**dict(table, **kw))
    energy = Energy(written)
    tau = ph.rho_al * 2000.0 * ph.dp ** 2 / (12.0 * ph.kg)
    case.particles(rho=ph.rho_p, u=ph.ug, v=0.0, w=0.0, T=300.0, n=ph.n)
    case.gas(rho=ph.rho_g, u=ph.ug, v=0.0, w=0.0, T=ph.Tg, R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=tau, drag='NoDrag', heat='Stokes', dt_max=tau / 1000.0)
    return case, energy, ph, tau


def run_heating(rep, label, case, energy, ph, tau):
    """Run a heating() case against its oracle. Returns part-field.tec's bytes."""
    sol = try_run(rep, case.run)
    if sol is None:
        return None
    t = sol['time']
    rate = 4.0 * math.pi * ph.kg * ph.rp * ph.n / ph.rho_p
    steps = 20000
    dt = t / steps
    e = energy.of(300.0)
    f = lambda x: rate * (ph.Tg - energy.T_of(x))                    # noqa: E731
    for _ in range(steps):
        k1 = f(e); k2 = f(e + 0.5 * dt * k1); k3 = f(e + 0.5 * dt * k2); k4 = f(e + dt * k3)
        e += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    ref = energy.T_of(e)
    temp = sol['var'][T]
    err = max(abs(x - ref) for x in temp) / (ph.Tg - 300.0)
    print('   %-3s t = %.4f tau_T: T_p ICE %.9f K, oracle %.9f K, error %.3e of the 500 K gap'
          % (label, t / tau, temp[0], ref, err))
    rep.check(err <= TOL_HEAT, '%s: T_p(t) matches the oracle to %.1e of the gap' % (label, err))
    return (case.dir / 'OUTPUT/part-field.tec').read_bytes()


def main():
    rep = Report('K. Property table')
    step = lambda T: 2000.0 if T <= 300 else 1000.0          # noqa: E731

    leg(rep, 'K1', 300.4, dict(Tmin=1, Tmax=5000, cp=CP, rho=step))
    leg(rep, 'K1c', 300.4, dict(Tmin=1, Tmax=5000, cp=CP, rho=1600.0))
    leg(rep, 'K6', 300.4, dict(Tmin=280, Tmax=400, cp=CP, rho=1500.0))
    leg(rep, 'K7', 250.0, dict(Tmin=280, Tmax=400, cp=CP, rho=lambda T: 1500.0 + 5.0 * (T - 280)))

    stationary(rep, 'K2', VARYING)
    stationary(rep, 'K2c', dict(VARYING, cp=CP))

    k3 = run_heating(rep, 'K3', *heating('K3', VARYING))
    run_heating(rep, 'K3c', *heating('K3c', dict(VARYING, cp=CP)))
    k4 = run_heating(rep, 'K4', *heating('K4', VARYING, datum='Enthalpy_abs', h0=-1.5e7))
    case, energy, ph, tau = heating('K5', VARYING)
    permute_columns(case.dir / 'INPUT/part-properties.dat', ['Enthalpy', 'Density', 'Cp'])
    k5 = run_heating(rep, 'K5', case, energy, ph, tau)
    if k3 is not None:
        rep.check(k4 == k3, 'K4: an absolute datum gives part-field.tec bit for bit')
        rep.check(k5 == k3, 'K5: permuted columns give part-field.tec bit for bit')

    rep.close(WORK)


if __name__ == '__main__':
    main()
