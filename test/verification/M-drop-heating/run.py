"""M. Drop heating and boiling - the parabolic temperature profile inside the drop (P2T) and the boiling branch (ZGR).

A uniform cloud at rest in a uniform gas at rest (eight cells, extrapolation boundaries, drag Stokes, heat-transfer
Stokes: Nu = 2, Re = 0), so every cell integrates its sources alone, at dt = dt-max = h. The probe writes the mean
temperature T, rho_p, n and the surface temperature Ts_p1 of the centre cell. IGLOO's gas (rho 1.2, R 287, gamma 1.4,
k 0.026, mu 1.8e-5) and vapour (Lv 2.26e6, Mv 18, cp_v 1900, Le 1), rho_p0 = 1 kg/m^3, n = 1/m0.

  M1  P2T without mass transfer (IGLOO's temp-relax-p2t drop: rho 2950, c 1250, d 11.89 um, k-liq 0.026, zeta = 1),
      T0 300 K in a 600 K gas: Ts = (5 T + T_g)/6 and T - T_g = (T0 - T_g)(1 - z + z^2/2)^k, z = h/tau, tau = 1.2 tau_ITC.
  M2  ZGR, the adiabatic flash (IGLOO's boiling-zgr-flash drop: water, d 50 um, T0 640.002 K, boiling-temperature
      551.7783324 K, T_B(206 640 Pa) = 600 K = T_g): Ts = T_B, T > T_B on every row, and the flash identity
      m (Lv - c (T - T_B)) = m0 (Lv - c (T0 - T_B)) (the heat-fed mass is 3.7e-11 of m0).
  M3  ZGR, heat-fed (IGLOO's boiling-zgr-hot drop: T0 490.0002 K, boiling-temperature 422.3191381, T_B 450 K).
  M4  P2T with mass transfer (water, d 50 um, k-liq 0.6, CEM, Y_inf 0) from 300 K in an 800 K gas: the wet bulb.
  M5  IGLOO's boiling-zgr-switch drop (d 20 um, T0 640.002 K, the flash gas), run 5 ms past the step whose mean
      temperature falls below T_B: Ts = T_B above it and Ts = T below, T and m decreasing below, no invalid state.

The oracle is ICE's own time scheme on the same ODEs, typed here (a discrete twin): RK2 (Heun) in (rho_p, rho_p e)
with n fixed, the sources at Re 0 from IGLOO's equations (the P2T surface balance, CEM, ZGR with the 1/3-rule film,
alpha_s and the heat-fed rate), each root found here by bisection to 4 eps; Ts at each row's state, as the probe writes it. Gates:
on every row T and Ts within 1e-9 of the leg's gap (T0 - T_g) of the twin, m/m0 within 1e-10 (M5: 1e-5 of the gap
after the switch, whose step multiplies the difference by about 1e4). The premise, dt = dt-max on every step (every
row's time a multiple of h), is checked first. The twin's end values beside IGLOO's printed outlet values are shown,
not asserted.
"""
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import Case, Report, read_solution, try_run                     # noqa: E402

WORK = HERE / 'work'
RHO_G, R_G, GAM, K_G, MU_G = 1.2, 287.0, 1.4, 0.026, 1.8e-5
CPG = GAM * R_G / (GAM - 1.0)
RU, PATM = 8314.46, 101325.0
LV, MV, CPV, LE = 2.26e6, 18.0, 1900.0, 1.0
NU = 2.0
XSCAP = 1.0e-12
EPS4 = 4.0 * sys.float_info.epsilon


# --- the sources at Re 0, typed from IGLOO's equations ---------------------------------------------------------------

def cem(T, Tg, d, yinf, invtb):
    """CEM's rate [kg/s] (negative) at the surface temperature T: Clausius-Clapeyron, the clamp, Sh = 2."""
    p = RHO_G * R_G * Tg
    psat = PATM * math.exp(-(LV * MV / RU) * (1.0 / T - invtb))
    xs = 1.0 - XSCAP if psat >= p else min(psat / p, 1.0 - XSCAP)
    mg = RU / R_G
    ys = xs * MV / (xs * MV + (1.0 - xs) * mg)
    if ys <= yinf:
        return 0.0
    bm = (ys - yinf) / (1.0 - ys)
    return -math.pi * d * RHO_G * (K_G / (RHO_G * CPG * LE)) * 2.0 * math.log(1.0 + bm)


def alpha_s(dT):
    t = max(dT, 0.5)
    if t < 5.0:
        return 760.0 * t ** 0.26
    if t < 25.0:
        return 27.0 * t ** 2.33
    return 13800.0 * t ** 0.39


def heat_fed(gf, a, b):
    """G = Gf/x, x the root of Gf (1 + x)/x = B ln(1 + A (1 + x)), by bisection in u = ln x on [-60, 60]."""
    h = lambda u: (gf / b) * (1.0 + math.exp(-u)) - math.log1p(a * (1.0 + math.exp(u)))      # noqa: E731
    lo, hi = -60.0, 60.0
    if h(lo) <= 0.0:
        return b * math.log1p(a)
    if h(hi) >= 0.0:
        return 0.0
    while hi - lo > EPS4 * (abs(lo) + 1.0):
        mid = 0.5 * (lo + hi)
        if h(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return gf * math.exp(-0.5 * (lo + hi))


def zgr(T, TB, d, Tg):
    """(mdot, Q_s) of the boiling branch, the flash fed by alpha_s."""
    dtf = max(T - TB, 0.5)
    gf = alpha_s(T - TB) * dtf * math.pi * d * d / LV
    a = CPG * (Tg - TB) / LV
    kf = K_G * ((TB + (Tg - TB) / 3.0) / Tg) ** 0.7
    b = math.pi * kf * d * NU / CPG
    g = 0.0 if a <= 0.0 else heat_fed(gf, a, b)
    s = g + gf
    qs = s * CPG * (Tg - TB) / math.expm1(s / b) if s / b > 0.0 else b * CPG * (Tg - TB)
    return -s, qs


def root(f, a, b):
    """The root of f on [a, b] (a sign change), by bisection to 4 eps."""
    fa = f(a)
    if fa == 0.0:
        return a
    while b - a > EPS4 * max(abs(a), abs(b)):
        mid = 0.5 * (a + b)
        if mid <= a or mid >= b:
            break
        fm = f(mid)
        if fm == 0.0:
            return mid
        if (fm > 0.0) == (fa > 0.0):
            a, fa = mid, fm
        else:
            b = mid
    return 0.5 * (a + b)


class Leg(object):
    """One leg's material, gas and models, its sources and its discrete twin."""

    def __init__(self, name, rho_l, c, d0, T0, Tg, model, yinf, tboil, p2t=None, zgr=False, h=1e-5, steps=1000,
                 every=10, tokens='', physics=None):
        self.name, self.rho_l, self.c, self.d0, self.T0, self.Tg = name, rho_l, c, d0, T0, Tg
        self.model, self.yinf, self.invtb, self.kl, self.zgr = model, yinf, 1.0 / tboil, p2t, zgr
        self.h, self.steps, self.every, self.tokens, self.physics = h, steps, every, tokens, physics or {}
        self.m0 = rho_l * math.pi * d0 ** 3 / 6.0
        self.n = float('%.15E' % (1.0 / self.m0))   # as the initial condition writes it
        p = RHO_G * R_G * Tg
        self.TB = 1.0 / (self.invtb - math.log(p / PATM) / (LV * MV / RU))

    def mdot(self, T, d):
        return 0.0 if self.model == 'none' else cem(T, self.Tg, d, self.yinf, self.invtb)

    def surface(self, rho, rhoe):
        """(T, Ts, mdot, Q) at the state: the regime, the P2T balance, or the mean temperature."""
        T = rhoe / (self.c * rho + 1e-25)
        d = 2.0 * (0.75 * rho / (self.n * math.pi * self.rho_l)) ** (1.0 / 3.0)
        if self.zgr and T >= self.TB:
            md, q = zgr(T, self.TB, d, self.Tg)
            return T, self.TB, md, q
        qg = lambda t: math.pi * d * K_G * NU * (self.Tg - t)                  # noqa: E731
        if self.kl is not None:
            k = 10.0 * math.pi * self.kl * d
            ts = root(lambda t: t - T - (qg(t) + self.mdot(t, d) * LV) / k, T, self.Tg)
        else:
            ts = T
        return T, ts, self.mdot(ts, d), qg(ts)

    def source(self, rho, rhoe):
        T, ts, md, q = self.surface(rho, rhoe)
        return self.n * md, self.n * (q + md * (self.c * ts + LV))

    def twin(self, steps):
        """{step: (T, Ts, m/m0)} of the twin at the probe's steps."""
        rho, rhoe = 1.0, 1.0 * self.c * self.T0
        out = {}
        h = self.h
        for k in range(1, steps + 1):
            s0 = self.source(rho, rhoe)
            r1, e1 = rho + h * s0[0], rhoe + h * s0[1]
            s1 = self.source(r1, e1)
            rho, rhoe = rho + ((r1 - rho) + h * s1[0]) * 0.5, rhoe + ((e1 - rhoe) + h * s1[1]) * 0.5
            if k in self.want:
                T, ts, _, _ = self.surface(rho, rhoe)
                out[k] = (T, ts, rho)
        return out

    def run(self, rep, steps):
        """ICE's probe rows {step: (time, T, rho, n, Ts)} and its end field, or None when ICE failed."""
        case = Case(WORK / self.name, nx=8, Lx=1.0)
        case.particles(rho=1.0, u=0.0, v=0.0, w=0.0, T=self.T0, n=self.n)
        case.gas(rho=RHO_G, u=0.0, v=0.0, w=0.0, T=self.Tg, R=R_G, gam=GAM, k=K_G, mu=MU_G)
        case.boundaries('extrapolation')
        phys = {'latent-heat': LV, 'vapour-molar-mass': MV, 'vapour-specific-heat': CPV,
                'boiling-temperature': 1.0 / self.invtb, 'lewis-number': LE, 'vapour-mass-fraction': self.yinf}
        phys.update(self.physics)
        case.ini(t_end=steps * self.h, rk='RK2', drag='Stokes', heat='Stokes', rho_al=self.rho_l, cs=self.c,
                 dt_max=self.h, physics=phys)
        case.phase('A 1 ' + self.tokens)
        ini_set(case.dir, 'ICE-Probes', probe1='centre')
        ini_set(case.dir, 'centre', variables='T_p1 rho_p1 n_p1 Ts_p1', dtime=repr(self.every * self.h),
                index_position='1 4 1 1')
        sol = try_run(rep, case.run)
        if sol is None:
            return None, None
        # list-directed rows of five values, which a compiler may wrap over two lines
        vals = [float(x) for x in (case.dir / 'OUTPUT/centre.txt').read_text().split()]
        rows = {}
        for i in range(0, len(vals) - 4, 5):
            rows[int(round(vals[i] / self.h))] = vals[i:i + 5]
        return rows, sol


def ini_set(case_dir, section, **keys):
    """Overwrite or add keys in a section of input.ini (keys spelled with '_' for '-'); a missing section is added."""
    path = Path(case_dir) / 'input.ini'
    lines = path.read_text().split('\n')
    head = '[%s]' % section
    if head not in lines:
        lines += ['', head]
    start = lines.index(head) + 1
    end = start
    while end < len(lines) and not lines[end].startswith('['):
        end += 1
    for key, value in keys.items():
        key = key.replace('_', '-')
        for i in range(start, end):
            if lines[i].split('=')[0].strip() == key:
                lines[i] = '%s = %s' % (key, value)
                break
        else:
            lines.insert(start, '%s = %s' % (key, value))
            end += 1
    path.write_text('\n'.join(lines) + '\n')


def compare(rep, leg, rows, sol, gap, tol_T=1e-9, tol_m=1e-10, switch=None, tol_after=1e-5):
    """The premise, uniformity, n, and every row against the twin; returns the twin's rows or None."""
    print('   %s: %d probe rows' % (leg.name, len(rows)))
    off = max(abs(v[0] / leg.h - k) for k, v in rows.items())
    if not rep.check(off <= 1e-6, '%s premise: every row at a multiple of h (worst %.1e of h)' % (leg.name, off)):
        print('   %s: dt is not dt-max on every step - the twin does not apply' % leg.name)
        return None
    for iv, label in ((4, 'T'), (0, 'rho_p')):
        col = sol['var'][iv]
        rep.check(max(col) - min(col) <= 1e-12 * abs(col[0]), '%s end: %s uniform' % (leg.name, label))
    rep.check(all(x == leg.n for x in sol['var'][5]) and all(abs(v[3] / leg.n - 1.0) <= 1e-14 for v in rows.values()),
              '%s: n unchanged (exactly in the solution, to the probe\'s digits in its rows)' % leg.name)
    leg.want = set(rows)
    tw = leg.twin(max(rows))
    wT = wTs = wm = 0.0
    wT2 = wTs2 = 0.0
    for k, v in rows.items():
        T, ts, rho = tw[k]
        dT, dTs, dm = abs(v[1] - T) / gap, abs(v[4] - ts) / gap, abs(v[2] - rho)
        if switch is not None and k > switch:
            wT2, wTs2 = max(wT2, dT), max(wTs2, dTs)
        else:
            wT, wTs = max(wT, dT), max(wTs, dTs)
        wm = max(wm, dm)
    kend = max(rows)
    print('   %s: end (step %d) T %.10f K (twin %.10f), Ts %.10f K (twin %.10f), m/m0 %.10f (twin %.10f)'
          % (leg.name, kend, rows[kend][1], tw[kend][0], rows[kend][4], tw[kend][1], rows[kend][2], tw[kend][2]))
    rep.check(wT <= tol_T, '%s: T within %.0e of the gap of the twin (worst %.2e)' % (leg.name, tol_T, wT))
    rep.check(wTs <= tol_T, '%s: Ts within %.0e of the gap of the twin (worst %.2e)' % (leg.name, tol_T, wTs))
    if switch is not None:
        rep.check(wT2 <= tol_after and wTs2 <= tol_after, '%s: after the switch T and Ts within %.0e of the gap '
                  '(worst %.2e, %.2e)' % (leg.name, tol_after, wT2, wTs2))
    else:
        rep.check(wm <= tol_m, '%s: m/m0 within %.0e of the twin (worst %.2e)' % (leg.name, tol_m, wm))
    return tw


def main():
    rep = Report('M. Drop heating and boiling')
    flash = dict(rho_l=958.0, c=4186.0, Tg=600.0, model='CEM', yinf=0.0, zgr=True, physics={'evaporation': 'CEM', 'boiling': 'ZGR'},
                 tokens='')

    # M1: P2T without mass transfer
    leg = Leg('M1', 2950.0, 1250.0, 11.89e-6, 300.0, 600.0, 'none', 1.0, 373.15, p2t=0.026, steps=1000, every=10,
              physics={'evaporation': 'd2-law'}, tokens='liquid-conduction=P2T k-liq=0.026 mu-liq=1e-3')
    rows, sol = leg.run(rep, leg.steps)
    if rows is not None and compare(rep, leg, rows, sol, 300.0) is not None:
        d = 2.0 * (0.75 / (leg.n * math.pi * leg.rho_l)) ** (1.0 / 3.0)
        tau = leg.m0 * leg.c * 6.0 / (5.0 * math.pi * d * K_G * NU)
        z = leg.h / tau
        wcf = max(abs((v[1] - 600.0) - (300.0 - 600.0) * (1.0 - z + 0.5 * z * z) ** k) / 300.0 for k, v in rows.items())
        wts = max(abs(v[4] / ((5.0 * v[1] + 600.0) / 6.0) - 1.0) for v in rows.values())
        print('   M1: tau = %.6e s (1.2 tau_ITC)' % tau)
        rep.check(wcf <= 1e-9, 'M1: the Heun closed form to 1e-9 of the gap (worst %.2e)' % wcf)
        rep.check(wts <= 1e-12, 'M1: Ts = (5 T + T_g)/6 on every row to 1e-12 (worst %.2e)' % wts)

    # M2: ZGR, the adiabatic flash
    leg = Leg('M2', d0=50e-6, T0=640.002, tboil=551.7783324, h=1e-6, steps=15000, every=500, **flash)
    rows, sol = leg.run(rep, leg.steps)
    if rows is not None and compare(rep, leg, rows, sol, 40.002) is not None:
        inv0 = leg.m0 * (LV - leg.c * (leg.T0 - leg.TB))
        wid = max(abs(v[2] * leg.m0 * (LV - leg.c * (v[1] - leg.TB)) / inv0 - 1.0) for v in rows.values())
        wtb = max(abs(v[4] / leg.TB - 1.0) for v in rows.values())
        print('   M2: T_B %.9f K, the twin with the clamp ends at %.6f K (IGLOO: 604.58 K, the clamp 468.07 K)'
              % (leg.TB, clamp_end(leg)))
        rep.check(wid <= 1e-10, 'M2: the flash identity to 1e-10 (worst %.2e)' % wid)
        rep.check(wtb <= 1e-12 and all(v[1] > leg.TB for v in rows.values()),
                  'M2: Ts = T_B and T > T_B on every row (worst %.2e)' % wtb)

    # M3: ZGR, heat-fed
    leg = Leg('M3', d0=50e-6, T0=490.0002, tboil=422.3191381, h=1e-6, steps=15000, every=500, **flash)
    rows, sol = leg.run(rep, leg.steps)
    if rows is not None and compare(rep, leg, rows, sol, 40.0002) is not None:
        wtb = max(abs(v[4] / leg.TB - 1.0) for v in rows.values())
        print('   M3: T_B %.9f K, m/m0 %.7f, the twin with the clamp ends at %.6f K (IGLOO: 454.84 K, m/m0 0.829109, '
              'the clamp 381.46 K)' % (leg.TB, rows[max(rows)][2], clamp_end(leg)))
        rep.check(wtb <= 1e-12 and all(v[1] > leg.TB for v in rows.values()),
                  'M3: Ts = T_B and T > T_B on every row (worst %.2e)' % wtb)

    # M4: P2T with mass transfer
    leg = Leg('M4', 1000.0, 4180.0, 50e-6, 300.0, 800.0, 'CEM', 0.0, 373.15, p2t=0.6, steps=3000, every=100,
              physics={'evaporation': 'CEM'}, tokens='liquid-conduction=P2T k-liq=0.6 mu-liq=1e-3')
    rows, sol = leg.run(rep, leg.steps)
    if rows is not None:
        compare(rep, leg, rows, sol, 500.0)

    # M5: the switch below T_B
    leg = Leg('M5', d0=20e-6, T0=640.002, tboil=551.7783324, h=1e-6, steps=66239, every=500, **flash)
    rows, sol = leg.run(rep, leg.steps)
    if rows is not None:
        below = sorted(k for k, v in rows.items() if v[1] < leg.TB)
        above = sorted(k for k, v in rows.items() if v[1] >= leg.TB)
        rep.check(len(below) > 0 and len(above) > 0, 'M5: rows on both sides of T_B (%d above, %d below)'
                  % (len(above), len(below)))
        if below and above:
            switch = max(above)
            tw = compare(rep, leg, rows, sol, 40.002, switch=switch)
            if tw is not None:
                ok_ts = all(abs(rows[k][4] / leg.TB - 1.0) <= 1e-12 for k in above) and \
                    all(rows[k][4] == rows[k][1] for k in below)
                rep.check(ok_ts, 'M5: Ts = T_B (1e-12) above T_B, Ts = T below')
                dec = all(rows[b][1] < rows[a][1] and rows[b][2] < rows[a][2] for a, b in zip(below, below[1:]))
                rep.check(dec, 'M5: T and m decreasing below T_B')
    rep.close(WORK)


def clamp_end(leg):
    """The twin's end temperature with the boiling clamp in place of ZGR (the regime bypassed)."""
    c = Leg(leg.name, leg.rho_l, leg.c, leg.d0, leg.T0, leg.Tg, leg.model, leg.yinf, 1.0 / leg.invtb, h=leg.h,
            steps=leg.steps, every=leg.every)
    c.want = {leg.steps}
    return c.twin(leg.steps)[leg.steps][0]


if __name__ == '__main__':
    main()
