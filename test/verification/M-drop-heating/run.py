"""M. Drop heating and boiling - the parabolic temperature profile inside the drop (P2T) and the boiling branch (ZGR).

A uniform cloud at rest in a uniform gas at rest (eight cells, extrapolation boundaries, drag Stokes, heat-transfer
Stokes: Nu = 2, Re = 0), so every cell integrates its sources alone, at dt = dt-max = h. The probe writes the mean
temperature T, rho_p, n and the surface temperature Ts_p1 of the centre cell. IGLOO's gas (rho 1.2, R 287, gamma 1.4,
k 0.026, mu 1.8e-5) and vapour (Lv 2.26e6, Mv 18, cp_v 1900, Le 1), rho_p0 = 1 kg/m^3, n = 1/m0.

  M1  P2T without mass transfer (IGLOO's temp-relax-p2t drop: rho 2950, c 1250, d 11.89 um, k-liq 0.026, zeta = 1),
      T0 300 K in a 600 K gas, 15 ms: Ts = (5 T + T_g)/6 and T - T_g = (T0 - T_g)(1 - z + z^2/2)^k, z = h/tau,
      tau = 1.2 tau_ITC.
  M2  ZGR, the adiabatic flash (IGLOO's boiling-zgr-flash drop: water, d 50 um, T0 640.002 K, boiling-temperature
      551.7783324 K, T_B(206 640 Pa) = 600 K = T_g): Ts = T_B, T > T_B on every row, and the flash identity
      m (Lv - c (T - T_B)) = m0 (Lv - c (T0 - T_B)) (the heat-fed mass is 3.7e-11 of m0).
  M3  ZGR, heat-fed (IGLOO's boiling-zgr-hot drop: T0 490.0002 K, boiling-temperature 422.3191381, T_B 450 K).
  M4  P2T with mass transfer (water, d 50 um, k-liq 0.6, CEM, Y_inf 0) from 300 K in an 800 K gas: the wet bulb.
  M5  IGLOO's boiling-zgr-switch drop (d 20 um, T0 640.002 K, the flash gas), 75 ms: the mean temperature falls below
      T_B near 61 ms and reaches the wet bulb; Ts = T_B above T_B and Ts = T below, m decreasing and T not increasing
      below, no invalid state.
  M5c M5 at dt-max 1e-4 s, which does not resolve the switch at T_B (the step lands in the clamp below T_B, whose
      evaporation removes the drop mass in less than a step): ICE warns, once; M2, M3 and M5 carry no warning.

The oracle is ICE's own time scheme on the same ODEs, typed here (a discrete twin): RK2 (Heun) in (rho_p, rho_p e)
with n fixed, the sources at Re 0 from IGLOO's equations (the P2T surface balance, CEM, ZGR with the 1/3-rule film,
alpha_s and the heat-fed rate), each root found here by bisection to 4 eps; Ts at each row's state, as the probe writes it. Gates:
on every row T and Ts within 1e-9 of the leg's gap (T0 - T_g) of the twin, m/m0 within 1e-10 (M5: 1e-5 of the gap
after the switch, whose step multiplies the difference by about 1e4). The premise, dt = dt-max on every step (every
row's time a multiple of h), is checked first.

The ODEs themselves are checked against IGLOO's trajectory rows of the same drops (its cases temp-relax-p2t,
boiling-zgr-flash, boiling-zgr-hot, boiling-zgr-switch at t = x/u_g; M4 has no IGLOO case): ICE's value at that time
(the twin carries ICE's nearest row to it) within IGLOO's own row tolerance against its oracle (the F12.6/E13.6 prints
and its integrator) plus ICE's time-step error, 2 |twin(h) - twin(h/2)| + h |dX/dt| times the steps of alpha_s at 5
and 25 K (measured Richardson ratios 4.0 on M1, 0.67 on M2 where a step of alpha_s is crossed, 22 on M3, 2.3 on M5
past the switch), plus |ICE - twin| at the row. A constant mistyped in this driver, which ICE's input and the twin
share, passes the twin and fails here.
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

# IGLOO develop's trajectory rows at OMP 1 (25 identical parcels): leg -> [(case, u_g, x, T, m, Ts,
# (T, m) of the row 0.0025 m upstream, (a, b) of its T tolerance a + b |dT/dx|)]
IGLOO_DX = 0.0025
IGLOO_ROWS = {
    'M1': [('temp-relax-p2t', 10.0, 0.15, 599.830919, None, 599.859099, (599.808466, None), (1.01e-6, 1e-6))],
    'M2': [('boiling-zgr-flash', 10.0, 0.15, 604.575815, 0.585515e-10, None, (604.615216, 0.585559e-10), (1.5e-6, 5e-7))],
    'M3': [('boiling-zgr-hot', 10.0, 0.15, 454.841467, 0.519859e-10, None, (454.875281, 0.520986e-10), (1.5e-6, 5e-7))],
    'M5': [('boiling-zgr-switch', 1.0, 0.06, 600.030155, 0.371575e-11, None, (600.091038, 0.371617e-11), (1.5e-6, 5e-7)),
           ('boiling-zgr-switch', 1.0, 0.075, 463.825379, 0.908282e-12, None, (463.825391, 0.122417e-11), (1.5e-6, 5e-7))],
}
IGLOO_TS_TOL = 1e-6


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


def compare(rep, leg, rows, sol, gap, tol_T=1e-9, tol_m=1e-10, switch=None, tol_after=1e-5, extra=()):
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
    leg.want = set(rows) | set(extra)
    tw = leg.twin(max(leg.want))
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


WARNING = '[WARNING] boiling = ZGR'


def warnings(leg):
    """The lines of ICE's log that warn of an unresolved step under ZGR."""
    return [l for l in (WORK / leg.name / 'log').read_text().split('\n') if WARNING in l]


def igloo_steps(leg):
    return [int(round(r[2] / r[1] / leg.h)) for r in IGLOO_ROWS[leg.name]]


def e13_half_ulp(v):
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def vs_igloo(rep, leg, rows, tw):
    """ICE at IGLOO's row times against IGLOO's printed T, m and Ts, within IGLOO's row tolerance plus ICE's step error."""
    ks = igloo_steps(leg)
    half = Leg.__new__(Leg)
    half.__dict__.update(leg.__dict__)
    half.h, half.want = 0.5 * leg.h, set(2 * k for k in ks)
    th = half.twin(max(half.want))
    jump = 0.0
    if leg.zgr:
        jump = sum(abs(alpha_s(t * (1.0 - 1e-12)) / alpha_s(t) - 1.0) for t in (5.0, 25.0))
    for (case, ug, x, Ti, mi, tsi, prev, tolab), k in zip(IGLOO_ROWS[leg.name], ks):
        kr = min(rows, key=lambda j: abs(j - k))
        T0, ts0, r0 = tw[kr]
        T1, ts1, r1 = tw[k]
        T2, ts2, r2 = th[2 * k]
        rT, rts, rm = rows[kr][1] - T0, rows[kr][4] - ts0, (rows[kr][2] - r0) / leg.n
        dTdx = abs(Ti - prev[0]) / IGLOO_DX
        tolT = tolab[0] + tolab[1] * dTdx + 2.0 * abs(T2 - T1) + leg.h * jump * dTdx * ug + abs(rT)
        eT = abs(T1 + rT - Ti)
        lab = '%s vs IGLOO %s at t = %g s' % (leg.name, case, x / ug)
        print('   %s (step %d from row %d): T %.6f K, m %.7e kg, Ts %.6f K' % (lab, k, kr, T1 + rT, (r1 + rm * leg.n) / leg.n,
                                                                          ts1 + rts))
        rep.check(eT <= tolT, '%s: T within %.2e K of %.6f (off %.2e, %.2f of it)' % (lab, tolT, Ti, eT, eT / tolT))
        if mi is not None:
            m1 = r1 / leg.n + rm
            dmdx = abs(mi - prev[1]) / IGLOO_DX
            tolm = e13_half_ulp(mi) + 5e-7 * dmdx + 1e-8 * mi + 2.0 * abs(r2 - r1) / leg.n \
                + leg.h * jump * dmdx * ug + abs(rm)
            em = abs(m1 - mi)
            rep.check(em <= tolm, '%s: m within %.2e kg of %.6e (off %.2e, %.2f of it; the E13.6 half-ULP %.0e)'
                      % (lab, tolm, mi, em, em / tolm, e13_half_ulp(mi)))
        if tsi is not None:
            tolts = IGLOO_TS_TOL + 2.0 * abs(ts2 - ts1) + abs(rts)
            ets = abs(ts1 + rts - tsi)
            rep.check(ets <= tolts, '%s: Ts within %.2e K of %.6f (off %.2e, %.2f of it)'
                      % (lab, tolts, tsi, ets, ets / tolts))


def main():
    rep = Report('M. Drop heating and boiling')
    flash = dict(rho_l=958.0, c=4186.0, Tg=600.0, model='CEM', yinf=0.0, zgr=True, physics={'evaporation': 'CEM', 'boiling': 'ZGR'},
                 tokens='')

    # M1: P2T without mass transfer
    leg = Leg('M1', 2950.0, 1250.0, 11.89e-6, 300.0, 600.0, 'none', 1.0, 373.15, p2t=0.026, steps=1500, every=10,
              physics={'evaporation': 'd2-law'}, tokens='liquid-conduction=P2T k-liq=0.026 mu-liq=1e-3')
    rows, sol = leg.run(rep, leg.steps)
    tw = None if rows is None else compare(rep, leg, rows, sol, 300.0, extra=igloo_steps(leg))
    if tw is not None:
        vs_igloo(rep, leg, rows, tw)
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
    tw = None if rows is None else compare(rep, leg, rows, sol, 40.002, extra=igloo_steps(leg))
    if tw is not None:
        vs_igloo(rep, leg, rows, tw)
        inv0 = leg.m0 * (LV - leg.c * (leg.T0 - leg.TB))
        wid = max(abs(v[2] * leg.m0 * (LV - leg.c * (v[1] - leg.TB)) / inv0 - 1.0) for v in rows.values())
        wtb = max(abs(v[4] / leg.TB - 1.0) for v in rows.values())
        print('   M2: T_B %.9f K, the twin with the clamp ends at %.6f K (IGLOO\'s clamp 468.07 K)'
              % (leg.TB, clamp_end(leg)))
        rep.check(wid <= 1e-10, 'M2: the flash identity to 1e-10 (worst %.2e)' % wid)
        rep.check(wtb <= 1e-12 and all(v[1] > leg.TB for v in rows.values()),
                  'M2: Ts = T_B and T > T_B on every row (worst %.2e)' % wtb)

    # M3: ZGR, heat-fed
    leg = Leg('M3', d0=50e-6, T0=490.0002, tboil=422.3191381, h=1e-6, steps=15000, every=500, **flash)
    rows, sol = leg.run(rep, leg.steps)
    tw = None if rows is None else compare(rep, leg, rows, sol, 40.0002, extra=igloo_steps(leg))
    if tw is not None:
        vs_igloo(rep, leg, rows, tw)
        wtb = max(abs(v[4] / leg.TB - 1.0) for v in rows.values())
        print('   M3: T_B %.9f K, m/m0 %.7f, the twin with the clamp ends at %.6f K (IGLOO\'s clamp 381.46 K)'
              % (leg.TB, rows[max(rows)][2], clamp_end(leg)))
        rep.check(wtb <= 1e-12 and all(v[1] > leg.TB for v in rows.values()),
                  'M3: Ts = T_B and T > T_B on every row (worst %.2e)' % wtb)

    # M4: P2T with mass transfer
    leg = Leg('M4', 1000.0, 4180.0, 50e-6, 300.0, 800.0, 'CEM', 0.0, 373.15, p2t=0.6, steps=3000, every=100,
              physics={'evaporation': 'CEM'}, tokens='liquid-conduction=P2T k-liq=0.6 mu-liq=1e-3')
    rows, sol = leg.run(rep, leg.steps)
    if rows is not None:
        compare(rep, leg, rows, sol, 500.0)

    # M5: the switch below T_B
    leg = Leg('M5', d0=20e-6, T0=640.002, tboil=551.7783324, h=1e-6, steps=75000, every=500, **flash)
    rows, sol = leg.run(rep, leg.steps)
    if rows is not None:
        below = sorted(k for k, v in rows.items() if v[1] < leg.TB)
        above = sorted(k for k, v in rows.items() if v[1] >= leg.TB)
        rep.check(len(below) > 0 and len(above) > 0, 'M5: rows on both sides of T_B (%d above, %d below)'
                  % (len(above), len(below)))
        if below and above:
            switch = max(above)
            tw = compare(rep, leg, rows, sol, 40.002, switch=switch, extra=igloo_steps(leg))
            if tw is not None:
                vs_igloo(rep, leg, rows, tw)
                ok_ts = all(abs(rows[k][4] / leg.TB - 1.0) <= 1e-12 for k in above) and \
                    all(rows[k][4] == rows[k][1] for k in below)
                rep.check(ok_ts, 'M5: Ts = T_B (1e-12) above T_B, Ts = T below')
                dec = all(rows[b][1] <= rows[a][1] and rows[b][2] < rows[a][2] for a, b in zip(below, below[1:]))
                rep.check(dec, 'M5: m decreasing and T not increasing below T_B')

    # M5c: the switch unresolved at dt-max 1e-4 s, and the warning
    for name in ('M2', 'M3', 'M5'):
        if (WORK / name / 'log').exists():
            rep.check(not warnings(Leg(name, d0=20e-6, T0=640.002, tboil=551.7783324, **flash)),
                      '%s: no warning of an unresolved step' % name)
    leg = Leg('M5c', d0=20e-6, T0=640.002, tboil=551.7783324, h=1e-4, steps=750, every=50, **flash)
    rows, sol = leg.run(rep, leg.steps)
    if (WORK / leg.name / 'log').exists():
        w = warnings(leg)
        print('   M5c: ' + (w[0].strip() if w else 'no warning'))
        rep.check(len(w) == 1, 'M5c: one warning of an unresolved step (%d)' % len(w))
    rep.close(WORK)


def clamp_end(leg):
    """The twin's end temperature with the boiling clamp in place of ZGR (the regime bypassed)."""
    c = Leg(leg.name, leg.rho_l, leg.c, leg.d0, leg.T0, leg.Tg, leg.model, leg.yinf, 1.0 / leg.invtb, h=leg.h,
            steps=leg.steps, every=leg.every)
    c.want = {leg.steps}
    return c.twin(leg.steps)[leg.steps][0]


if __name__ == '__main__':
    main()
