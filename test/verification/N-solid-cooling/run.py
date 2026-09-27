"""N. Solidification in a closed cell - supercooling, recalescence, plateau, solid, and melting at T-melt.

IGLOO's solid-box material (tests/solidification/solid-box): droplets of d = 30 um, rho 2950 kg/m^3, c_l 1250,
c_s 600 J/(kg K), h_fus 1.07e6 J/kg, T_m 2327 K, T_nuc = 0.8 T_m = 1861.6 K, phase line
"A 1 solidification=on h-fus=1.07e6 cp-solid=600" (T-melt and T-nuc at their defaults). Eight uniform cells move with
the gas (no slip, drag = NoDrag), heat-transfer = Stokes (Nu = 2 exactly), emissivity 0, extrapolation boundaries, so
every cell integrates the heat source alone: m = rho pi d^3/6, A = pi d k_g Nu, tau = m c/A, and

  liquid      T = T_g + (T_0 - T_g) exp(-t/tau_l)               until T = T_nuc at t_n
  nucleation  f0 = c_l (T_m - T_nuc)/h_fus: plateau at T_m with f = f0 (f0 < 1), or the solid at
              T_m - (c_l (T_m - T_nuc) - h_fus)/c_s with f = 1 (the energy-exact jump)
  plateau     T = T_m, h_fus df/dt = A (T_m - T_g)/m, freezing; melting when the gas is hotter
  solid       T = T_g + (T_s - T_g) exp(-(t - t_s)/tau_s)       up to T_m when heated, which starts the plateau
  melted      f = 0 at T_m, then the liquid again, with chi cleared

The nucleated fraction chi is 1 on every nucleated state and 0 on every liquid one (a closed cell has no mixing).

  N1  cooling from 2400 K in a 600 K gas: the liquid, three plateau times, two solid times.
  N2  h-fus = 4e5 (f0 = 1.4544 > 1): the jump lands on the solid at 2024.083 K, then tau_s.
  N3  a solid at 1500 K in a 3000 K gas: it reaches T_m, stays there while f falls by the heat over h_fus.
  N4  a mush (T_m, f = 0.8, chi = 1) in a 3000 K gas: f falls to 0, then a liquid heating, chi = 0.
  N5  N1's plateau midpoint and 0.015 s, and N3, with the IG and AG closures (P = 1e-6).
  N6  N1 to the plateau midpoint, then a restart in place (newrun = false) to 0.015 s: the straight run's state.
  N7  two materials, "A 1 solidification=on ..." and "B 1" (cp 1000): each its own history; B's columns are
      those of a run with B alone, byte for byte.

Tolerances (time-scheme RK2, dt-max = 1e-5 binding). The liquid and melting legs have no event: 1e-5 of the gap,
|f - f_cf| <= 1e-5. After a nucleation the switch happens at the end of a stage: the plateau is shifted by at most
pi d k_g Nu (T_m - T_n) dt/(m h_fus) = 51.12 dt, the solid by (T_m - T_n)/(T_m - T_g) dt = 0.2695 dt; the plateau
holds T = T_m to 1e-9 K and f to 51.12 dt + 1e-6, the solid T to (T - T_g)/tau_s 0.2695 dt 1.2 + 0.018 K.
"""
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import Case, Report, write_ini, run_ice, read_solution            # noqa: E402

WORK = HERE / 'work'

# --- IGLOO's solid-box material and gas --------------------------------------------------
RHO_M, CL, CS, HF, TM = 2950.0, 1250.0, 600.0, 1.07e6, 2327.0
TN = 0.8 * TM
DP, KG, NU = 30.0e-6, 0.026, 2.0
RP = 0.5 * DP
M = RHO_M * math.pi / 6.0 * DP ** 3
A = math.pi * DP * KG * NU
TAU_L, TAU_S = M * CL / A, M * CS / A
RHO_G, UG, R_G, GAM, MU = 1.2, 10.0, 287.05, 1.4, 1.8e-5
ALPHA = 1.0e-3
RHO_P = ALPHA * RHO_M
NDENS = RHO_P / (RHO_M * 4.0 / 3.0 * math.pi * RP ** 3)
DT = 1.0e-5
PHASE = 'A 1 solidification=on h-fus=1.07e6 cp-solid=600'
LAYOUT = {'MK': ('rp', 'up', 'vp', 'wp', 'Tp', 'np'),
          'IG': ('rp', 'up', 'vp', 'wp', 'Pp', 'Tp', 'np'),
          'AG': ('rp', 'up', 'vp', 'wp', 'P11', 'P12', 'P13', 'P22', 'P23', 'P33', 'Tp', 'np')}
PSEUDO = 1.0e-6


# ---------------------------------------------------------------------------
#  Closed forms (0-D): T and f at time t
# ---------------------------------------------------------------------------

def cooling(t, T0, Tg, h=HF):
    """A liquid from T0 in a colder gas: supercooling, the jump, the plateau (or the hypercooled solid), the solid."""
    tn = TAU_L * math.log((T0 - Tg) / (TN - Tg))
    if t < tn:
        return Tg + (T0 - Tg) * math.exp(-t / TAU_L), 0.0
    f0 = CL * (TM - TN) / h
    if f0 >= 1.0:
        Ta = TM - (CL * (TM - TN) - h) / CS
        return Tg + (Ta - Tg) * math.exp(-(t - tn) / TAU_S), 1.0
    rate = A * (TM - Tg) / (M * h)
    ts = tn + (1.0 - f0) / rate
    if t < ts:
        return TM, f0 + (t - tn) * rate
    return Tg + (TM - Tg) * math.exp(-(t - ts) / TAU_S), 1.0


def heating(t, T0, Tg, h=HF, f_start=None):
    """A solid from T0 (or a mush at T_m with f_start) in a hotter gas: to T_m, melting there, then the liquid."""
    rate = A * (Tg - TM) / (M * h)
    if f_start is None:
        tm = TAU_S * math.log((Tg - T0) / (Tg - TM))
        if t < tm:
            return Tg - (Tg - T0) * math.exp(-t / TAU_S), 1.0
        f_start, t0 = 1.0, tm
    else:
        t0 = 0.0
    tl = t0 + f_start / rate
    if t < tl:
        return TM, f_start - (t - t0) * rate
    return Tg - (Tg - TM) * math.exp(-(t - tl) / TAU_L), 0.0


# ---------------------------------------------------------------------------
#  Runs
# ---------------------------------------------------------------------------

def attempt(rep, label, case):
    """case.run, or one FAIL line quoting the refusal the log holds (decision 8: the reason is read here). The columns
    are read by position, so a solidifying first material must write f_p1 and chi_p1 right after n_p1."""
    try:
        sol = case.run()
    except RuntimeError as exc:
        log = (case.dir / 'log').read_text(errors='replace').splitlines() if (case.dir / 'log').exists() else []
        hits = [s.strip() for s in log if 'Wrong solidification input' in s or '[ERROR]' in s]
        rep.check(False, '%s: %s%s' % (label, str(exc).split(' in ')[0], (': ' + hits[0]) if hits else ''))
        return None
    if 'solidification=on' in (case.dir / 'INPUT/part-phase.txt').read_text().splitlines()[1]:
        names = sol['names']
        at = names.index('n_p1') + 1 if 'n_p1' in names else len(names)
        if names[at:at + 2] != ['f_p1', 'chi_p1']:
            rep.check(False, '%s: the output does not name f_p1, chi_p1 after n_p1 (%s)' % (label, ' '.join(names)))
            return None
    return sol


def ic_state(closure, T, f=None, chi=None):
    """The IC fields of one family in its closure's layout, with f and chi appended when given."""
    P = dict(IG=[PSEUDO], AG=[PSEUDO, 0.0, 0.0, PSEUDO, 0.0, PSEUDO]).get(closure, [])
    fields = [RHO_P, UG, 0.0, 0.0] + P + [T, NDENS]
    if f is not None:
        fields += [f, chi]
    return fields


def closed_cell(label, T0, Tg, t_end, closure='MK', phase=(PHASE,), full=None, zones=None, extra=()):
    """Eight uniform cells; full = (f, chi) writes the IC with the fractions (the full layout)."""
    case = Case(WORK / label.lower(), nx=8, Lx=1.0)
    if zones is not None:
        case.properties_zones(250, 3000, zones)
    names, fields = [], []
    nfam = len(phase) + len(extra)
    for p in range(1, nfam + 1):
        cl = closure if p == 1 else 'MK'
        base = LAYOUT[cl]
        f_chi = full if (p == 1 and full is not None) else (None, None)
        names += ['%s%d' % (s, p) for s in base] + (['fp%d' % p, 'chip%d' % p] if f_chi[0] is not None else [])
        fields += ic_state(cl, T0, *f_chi)
    case.write_ic(names, fields)
    case.gas(rho=RHO_G, u=UG, v=0.0, w=0.0, T=Tg, R=R_G, gam=GAM, k=KG, mu=MU)
    case.boundaries('extrapolation')
    case.phase(*(tuple(phase) + tuple(extra)))
    physics = {'emissivity': ' '.join(['0.0'] * nfam)} if nfam > 1 else None
    case.ini(t_end=t_end, rk='RK2', drag='NoDrag', heat='Stokes', rho_al=RHO_M, cs=CL, dt_max=DT,
             physics=physics, closures=(closure,) + ('MK',) * (nfam - 1))
    return case


def read_ini(path):
    """{section: {key: value}} of an input.ini, values as written."""
    sections, cur = {}, None
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if s.startswith('[') and s.endswith(']'):
            cur = sections.setdefault(s[1:-1], {})
        elif '=' in s and cur is not None:
            k, v = s.split('=', 1)
            cur[k.strip()] = v.strip()
    return sections


def restart_in_place(case, parameters, gas=None):
    """Rewrite the case's input.ini with write_ini from its own sections plus the [ICE-Parameters] given (newrun =
    false among them), and its gas file when gas is given (the keyword arguments of Case.gas)."""
    sections = read_ini(case.dir / 'input.ini')
    sections['ICE-Parameters'].update(parameters)
    write_ini(case.dir / 'input.ini', sections)
    if gas is not None:
        case.gas(**gas)


def columns(closure):
    """Indices of T, f, chi, rho, u, n in the output of a solidifying family of this closure."""
    nb = len(LAYOUT[closure])
    return dict(T=nb - 2, n=nb - 1, f=nb, chi=nb + 1, rho=0, u=1)


def check_state(rep, label, sol, closure, T_ref, f_ref, regime, gap):
    c = columns(closure)
    temp, frac, chi = sol['var'][c['T']], sol['var'][c['f']], sol['var'][c['chi']]
    t = sol['time']
    spread = max(max(v) - min(v) for v in (temp, frac, chi))
    dT, df = abs(temp[0] - T_ref), abs(frac[0] - f_ref)
    if regime == 'plateau':
        tolT, tolf = 1.0e-9, 51.12 * DT + 1.0e-6
    elif regime == 'solid-after-jump':
        tolT, tolf = abs(T_ref - gap[1]) / TAU_S * 0.2695 * DT * 1.2 + 0.018, 0.0
    elif regime == 'hypercooled':
        tolT, tolf = 0.35, 0.0
    else:
        tolT, tolf = 1.0e-5 * gap[0], 1.0e-5
    chi_ref = 1.0 if f_ref > 0.0 else 0.0
    print('       t = %.6e s: T %.4f K (oracle %.4f), f %.6f (%.6f), chi %.1f; off %.2e K, %.2e'
          % (t, temp[0], T_ref, frac[0], f_ref, chi[0], dT, df))
    rep.check(dT <= tolT and df <= tolf and chi[0] == chi_ref and spread <= 1.0e-10 * max(1.0, abs(temp[0])),
              '%s: %s, T within %.1e K, f within %.1e, chi = %g, 8 cells equal' % (label, regime, tolT, tolf, chi_ref))
    rho, u, n = sol['var'][c['rho']], sol['var'][c['u']], sol['var'][c['n']]
    rep.check(max(abs(x - UG) for x in u) <= 1.0e-12 * UG and max(abs(x - RHO_P) for x in rho) <= 1.0e-12 * RHO_P
              and max(abs(x - NDENS) for x in n) <= 1.0e-9 * NDENS,
              '%s: velocity, density and number density unchanged' % label)


def n1(rep):
    print('   N1  cooling from 2400 K in a 600 K gas (tau_l %.6e s, tau_s %.6e s)' % (TAU_L, TAU_S))
    tn = TAU_L * math.log((2400.0 - 600.0) / (TN - 600.0))
    out = {}
    for t_end, regime in ((0.5 * tn, 'liquid'), (4.381910e-3, 'plateau'), (4.983360e-3, 'plateau'),
                          (5.945681e-3, 'plateau'), (0.010, 'solid-after-jump'), (0.015, 'solid-after-jump')):
        case = closed_cell('N1-%s' % ('%.6e' % t_end), 2400.0, 600.0, t_end)
        sol = attempt(rep, 'N1 t = %.6e' % t_end, case)
        if sol is None:
            continue
        T_ref, f_ref = cooling(sol['time'], 2400.0, 600.0)
        check_state(rep, 'N1 %s' % regime, sol, 'MK', T_ref, f_ref, regime, (1800.0, 600.0))
        out[t_end] = (case, sol)
    return out


def n2(rep):
    h2 = 4.0e5
    print('   N2  h-fus = 4e5 (f0 = %.4f): the hypercooled jump' % (CL * (TM - TN) / h2))
    tn = TAU_L * math.log((2400.0 - 600.0) / (TN - 600.0))
    case = closed_cell('N2', 2400.0, 600.0, tn + 50 * DT, phase=('A 1 solidification=on h-fus=4e5 cp-solid=600',))
    sol = attempt(rep, 'N2', case)
    if sol is not None:
        T_ref, f_ref = cooling(sol['time'], 2400.0, 600.0, h=h2)
        check_state(rep, 'N2', sol, 'MK', T_ref, f_ref, 'hypercooled', (1800.0, 600.0))


def n3(rep, closure='MK'):
    print('   N3  %s: a solid at 1500 K in a 3000 K gas melts at T_m' % closure)
    case = closed_cell('N3-%s' % closure, 1500.0, 3000.0, 3.0 * TAU_S, closure=closure)
    sol = attempt(rep, 'N3 %s' % closure, case)
    if sol is not None:
        T_ref, f_ref = heating(sol['time'], 1500.0, 3000.0)
        check_state(rep, 'N3 %s melting' % closure, sol, closure, T_ref, f_ref, 'melting', (1500.0, 3000.0))


def n4(rep):
    print('   N4  a mush at T_m, f = 0.8, chi = 1 in a 3000 K gas')
    rate = A * (3000.0 - TM) / (M * HF)
    for t_end in (5.0e-3, 0.8 / rate + 2.0e-3):
        case = closed_cell('N4-%.6e' % t_end, TM, 3000.0, t_end, full=(0.8, 1.0))
        sol = attempt(rep, 'N4 t = %.6e' % t_end, case)
        if sol is not None:
            T_ref, f_ref = heating(sol['time'], TM, 3000.0, f_start=0.8)
            check_state(rep, 'N4 %s' % ('melting' if f_ref > 0 else 'melted'), sol, 'MK', T_ref, f_ref,
                        'melting', (673.0, 3000.0))


def n5(rep):
    for closure in ('IG', 'AG'):
        print('   N5  %s closure' % closure)
        for t_end, regime in ((4.983360e-3, 'plateau'), (0.015, 'solid-after-jump')):
            case = closed_cell('N5-%s-%.6e' % (closure, t_end), 2400.0, 600.0, t_end, closure=closure)
            sol = attempt(rep, 'N5 %s t = %.6e' % (closure, t_end), case)
            if sol is not None:
                T_ref, f_ref = cooling(sol['time'], 2400.0, 600.0)
                check_state(rep, 'N5 %s %s' % (closure, regime), sol, closure, T_ref, f_ref, regime, (1800.0, 600.0))
        n3(rep, closure)


def n6(rep, straight):
    print('   N6  restart in place at the plateau midpoint, on to 0.015 s')
    case = closed_cell('N6', 2400.0, 600.0, 4.983360e-3)
    if attempt(rep, 'N6 first part', case) is None:
        return
    # the same case directory: the restart reads OUTPUT/part-field.tec (a new Case would wipe it)
    restart_in_place(case, {'newrun': False, 'time-threshold': 0.015})
    try:
        run_ice(case.dir)
    except RuntimeError as exc:
        rep.check(False, 'N6 restart: %s' % exc)
        return
    sol = read_solution(case.dir / 'OUTPUT/part-field.tec')
    ref = straight.get(0.015)
    if ref is None:
        rep.check(False, 'N6: the straight N1 run to 0.015 s is missing')
        return
    ref_sol = ref[1]
    worst = max(abs(a - b) / max(abs(b), 1.0e-300) for v in range(8) for a, b in zip(sol['var'][v], ref_sol['var'][v]))
    print('       restarted T %.10f K, straight %.10f K; worst relative difference %.2e'
          % (sol['var'][4][0], ref_sol['var'][4][0], worst))
    rep.check(worst <= 1.0e-12, 'N6: the restart reaches the straight run to %.1e relative' % worst)


def n7(rep):
    print('   N7  two materials: A solidifies, B (cp 1000) does not')
    zones = [dict(cp=CL, rho=RHO_M, zone='A'), dict(cp=1000.0, rho=RHO_M, zone='B')]
    case = closed_cell('N7', 2400.0, 600.0, 0.015, zones=zones, extra=('B 1',))
    sol = attempt(rep, 'N7', case)
    solo = closed_cell('N7-B-alone', 2400.0, 600.0, 0.015, phase=('B 1',), zones=[dict(cp=1000.0, rho=RHO_M, zone='B')])
    sol_b = attempt(rep, 'N7 B alone', solo)
    if sol is None or sol_b is None:
        return
    T_a = sol['var'][4][0]
    T_ref_a, _ = cooling(sol['time'], 2400.0, 600.0)
    tau_b = M * 1000.0 / A
    T_ref_b = 600.0 + 1800.0 * math.exp(-sol['time'] / tau_b)
    T_b = sol['var'][8 + 4][0]
    print('       A %.4f K (oracle %.4f), B %.4f K (oracle %.4f)' % (T_a, T_ref_a, T_b, T_ref_b))
    rep.check(abs(T_a - T_ref_a) <= abs(T_ref_a - 600.0) / TAU_S * 0.2695 * DT * 1.2 + 0.018,
              'N7: family A follows its solidification history')
    rep.check(abs(T_b - T_ref_b) <= 1.0e-5 * 1800.0, 'N7: family B relaxes as a liquid with cp 1000')
    same = all(a == b for v in range(6) for a, b in zip(sol['var'][8 + v], sol_b['var'][v]))
    rep.check(same, 'N7: family B\'s columns equal those of a run with B alone, value for value')


def selected(leg):
    """ICE_SOLID_LEGS=N1,N3 runs those legs only (the mutation probes use it); every leg by default."""
    legs = os.environ.get('ICE_SOLID_LEGS')
    return legs is None or leg in legs.split(',')


def main():
    rep = Report('N. Solidification, closed cell')
    straight = n1(rep) if selected('N1') or selected('N6') else {}
    for leg, fn in (('N2', n2), ('N3', n3), ('N4', n4), ('N5', n5), ('N7', n7)):
        if selected(leg):
            fn(rep)
    if selected('N6'):
        n6(rep, straight)
    rep.close(WORK)


if __name__ == '__main__':
    main()
