"""P. Direct combustion in a closed cell - Gallier and Godfroy's burning drop (JPP 25(2) 2009, section II.A).

Their aluminium (Tables 1 and 2): rho_p 2000 kg/m^3, c_p 1177 J/(kg K), Q 9.53e6 J/kg, L_v 10.8e6 J/kg, T_sat 2791 K,
D_inj 120 um, D_res 60 um; a gas of mu 9.1e-5 Pa s, Pr 0.4, gamma 1.16, cp 2000 J/(kg K) at 3450 K and 8.71 MPa,
whose one species is the vapour that burns (INPUT/phase.txt). Phase line "A 1 combustion=direct q-comb=9.53e6
d-res=6e-5"; latent-heat and boiling-temperature (the burning temperature) in [ICE-Physics]. Eight uniform cells move
with the gas (drag = NoDrag), heat-transfer = Stokes (Nu = 2 exactly), emissivity 0, so every cell is one drop at rest
in the gas: m = rho pi D^3/6, k = mu cp/Pr, B = (cp (T_g - T_sat) + Q)/L_v, and

  burning   T = T_sat, D^2 = D_0^2 - K t, K = 8 k ln(1 + B)/(cp rho_p), down to D_res at t_res; with no residue the
            drop would burn out at D_0^2/K = 22.8 ms (their eq. 10 gives 23 ms)
  residue   D = D_res, T = T_g - (T_g - T_sat) exp(-(t - t_res)/tau), tau = m_res c_p/(pi D_res k Nu)
  heating   a drop below T_sat heats inertly, T = T_g - (T_g - T_0) exp(-t/tau_0), until T_sat at t_ign, then burns

  P1  a drop at T_sat: D at a quarter, a half and three quarters of t_res, T = T_sat to 1e-9 K
  P2  the residue: D = D_res at t_res + tau and t_res + 3 tau, T on its relaxation
  P3  a drop at 2000 K: heats to T_sat, then burns: T and D at t_ign + t_res/2
  P4  P1's midpoint with the IG and AG closures (P = 1e-6)
  P5  gas-phase = gas reads INPUT/gas-phase.txt; the refusals: two gas species, no gas phase file, two materials,
      evaporation, no d-res, no boiling-temperature
  P6  the [ICE-Physics] keys combustion, q-comb and d-res: alone on a bare "A 1" line they burn as P1; a d-res token
      overrides the key; a combustion=none token switches it off; d-res = 0 is refused

Tolerances (time-scheme RK2, dt-max = 1e-5 binding). Burning has no event: D to 1e-6 relative. The residue switch
happens at the end of a stage, so D lies up to 1.2 K dt/(2 D_res) below D_res, and T is late by up to dt on its
relaxation; the ignition switch likewise leaves T up to dT/dt dt above T_sat and D^2 off by up to K dt.
"""
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import Case, Report                                               # noqa: E402

WORK = HERE / 'work'

# --- Gallier and Godfroy's aluminium and gas ----------------------------------------------
RHO_M, CP, Q, LV, TSAT = 2000.0, 1177.0, 9.53e6, 10.8e6, 2791.0
D0, DRES = 120.0e-6, 60.0e-6
MU, PR, GAM, CPG, TG, PG = 9.1e-5, 0.4, 1.16, 2000.0, 3450.0, 8.71e6
R_G = CPG * (GAM - 1.0) / GAM
KG = MU * CPG / PR
RHO_G = PG / (R_G * TG)
NU, UG = 2.0, 10.0
B = (CPG * (TG - TSAT) + Q) / LV
K = 8.0 * KG * math.log(1.0 + B) / (CPG * RHO_M)
T_RES = (D0 ** 2 - DRES ** 2) / K
TAU_RES = RHO_M * math.pi / 6.0 * DRES ** 3 * CP / (math.pi * DRES * KG * NU)
ALPHA = 1.0e-4
RHO_P = ALPHA * RHO_M
NDENS = RHO_P / (RHO_M * math.pi / 6.0 * D0 ** 3)
DT = 1.0e-5
PHASE = 'A 1 combustion=direct q-comb=9.53e6 d-res=6e-5'
PHYSICS = {'latent-heat': LV, 'boiling-temperature': TSAT}
LAYOUT = {'MK': 6, 'IG': 7, 'AG': 12}
PSEUDO = 1.0e-6


def diameter(rho, n):
    return (6.0 * rho / (math.pi * RHO_M * n)) ** (1.0 / 3.0)


def tau_of(d):
    return RHO_M * math.pi / 6.0 * d ** 3 * CP / (math.pi * d * KG * NU)


def closed_cell(label, T0, t_end, closure='MK', phase=(PHASE,), physics=None, gas_file='phase.txt',
                gas_lines=('GAS 30.14',), parameters=None, families=1):
    case = Case(WORK / label.lower(), nx=8, Lx=1.0)
    P = dict(IG=[PSEUDO], AG=[PSEUDO, 0.0, 0.0, PSEUDO, 0.0, PSEUDO]).get(closure, [])
    names, fields = [], []
    for p in range(1, families + 1):
        names += ['%s%d' % (s, p) for s in ['rp', 'up', 'vp', 'wp'] + ['P%d_' % i for i in range(len(P))] + ['Tp', 'np']]
        fields += [RHO_P, UG, 0.0, 0.0] + P + [T0, NDENS]
    case.write_ic(names, fields)
    case.gas(rho=RHO_G, u=UG, v=0.0, w=0.0, T=TG, R=R_G, gam=GAM, k=KG, mu=MU)
    if gas_file is not None:
        (case.dir / 'INPUT' / gas_file).write_text('\n'.join(('ideal-gas phase',) + tuple(gas_lines)) + '\n')
    case.boundaries('extrapolation')
    case.phase(*phase)
    phys = dict(PHYSICS)
    phys.update(physics or {})
    case.ini(t_end=t_end, rk='RK2', drag='NoDrag', heat='Stokes', rho_al=RHO_M, cs=CP, dt_max=DT,
             physics=phys, closures=(closure,) * families)
    if parameters:
        text = (case.dir / 'input.ini').read_text()
        extra = ''.join('%s = %s\n' % kv for kv in parameters.items())
        (case.dir / 'input.ini').write_text(text.replace('[ICE-Parameters]\n', '[ICE-Parameters]\n' + extra, 1))
    return case


def attempt(rep, label, case):
    try:
        return case.run()
    except RuntimeError as exc:
        log = (case.dir / 'log').read_text(errors='replace').splitlines() if (case.dir / 'log').exists() else []
        hits = [s.strip() for s in log if '[ERROR]' in s]
        rep.check(False, '%s: %s%s' % (label, str(exc).split(' in ')[0], (': ' + hits[0]) if hits else ''))
        return None


def state(sol, closure):
    nb = LAYOUT[closure]
    rho, u, T, n = (sol['var'][i] for i in (0, 1, nb - 2, nb - 1))
    spread = max(max(v) - min(v) for v in (rho, T))
    return rho, u, T, n, spread


def check_cell(rep, label, sol, closure, D_ref, T_ref, tolD, tolT):
    rho, u, T, n, spread = state(sol, closure)
    D = diameter(rho[0], n[0])
    print('       t = %.6e s: D %.6e m (oracle %.6e), T %.6f K (oracle %.6f); off %.2e relative, %.2e K'
          % (sol['time'], D, D_ref, T[0], T_ref, abs(D - D_ref) / D_ref, abs(T[0] - T_ref)))
    rep.check(abs(D - D_ref) <= tolD and abs(T[0] - T_ref) <= tolT and spread <= 1.0e-10 * max(1.0, T[0]),
              '%s: D within %.1e m, T within %.1e K, 8 cells equal' % (label, tolD, tolT))
    rep.check(max(abs(x - UG) for x in u) <= 1.0e-12 * UG and max(abs(x - NDENS) for x in n) <= 1.0e-12 * NDENS,
              '%s: velocity and number density unchanged' % label)
    return D, T[0]


def p1(rep, closure='MK', fractions=(0.25, 0.5, 0.75)):
    for frac in fractions:
        t_end = frac * T_RES
        case = closed_cell('P1-%s-%.2f' % (closure, frac), TSAT, t_end, closure=closure)
        sol = attempt(rep, 'P1 %s t = %.6e' % (closure, t_end), case)
        if sol is None:
            continue
        D_ref = math.sqrt(D0 ** 2 - K * sol['time'])
        check_cell(rep, 'P%s %s burning at %.2f t_res' % ('1' if closure == 'MK' else '4', closure, frac), sol,
                   closure, D_ref, TSAT, 1.0e-6 * D_ref, 1.0e-9)


def p2(rep):
    tolD = 1.2 * K * DT / (2.0 * DRES)
    for m in (1.0, 3.0):
        t_end = T_RES + m * TAU_RES
        case = closed_cell('P2-%g' % m, TSAT, t_end)
        sol = attempt(rep, 'P2 t = %.6e' % t_end, case)
        if sol is None:
            continue
        dt_rel = sol['time'] - T_RES
        T_ref = TG - (TG - TSAT) * math.exp(-dt_rel / TAU_RES)
        rate = (TG - T_ref) / TAU_RES
        rho, u, T, n, _ = state(sol, 'MK')
        D = diameter(rho[0], n[0])
        print('       t = %.6e s: D %.9e m (D_res %.9e), T %.6f K (oracle %.6f)' % (sol['time'], D, DRES, T[0], T_ref))
        rep.check(DRES - tolD <= D <= DRES and abs(T[0] - T_ref) <= 1.2 * rate * DT + 1.0e-6,
                  'P2 the residue at t_res + %g tau: D in [D_res - %.1e, D_res], T on its relaxation within %.1e K'
                  % (m, tolD, 1.2 * rate * DT + 1.0e-6))


def p3(rep):
    T0 = 2000.0
    tau0 = tau_of(D0)
    t_ign = tau0 * math.log((TG - T0) / (TG - TSAT))
    t_end = t_ign + 0.5 * T_RES
    print('   P3  a drop at %.0f K: tau_0 %.6e s, t_ign %.6e s' % (T0, tau0, t_ign))
    case = closed_cell('P3', T0, t_end)
    sol = attempt(rep, 'P3', case)
    if sol is None:
        return
    D_ref = math.sqrt(D0 ** 2 - K * (sol['time'] - t_ign))
    rate = (TG - TSAT) / tau0
    tolD = 1.2 * K * DT / (2.0 * D_ref) + 1.0e-6 * D_ref
    check_cell(rep, 'P3 ignited and burning', sol, 'MK', D_ref, TSAT, tolD, 1.2 * rate * DT)
    rho, u, T, n, _ = state(sol, 'MK')
    rep.check(T[0] >= TSAT, 'P3 the burning drop holds at or just above T_sat')


def refusal(rep, label, needle, leg='P5', **kw):
    case = closed_cell('%s-%s' % (leg, label), TSAT, 10 * DT, **kw)
    try:
        case.run()
        rep.check(False, '%s %s: the run was not refused' % (leg, label))
        return
    except RuntimeError:
        pass
    log = (case.dir / 'log').read_text(errors='replace')
    hit = needle in log
    line = next((s.strip() for s in log.splitlines() if '[ERROR]' in s), '')
    print('       %s: %s' % (label, line))
    rep.check(hit, '%s %s: refused, naming "%s"' % (leg, label, needle))


def p5(rep):
    case = closed_cell('P5-gas-phase', TSAT, 0.25 * T_RES, gas_file='gas-phase.txt', parameters={'gas-phase': 'gas'})
    sol = attempt(rep, 'P5 gas-phase = gas', case)
    if sol is not None:
        rho, u, T, n, _ = state(sol, 'MK')
        D_ref = math.sqrt(D0 ** 2 - K * sol['time'])
        rep.check(abs(diameter(rho[0], n[0]) - D_ref) <= 1.0e-6 * D_ref,
                  'P5 gas-phase = gas: INPUT/gas-phase.txt read, the drop burns as P1')
    refusal(rep, 'two-species', 'needs one gas species', gas_lines=('O2 32.0', 'N2 28.0'))
    refusal(rep, 'no-gas-file', 'which is missing', gas_file=None)
    refusal(rep, 'two-materials', 'needs one material', phase=(PHASE, 'B 1'), families=2,
            physics={'emissivity': '0 0', 'density': '%g %g' % (RHO_M, RHO_M), 'specific-heat': '%g %g' % (CP, CP),
                     'latent-heat': '%g %g' % (LV, LV), 'boiling-temperature': '%g %g' % (TSAT, TSAT)})
    refusal(rep, 'evaporation', 'exclusive with evaporation', phase=(PHASE + ' evaporation=d2-law',),
            physics={'vapour-molar-mass': 26.98, 'vapour-specific-heat': 0, 'lewis-number': 1,
                     'vapour-mass-fraction': 0})
    refusal(rep, 'no-d-res', 'requires d-res > 0', phase=('A 1 combustion=direct q-comb=9.53e6',))
    case = closed_cell('P5-no-tboil', TSAT, 10 * DT)
    ini = (case.dir / 'input.ini').read_text()
    (case.dir / 'input.ini').write_text('\n'.join(s for s in ini.splitlines() if not s.startswith('boiling-temperature')))
    try:
        case.run()
        rep.check(False, 'P5 no-boiling-temperature: the run was not refused')
    except RuntimeError:
        log = (case.dir / 'log').read_text(errors='replace')
        rep.check('requires boiling-temperature' in log,
                  'P5 no-boiling-temperature: refused, naming "requires boiling-temperature"')


INI_KEYS = {'combustion': 'direct', 'q-comb': Q, 'd-res': DRES}


def p6(rep):
    t_end = 0.5 * T_RES
    for label, phase, extra, burns in (('keys alone', 'A 1', {}, True),
                                       ('a d-res token over d-res = 1e-4', 'A 1 d-res=6e-5', {'d-res': 1.0e-4}, True),
                                       ('a combustion=none token', 'A 1 combustion=none', {}, False)):
        phys = dict(INI_KEYS)
        phys.update(extra)
        case = closed_cell('P6-%s' % label.replace(' ', '-').replace('=', ''), TSAT, t_end, phase=(phase,), physics=phys)
        sol = attempt(rep, 'P6 %s' % label, case)
        if sol is None:
            continue
        rho, u, T, n, _ = state(sol, 'MK')
        D = diameter(rho[0], n[0])
        if burns:
            D_ref = math.sqrt(D0 ** 2 - K * sol['time'])
            rep.check(abs(D - D_ref) <= 1.0e-6 * D_ref and abs(T[0] - TSAT) <= 1.0e-9,
                      'P6 %s: burns as P1 (D off %.1e relative)' % (label, abs(D - D_ref) / D_ref))
        else:
            rep.check(abs(D - D0) <= 1.0e-12 * D0 and T[0] > TSAT + 1.0,
                      'P6 %s: inert, D unchanged and the drop heating (T %.3f K)' % (label, T[0]))
    phys = dict(INI_KEYS)
    phys['d-res'] = 0.0
    refusal(rep, 'ini-d-res-0', 'd-res must be > 0', leg='P6', phase=('A 1',), physics=phys)


def selected(leg):
    legs = os.environ.get('ICE_COMB_LEGS')
    return legs is None or leg in legs.split(',')


def main():
    rep = Report('P. Direct combustion, closed cell')
    print('   B %.6f, K %.6e m^2/s, t_res %.6e s, burn-out without residue %.6e s, tau_res %.6e s'
          % (B, K, T_RES, D0 ** 2 / K, TAU_RES))
    if selected('P1'):
        print('   P1  a drop at T_sat')
        p1(rep)
    if selected('P2'):
        print('   P2  the residue')
        p2(rep)
    if selected('P3'):
        p3(rep)
    if selected('P4'):
        for closure in ('IG', 'AG'):
            print('   P4  %s closure' % closure)
            p1(rep, closure, fractions=(0.5,))
    if selected('P5'):
        print('   P5  the gas phase file and the refusals')
        p5(rep)
    if selected('P6'):
        print('   P6  the [ICE-Physics] keys')
        p6(rep)
    rep.close(WORK)


if __name__ == '__main__':
    main()
