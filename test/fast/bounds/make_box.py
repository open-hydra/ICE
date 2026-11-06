"""Generated inputs of the bounds-checked leg (run.sh), built with the drivers of the verification cases that own them.

  python3 make_box.py solid <dir>   N4's mush (N-solid-cooling): T_m, f = 0.8, chi = 1 in a 3000 K gas, so the solid
                                    reconstruction runs from the first step
  python3 make_box.py drop <dir>    the M2 drop (M-drop-heating: water, d 50 um, 600 K gas, evaporation = CEM,
                                    boiling = ZGR) with liquid-conduction = P2T: cells 1-4 above the boiling
                                    temperature (600 K: the ZGR branch), cells 5-8 below it (the P2T surface balance),
                                    a Ts_p1 probe in each half, written every iteration
  python3 make_box.py check <dir>   after the drop's run: each probe holds five records (five values each, however the
                                    compiler wraps them), the hot cell's surface sits at the boiling temperature
                                    (Ts = T_B to 1e-9, T > T_B: ZGR), the cool cell's below its mean temperature
                                    (Ts < T - 1 K: the P2T profile; Ts = T without it), and the banner reports both models

Both runs stop after five iterations. <dir> is created anew; its name must be lower case (closed_cell lowers it).
"""
import importlib.util
import sys
from pathlib import Path

VERIFICATION = Path(__file__).resolve().parents[2] / 'verification'
TOKENS = 'liquid-conduction=P2T k-liq=0.6 mu-liq=1e-3'


def driver(name):
    """The run.py of a verification case, as a module (its main() is not run)."""
    spec = importlib.util.spec_from_file_location(name.replace('-', '_'), str(VERIFICATION / name / 'run.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def solid(target):
    N = driver('N-solid-cooling')
    N.WORK = target.parent
    return N.closed_cell(target.name, N.TM, 3000.0, 5 * N.DT, full=(0.8, 1.0))


def drop_leg(M, name):
    return M.Leg(name, 958.0, 4186.0, 50e-6, 640.002, 600.0, 'CEM', 0.0, 551.7783324, p2t=0.6, zgr=True,
                 h=1e-6, steps=5, every=1, physics={'evaporation': 'CEM', 'boiling': 'ZGR'}, tokens=TOKENS)


def drop(target, M):
    M.WORK = target.parent
    leg = drop_leg(M, target.name)
    # Leg.run's case, without the run, and with T above T_B (600 K) in cells 1-4, below it in cells 5-8
    case = M.Case(M.WORK / leg.name, nx=8, Lx=1.0)
    case.particles(rho=1.0, u=0.0, v=0.0, w=0.0, T=lambda x: leg.T0 if x < 0.5 else 560.0, n=leg.n)
    case.gas(rho=M.RHO_G, u=0.0, v=0.0, w=0.0, T=leg.Tg, R=M.R_G, gam=M.GAM, k=M.K_G, mu=M.MU_G)
    case.boundaries('extrapolation')
    phys = {'latent-heat': M.LV, 'vapour-molar-mass': M.MV, 'vapour-specific-heat': M.CPV,
            'boiling-temperature': 1.0 / leg.invtb, 'lewis-number': M.LE, 'vapour-mass-fraction': leg.yinf}
    phys.update(leg.physics)
    case.ini(t_end=leg.steps * leg.h, rk='RK2', drag='Stokes', heat='Stokes', rho_al=leg.rho_l, cs=leg.c,
             dt_max=leg.h, physics=phys)
    case.phase('A 1 ' + leg.tokens)
    M.ini_set(case.dir, 'ICE-Probes', probe1='hot', probe2='cool')
    M.ini_set(case.dir, 'hot', variables='T_p1 rho_p1 n_p1 Ts_p1', diter=1, index_position='1 4 1 1')
    M.ini_set(case.dir, 'cool', variables='T_p1 rho_p1 n_p1 Ts_p1', diter=1, index_position='1 5 1 1')
    return case


def records(path):
    """The probe's records of five values (time, T, rho, n, Ts), which a compiler may wrap over two lines."""
    vals = [float(x) for x in Path(path).read_text().split()]
    return [vals[i:i + 5] for i in range(0, len(vals), 5)], len(vals)


def check(target, M):
    """The witnesses of the drop's run; returns the failures."""
    TB = drop_leg(M, target.name).TB
    bad = []
    for name in ('hot', 'cool'):
        rows, nvals = records(target / 'OUTPUT' / ('%s.txt' % name))
        if nvals != 25:
            bad.append('%s.txt holds %d values, not five records of five' % (name, nvals))
            continue
        for t, T, rho, n, ts in rows:
            if name == 'hot' and not (abs(ts / TB - 1.0) <= 1e-9 and T > TB):
                bad.append('hot cell at %g: Ts %.12g, T %.12g, not Ts = T_B %.12g < T (ZGR)' % (t, ts, T, TB))
            if name == 'cool' and not (ts < T - 1.0 and T < TB):
                bad.append('cool cell at %g: Ts %.12g, T %.12g, not Ts < T - 1 K < T_B (P2T)' % (t, ts, T))
        print('[make_box] %s: Ts %.12g, T %.12g at the fifth record (T_B %.12g)' % (name, rows[-1][4], rows[-1][1], TB))
    log = (target / 'log').read_text(errors='replace')
    for text in ('- Liquid  --> P2T', '- Boiling --> ZGR'):
        if text not in log:
            bad.append('the banner does not report "%s"' % text)
    return bad


if __name__ == '__main__':
    kind, target = sys.argv[1], Path(sys.argv[2]).resolve()
    M = driver('M-drop-heating')
    if kind == 'check':
        failures = check(target, M)
        for f in failures:
            print('[make_box] %s' % f)
        sys.exit(1 if failures else 0)
    if target.name != target.name.lower():
        raise SystemExit('make_box: the directory name must be lower case: %s' % target.name)
    if kind == 'solid':
        case = solid(target)
    elif kind == 'drop':
        case = drop(target, M)
    else:
        raise SystemExit('make_box: solid, drop or check, not %s' % kind)
    if case.dir.resolve() != target:
        raise SystemExit('make_box: the case went to %s, not %s' % (case.dir, target))
    M.ini_set(case.dir, 'ICE-Parameters', iter_threshold=5)
