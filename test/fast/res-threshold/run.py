"""Fast test: res-threshold is relative to the first non-zero density residual of a steady run.

  scaling   a 1-D stream that carries a density wave out through its outlet, steady, MUSCL,
            res-threshold = 1e-5, at the loading 1 and at 2^-16 times it, an exact scaling of the MK
            equations at fixed velocity (every density stays above rho_empty). Both runs must end at the
            same iteration, before iter-threshold, on the first residual below 1e-5 times their first
            one. An absolute test ends the scaled run at its first iteration, whose residual is already
            below the threshold, and the full-loading run tens of iterations later.
  at rest   a cloud at rest with a density wave in a uniform gas stream (Stokes drag), one-stage Euler
            steps: the first iteration moves no mass (every velocity is zero), so its density residual
            is exactly zero. The run must not end there: the reference is the first non-zero residual.
  frozen    a uniform cloud moving with its outlet: every residual is exactly zero, nothing ends the
            run but iter-threshold, which must run all its iterations.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import Case, Physics, GREEN, RED, RESET                    # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
THR, ITERS = 1.0e-5, 20000
VINJ, TINJ, RINJ, RHO_AL = 10.0, 300.0, 1.0e-5, 1000.0
MASS = RHO_AL * 4.0 / 3.0 * math.pi * RINJ**3


def set_keys(case, **keys):
    """Overwrite keys of the case's input.ini (written by Case.ini): res-threshold, res-diter."""
    path = case.dir / 'input.ini'
    lines = path.read_text().split('\n')
    for i, line in enumerate(lines):
        key = line.split('=')[0].strip()
        if key in keys:
            lines[i] = '%s = %s' % (key, keys[key])
    path.write_text('\n'.join(lines))


def history(case):
    """[(iteration, density residual)] of the run, one row per iteration (res-diter = 1)."""
    rows = []
    for line in (case.dir / 'OUTPUT' / 'part-residual-history.dat').read_text().split('\n'):
        f = line.split()
        if len(f) >= 3:
            rows.append((int(f[0]), float(f[2])))
    return rows


def report(name, ok, detail):
    print('  %-9s %s  %s' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET), detail))
    return 0 if ok else 1


def ended_on_residual(rows):
    """(ok, text): the last row is the first, from the first non-zero residual on, whose residual is <= THR x it."""
    k0 = next((k for k, (_, r) in enumerate(rows) if r > 0.0), None)
    if k0 is None:
        return False, 'no non-zero residual'
    ref = rows[k0][1]
    below = [k for k, (_, r) in enumerate(rows) if k > k0 and r <= THR * ref]
    ok = bool(below) and below[0] == len(rows) - 1
    return ok, 'reference %.3e, last %.3e (ratio %.2e)' % (ref, rows[-1][1], rows[-1][1] / ref)


def stream(work, scale):
    case = Case(work, nx=32, Lx=0.32)
    rho = lambda x: scale * (1.0 + 0.5 * math.sin(2.0 * math.pi * x / 0.32))         # noqa: E731
    case.particles(rho, VINJ, 0.0, 0.0, TINJ, lambda x: rho(x) / MASS)
    case.boundaries(inlets=[{'gp': scale * VINJ, 'v': VINJ, 'T': TINJ, 'r': RINJ}])
    case.ini(iters=ITERS, cfl=0.8, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=RHO_AL,
             numerics={'time-accurate': False})
    set_keys(case, **{'res-threshold': THR, 'res-diter': 1})
    case.run()
    return history(case)


def scaling():
    full = stream(WORK / 'scaling-1', 1.0)
    small = stream(WORK / 'scaling-2^-16', 2.0**-16)
    ok_f, txt_f = ended_on_residual(full)
    ok_s, txt_s = ended_on_residual(small)
    same = full[-1][0] == small[-1][0]
    ok = ok_f and ok_s and same and 1 < full[-1][0] < ITERS
    return report('scaling', ok, 'iterations %d and %d; loading 1: %s; 2^-16: %s'
                  % (full[-1][0], small[-1][0], txt_f, txt_s))


def at_rest():
    ph = Physics()
    case = Case(WORK / 'at-rest', nx=32, Lx=0.32)
    rho = lambda x: 1.0 + 0.5 * math.sin(2.0 * math.pi * x / 0.32)                  # noqa: E731
    case.particles(rho, 0.0, 0.0, 0.0, ph.Tg, lambda x: rho(x) / MASS)
    case.gas(rho=ph.rho_g, u=ph.ug, v=ph.vg, w=0.0, T=ph.Tg, R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries()
    case.ini(iters=ITERS, cfl=0.8, rk='euler', drag='Stokes', heat='Stokes', rho_al=RHO_AL,
             numerics={'time-accurate': False})
    set_keys(case, **{'res-threshold': THR, 'res-diter': 1})
    case.run()
    rows = history(case)
    ok_r, txt = ended_on_residual(rows)
    ok = rows[0][1] == 0.0 and len(rows) > 1 and ok_r and rows[-1][0] < ITERS
    return report('at rest', ok, 'first residual %g, ended at iteration %d; %s' % (rows[0][1], rows[-1][0], txt))


def frozen():
    case = Case(WORK / 'frozen', nx=16, Lx=0.16)
    case.particles(1.0, VINJ, 0.0, 0.0, TINJ, 1.0 / MASS)
    case.boundaries()
    case.ini(iters=40, cfl=0.8, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=RHO_AL,
             numerics={'time-accurate': False})
    set_keys(case, **{'res-threshold': THR, 'res-diter': 1})
    case.run()
    rows = history(case)
    ok = len(rows) == 40 and rows[-1][0] == 40 and all(r == 0.0 for _, r in rows)
    return report('frozen', ok, '%d iterations of 40, largest residual %g' % (rows[-1][0], max(r for _, r in rows)))


if __name__ == '__main__':
    print('res-threshold relative to the first non-zero residual')
    rc = scaling() + at_rest() + frozen()
    sys.exit(1 if rc else 0)
