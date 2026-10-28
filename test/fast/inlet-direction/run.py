"""Fast test: an inlet injects into the domain on a face whose normal points out of it.

The face normal N of an i-face points to growing i: into the domain on face 1, out of it on face 2.
A 1-D duct, MK, the inlet on face 2 (x = L), particles at rest elsewhere, 200 time-accurate steps:

  normal   402, direction `normal`: the particles must enter along -x at the inlet speed, the cell at
           the inlet carrying the injected density gp / v (with the floor of a cloud at rest elsewhere).
  angles   401 (speed kV times the gas velocity into the domain) with explicit angles alpha = 180 deg,
           beta = 0 (the direction -x), in a gas moving along -x: the particles must enter along -x.
           The speed is a ratio of the gas velocity INTO the domain, positive for a gas that enters.
  mixed    `normal` for alpha and a number for beta: refused at set-up, with the record named.
"""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import Case, Physics, GREEN, RED, RESET, ICE_BIN, inlet_payload      # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
NX, LX = 16, 0.16
GP, VINJ, TINJ, RINJ, RHO_AL = 1.0, 10.0, 300.0, 1.0e-5, 1000.0
MASS = RHO_AL * 4.0 / 3.0 * math.pi * RINJ**3
RHO0 = 1.0e-3                      # the cloud at rest everywhere at the start (above rho_empty)
RHO, U = 0, 1


def write_bc(case, payload, code):
    """The duct's table: 400 everywhere, `code` with `payload` on face 2, 0 on faces 5 and 6."""
    h = lambda i, j, k, f, c: '%8d%8d%8d%8d%8d%8d' % (1, i, j, k, f, c)          # noqa: E731
    rows = [h(1, 1, 1, 1, 400), h(NX, 1, 1, 2, code), '   ' + payload]
    rows += [h(i, 1, 1, 3, 400) for i in range(1, NX + 1)] + [h(i, 1, 1, 4, 400) for i in range(1, NX + 1)]
    rows += [h(i, 1, 1, 5, 0) for i in range(1, NX + 1)] + [h(i, 1, 1, 6, 0) for i in range(1, NX + 1)]
    (case.dir / 'INPUT' / 'part-bc.txt').write_text('\n'.join(rows) + '\n')


def duct(name, payload, code, gas_u=None):
    case = Case(WORK / name, nx=NX, Lx=LX)
    case.particles(RHO0, 0.0, 0.0, 0.0, TINJ, RHO0 / MASS)
    if gas_u is not None:
        ph = Physics(ug=gas_u)
        case.gas(rho=ph.rho_g, u=ph.ug, v=0.0, w=0.0, T=ph.Tg, R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    write_bc(case, payload, code)
    case.ini(iters=200, cfl=0.5, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=RHO_AL)
    return case


def report(name, ok, detail):
    print('  %-7s %s  %s' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET), detail))
    return 0 if ok else 1


def normal():
    sol = duct('normal', inlet_payload(GP, VINJ, TINJ, RINJ), 402).run()
    rho, u = sol['var'][RHO], sol['var'][U]
    ok = abs(u[-1] + VINJ) <= 1e-9 * VINJ and abs(rho[-1] - GP / VINJ) <= 1e-9 * GP / VINJ
    return report('normal', ok, 'inlet cell: u %+.6f m/s (expected %+.1f), rho %.6f (expected %.6f)'
                  % (u[-1], -VINJ, rho[-1], GP / VINJ))


def angles():
    krho, kv, kt = 0.5, 1.0, 1.0
    payload = '%.10e %.10e %.10e %.10e %.10e %.10e 0.0 Dirac 0.0' % (krho, kv, math.pi, 0.0, kt, RINJ)
    sol = duct('angles', payload, 401, gas_u=-VINJ).run()
    u = sol['var'][U]
    ok = abs(u[-1] + kv * VINJ) <= 1e-9 * VINJ
    return report('angles', ok, 'inlet cell: u %+.6f m/s (expected %+.1f: kV times the gas entering)'
                  % (u[-1], -kv * VINJ))


def mixed():
    payload = '%.10e %.10e normal 0.0 %.10e %.10e 0.0 Dirac 0.0' % (GP, VINJ, TINJ, RINJ)
    case = duct('mixed', payload, 402)
    out = subprocess.run([str(ICE_BIN)], cwd=str(case.dir), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    text = out.stdout.decode('utf-8', 'replace')
    said = 'both normal or both numbers' in text
    return report('mixed', out.returncode != 0 and said,
                  'exit %d, %s' % (out.returncode, 'reported' if said else 'no message'))


if __name__ == '__main__':
    print('Inlet direction on a face whose normal points out of the domain')
    rc = normal() + angles() + mixed()
    sys.exit(1 if rc else 0)
