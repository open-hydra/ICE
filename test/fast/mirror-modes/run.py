"""Fast test: the components odd across mirror faces one cell apart are not integrated.

Between two faces that mirror the flow one cell apart -- the k faces of a 2-D mesh,
the j and k faces of a 1-D one -- the components of the state odd across them have
the mode dw/dt = -2 a w / dz (a the signal speed of the closure, dz the width of the
cell across the faces). The time step follows the directions the mesh resolves only,
so on a layer thin against its other cells that mode is far outside the stability
interval of the scheme, and round-off in w grows until the state is unphysical.
ICE sets those components to zero at set-up and does not integrate them
(setup_freeze). Pinned here, each on a cloud at rest (rho = 1, P = 1) whose odd
components start at 1e-8 of the signal speed (the AG tensor components at 1e-8 P):

  slab        a plane layer 1e-5 m thick under cells 2.5 mm wide, k faces 300, IG
              and AG, RK2 at cfl 0.8 with dt-max out of the way, 200 steps: w (AG:
              also P13, P23) zero in every cell, the cloud at rest, the warning.
  tube        a 1-D tube 1e-5 m across, j and k faces 300, IG: v and w zero.
  wedge       a wedge of one degree whose side faces carry 300 instead of 200:
              the setup stops and says why.
"""
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import ICE_BIN, Report, linspace, read_solution, write_ini, write_tec  # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
H = 1e-5                 # thickness across the mirror faces [m]
NX, NY, L = 16, 8, 0.04  # cells along x and y, length along x [m]: cells of 2.5 mm
EPS, STEPS = 1e-8, 200


def bc(path, nx, ny, codes):
    """ATLAS bc.txt of one block with one cell in k: faces 1 to 6 carry codes[f]."""
    rows = []
    for f, ii, jj in ((1, [1], range(1, ny + 1)), (2, [nx], range(1, ny + 1)),
                      (3, range(1, nx + 1), [1]), (4, range(1, nx + 1), [ny]),
                      (5, range(1, nx + 1), range(1, ny + 1)), (6, range(1, nx + 1), range(1, ny + 1))):
        rows += ['%8d%8d%8d%8d%8d%8d' % (1, i, j, 1, f, codes[f]) for j in jj for i in ii]
    Path(path).write_text('\n'.join(rows) + '\n')


def build(name, closure, xn, yn, zn, nodes, codes, odd):
    """A cloud at rest on the mesh of node lists (xn, yn, zn), or of explicit node arrays (nodes), perturbed in the
    slots of odd: 'w' (AG also P13, P23) and 'v' (AG also P12, P23)."""
    work = WORK / name
    if work.exists():
        shutil.rmtree(str(work))
    (work / 'INPUT').mkdir(parents=True)
    (work / 'OUTPUT').mkdir()
    nx, ny = len(xn) - 1, len(yn) - 1
    xc = [0.5 * (xn[i] + xn[i + 1]) for i in range(nx)]
    yc = [0.5 * (yn[j] + yn[j + 1]) for j in range(ny)]
    s = [EPS * math.sin(2 * math.pi * x / xn[-1]) * math.cos(math.pi * (y - yn[0]) / (yn[-1] - yn[0]))
         for y in yc for x in xc]
    c = math.sqrt(3.0)
    zero, one = [0.0] * (nx * ny), [1.0] * (nx * ny)
    w = [c * a for a in s] if 'w' in odd else zero
    v = [c * a for a in s] if 'v' in odd else zero
    n = [1.0 / (4.0 / 3.0 * math.pi * 1e-15 * 2000.0)] * (nx * ny)
    if closure == 'IG':
        names, cells = ['rho', 'u', 'v', 'w', 'P', 'T', 'n'], [one, zero, v, w, one, [300.0] * (nx * ny), n]
    else:
        p12 = s if 'v' in odd else zero
        p13 = s if 'w' in odd else zero
        names = ['rho', 'u', 'v', 'w', 'P11', 'P12', 'P13', 'P22', 'P23', 'P33', 'T', 'n']
        p23 = s if odd else zero
        cells = [one, zero, v, w, one, p12, p13, one, p23, one, [300.0] * (nx * ny), n]
    if nodes is None:
        write_tec(work / 'INPUT/part-ic.tec', names, xn, yn, zn, cells)
    else:
        write_nodes(work / 'INPUT/part-ic.tec', names, nodes, cells)
    bc(work / 'INPUT/part-bc.txt', nx, ny, codes)
    (work / 'INPUT/part-phase.txt').write_text('condensed-dispersed phase\nA 1\n')
    write_ini(work / 'input.ini', {
        'ICE-Parameters': {'iter-threshold': STEPS, 'res-threshold': 0},
        'ICE-Numerics': {'time-scheme': 'RK2', 'cfl': 0.8, 'dt-max': 1e3, 'time-accurate': False,
                         'space-reconstruction': 'MUSCL', 'flux-limiter': 'vanleer'},
        'ICE-Physics': {'density': 2000.0, 'specific-heat': 1000.0},
        'ICE-Family1': {'closure': closure},
        'ICE-IO': {'ic-format': 'tecplot ascii', 'sol-format': 'tecplot ascii', 'shell-diter': STEPS}})
    return work


def write_nodes(path, names, nodes, cells):
    """write_tec for explicit node arrays nodes[c][k][j][i] (the wedge)."""
    X = nodes
    ni, nj = len(X[0][0][0]), len(X[0][0])
    lines = [' VARIABLES ="x" "y" "z" ' + ' '.join('"%s"' % v for v in names),
             ' ZONE  T = B1, I=%d, J=%d, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED)'
             % (ni, nj, 3 + len(names))]
    for comp in X:
        lines += [' %.15E' % comp[k][j][i] for k in range(2) for j in range(nj) for i in range(ni)]
    for cell in cells:
        lines += [' %.15E' % a for a in cell]
    Path(path).write_text('\n'.join(lines) + '\n')


def run(work):
    """Exit code and log of ICE in work."""
    if not ICE_BIN.exists():
        raise SystemExit('[fast] no ICE binary at %s -- build first' % ICE_BIN)
    with open(str(work / 'log'), 'w') as out:
        rc = subprocess.call([str(ICE_BIN)], cwd=str(work), stdout=out, stderr=subprocess.STDOUT,
                             env=dict(os.environ, OMP_NUM_THREADS='1', KMP_STACKSIZE='100M'))
    return rc, (work / 'log').read_text(errors='replace')


def largest(sol, names):
    return max(abs(a) for nm in names for a in sol['var'][sol['names'].index(nm + '_p1')])


def main():
    r = Report('Mirror modes')
    xn, yn = linspace(0.0, L, NX + 1), linspace(0.0, L / 2, NY + 1)
    for closure, slots in (('IG', ['w']), ('AG', ['w', 'P13', 'P23'])):
        work = build('slab-' + closure, closure, xn, yn, [-H / 2, H / 2], None,
                     {1: 300, 2: 300, 3: 300, 4: 300, 5: 300, 6: 300}, 'w')
        rc, log = run(work)
        if not r.check(rc == 0, 'slab %s: the run ends normally (exit %d)' % (closure, rc)):
            continue
        sol = read_solution(work / 'OUTPUT/part-field.tec')
        r.check(largest(sol, slots) == 0.0, 'slab %s: %s zero in every cell' % (closure, ', '.join(slots)))
        r.check(largest(sol, ['u', 'v']) < 1e-10, 'slab %s: the cloud stays at rest (|u|, |v| %.1e m/s)'
                % (closure, largest(sol, ['u', 'v'])))
        r.check('2-D plane mesh with mirror k faces: w' in log, 'slab %s: the warning names w' % closure)

    work = build('tube', 'IG', linspace(0.0, L, 2 * NX + 1), [-H / 2, H / 2], [-H / 2, H / 2], None,
                 {1: 300, 2: 300, 3: 300, 4: 300, 5: 300, 6: 300}, 'vw')
    rc, log = run(work)
    if r.check(rc == 0, 'tube: the run ends normally (exit %d)' % rc):
        sol = read_solution(work / 'OUTPUT/part-field.tec')
        r.check(largest(sol, ['v', 'w']) == 0.0, 'tube: v and w zero in every cell')
        r.check('1-D mesh with mirror j and k faces: v and w' in log, 'tube: the warning names v and w')

    th = [math.radians(-0.5), math.radians(0.5)]
    rn = linspace(1e-8, L / 2, NY + 1)
    nodes = [[[[x for x in xn] for _ in rn] for _ in th],
             [[[rr * math.cos(t) for _ in xn] for rr in rn] for t in th],
             [[[rr * math.sin(t) for _ in xn] for rr in rn] for t in th]]
    work = build('wedge-300', 'IG', xn, rn, None, nodes, {1: 300, 2: 300, 3: 300, 4: 300, 5: 300, 6: 300}, '')
    rc, log = run(work)
    r.check(rc != 0 and 'the k faces of an axisymmetric wedge must carry 200' in log,
            'wedge with side faces 300: the setup stops and says why (exit %d)' % rc)

    r.close(WORK)


if __name__ == '__main__':
    main()
