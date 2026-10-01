"""Fast test: the multigrid path runs, and neither threads nor ranks change it.

Grid sequencing starts the solve on the coarsest level and prolongates onto the
next one each time a level exhausts its iteration budget, so a run with
MG-levels > 1 touches code no single-level case reaches: the restriction and
prolongation operators, a boundary table per level, and one halo schedule per
level. Three things are pinned here.

  runs        MG-levels 1, 2 and 3 all reach the fine level and produce a finite
              field, and 2 and 3 differ from 1 by a bounded, non-zero amount.
              Zero would mean the coarse levels never reached the fine solution --
              which is what happened while fine2coarse_prim and coarse2fine_prim
              declared their prim dummies with one ghost layer instead of gc, so
              every index in them named a different cell than it reached.
  threads     1 vs ICE_MG_THREADS threads, bit-identical, at MG-levels 3. The
              restriction, prolongation and coarse-grid loops each sit inside
              their own parallel region; before that they carried a bare !$omp do
              with no region to bind to and ran on one thread.
  ranks       1 vs ICE_MG_RANKS ranks, bit-identical, at MG-levels 3. Needs an
              MPI build; skipped when ICE_MG_RANKS is 1.

The case is two blocks closed onto each other in x, so every i-face is a
connection and at two ranks every one of them crosses a rank boundary.
"""
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import ICE_BIN, Report, _fmt, _nodal, linspace, read_solution, write_ini  # noqa: E402

NB, NXB, NY = 2, 32, 8          # divisible by 4, so three levels are legal
ITERS = [40, 20, 10]            # per level, fine first
THREADS = int(os.environ.get('ICE_MG_THREADS', '4'))
RANKS = int(os.environ.get('ICE_MG_RANKS', '1'))

# The serial and the MPI registration run this same directory, so they must not
# share a scratch tree or ctest -j will have them overwrite each other's cases.
WORK = Path(__file__).resolve().parent / ('work' if RANKS == 1 else 'work-r%d' % RANKS)


def write_tec_multi(path, names, zones):
    """Multi-zone structured BLOCK file: nodal x,y,z then cell-centred data."""
    nvar = 3 + len(names)
    out = [' VARIABLES ="x" "y" "z" ' + ' '.join('"%s"' % v for v in names)]
    for zone, xn, yn, zn, cellvars in zones:
        out.append(' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, '
                   'VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED)'
                   % (zone, len(xn), len(yn), len(zn), nvar))
        out += [_fmt(_nodal(xn, yn, zn, c)) for c in (0, 1, 2)]
        out += [_fmt(v) for v in cellvars]
    Path(path).write_text('\n'.join(out) + '\n')


def write_bc_chain(path, nb, nx, ny, nz):
    """Blocks closed into a ring along i: faces 1 and 2 connect to the neighbour,
    faces 3/4 extrapolate and the degenerate k faces carry the null code."""
    rows = []

    def head(b, i, j, k, f, code):
        return '%8d%8d%8d%8d%8d%8d' % (b, i, j, k, f, code)

    def conn(bs, i_s, j, k, fs):
        return '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (bs, i_s, j, k, fs, 1, 0, 0, 1)

    for b in range(1, nb + 1):
        across = ((1, 1, (b - 2) % nb + 1, nx, 2),      # face 1 <- neighbour's face 2
                  (2, nx, b % nb + 1, 1, 1))            # face 2 <- neighbour's face 1
        for f, i0, bs, i_s, fs in across:
            for k in range(1, nz + 1):
                for j in range(1, ny + 1):
                    rows.append(head(b, i0, j, k, f, 201))
                    rows.append(conn(bs, i_s, j, k, fs))
        for f, j0 in ((3, 1), (4, ny)):
            for k in range(1, nz + 1):
                for i in range(1, nx + 1):
                    rows.append(head(b, i, j0, k, f, 400))
        for f, k0 in ((5, 1), (6, nz)):
            for j in range(1, ny + 1):
                for i in range(1, nx + 1):
                    rows.append(head(b, i, j, k0, f, 0))

    Path(path).write_text('\n'.join(rows) + '\n')


def sine(amp, base):
    return lambda x: base * (1.0 + amp * math.sin(2.0 * math.pi * x))


def build(work, levels):
    """A two-block ring carrying a sine, with a boundary table for every level."""
    work = Path(work)
    if work.exists():
        shutil.rmtree(str(work))
    (work / 'INPUT').mkdir(parents=True)
    (work / 'OUTPUT').mkdir()

    Lx, Ly = 1.0, 0.125
    names = ['rp1', 'up1', 'vp1', 'wp1', 'Tp1', 'np1']
    fields = [sine(0.30, 1.0e-3), sine(0.20, 10.0), lambda x: 0.0,
              lambda x: 0.0, sine(0.10, 300.0), sine(0.30, 1.0e6)]
    zones = []
    for b in range(NB):
        xn = linspace(b * Lx / NB, (b + 1) * Lx / NB, NXB + 1)
        yn = linspace(0.0, Ly, NY + 1)
        xc = [0.5 * (xn[i] + xn[i + 1]) for i in range(NXB)]
        zones.append(('B%d' % (b + 1), xn, yn, [0.0, Ly / NY],
                      [[f(x) for _ in range(NY) for x in xc] for f in fields]))
    write_tec_multi(work / 'INPUT/part-ic.tec', names, zones)

    for m in range(1, levels + 1):
        r = 2 ** (m - 1)
        name = 'part-bc.txt' if m == 1 else 'part-bc%d.txt' % m
        write_bc_chain(work / 'INPUT' / name, NB, NXB // r, NY // r, 1)

    sections = {
        'ICE-Parameters': {'iter-threshold': ITERS[0], 'time-threshold': 1e30,
                           'res-threshold': 0.0},
        'ICE-Numerics': {'time-scheme': 'RK2', 'cfl': 0.8, 'time-accurate': False,
                         'space-reconstruction': 'MUSCL', 'flux-limiter': 'vanleer'},
        'ICE-Family1': {'closure': 'MK'},
        'ICE-Physics': {'drag': 'Stokes', 'heat-transfer': 'Stokes',
                        'density': 1000.0, 'specific-heat': 900.0, 'emissivity': 0.0},
        'ICE-IO': {'ic-format': 'tecplot ascii', 'sol-format': 'tecplot ascii',
                   'shell-diter': 1000000000, 'sol-diter': 1000000000,
                   'res-diter': 1000000000},
    }
    if levels > 1:
        sections['ICE-Multigrid'] = dict(
            [('levels', levels)] +
            [('level%d-iter' % (m + 1), ITERS[m]) for m in range(levels)])
    write_ini(work / 'input.ini', sections)
    return work


def run(work, ranks=1, threads=1):
    if not ICE_BIN.exists():
        raise SystemExit('[fast] no ICE binary at %s -- build first' % ICE_BIN)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), KMP_STACKSIZE='100M',
               OMP_STACKSIZE='100M')
    cmd = [str(ICE_BIN)] if ranks == 1 else ['mpirun', '-n', str(ranks), str(ICE_BIN)]
    with open(str(work / 'log'), 'w') as out, open(str(work / 'err'), 'w') as err:
        rc = subprocess.call(cmd, cwd=str(work), stdout=out, stderr=err, env=env)
    if rc != 0:
        raise RuntimeError('ICE failed (exit %d) in %s -- see log and err there' % (rc, work))
    return read_solution(work / 'OUTPUT/part-field.tec')


def finite(sol):
    return all(v == v and abs(v) != float('inf') for var in sol['var'] for v in var)


def worst_rel(a, b):
    """Largest |a-b| over all variables, each scaled by that variable's range."""
    out = 0.0
    for va, vb in zip(a['var'], b['var']):
        span = max(abs(x) for x in va) or 1.0
        out = max(out, max(abs(x - y) for x, y in zip(va, vb)) / span)
    return out


def main():
    r = Report('Multigrid')
    fields = {}
    for levels in (1, 2, 3):
        sol = run(build(WORK / ('L%d' % levels), levels))
        fields[levels] = sol
        r.check(finite(sol), 'MG-levels=%d reached the fine level, field finite' % levels)

    for levels in (2, 3):
        d = worst_rel(fields[1], fields[levels])
        r.check(0.0 < d < 0.2,
                'MG-levels=%d moved the answer by %.3g of the field range '
                '(0 would mean the coarse levels never reached it)' % (levels, d))

    base = fields[3]
    sol = run(build(WORK / 'threads', 3), threads=THREADS)
    r.check(base['var'] == sol['var'],
            '%d threads give a bit-identical MG-levels=3 result' % THREADS)

    if RANKS > 1:
        work = build(WORK / 'ranks', 3)
        sol = run(work, ranks=RANKS)
        log = (work / 'log').read_text()
        if 'Number of ranks   --> %4d' % RANKS not in log:
            r.check(False, '%s did not run as %d MPI ranks -- is it the MPI build?'
                    % (ICE_BIN, RANKS))
        else:
            r.check(base['var'] == sol['var'],
                    '%d ranks give a bit-identical MG-levels=3 result' % RANKS)

    r.close(WORK)


if __name__ == '__main__':
    main()
