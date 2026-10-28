"""Fast test: a restart goes on as the uninterrupted run would.

Each leg runs a case twice: A, 2N iterations in one go; B, N iterations, then a restart (newrun = false) of N more.
The OUTPUT directories of A and B must hold the same files with the same bytes: the final solution and the numbered
series of sol-overwrite = false, the diagnostic files, the residual history and the probe rows. The legs reach every
counter that a solution stores for its restart (CYCLE, ICE_ITERATION, ICE_OUTPUT_NUMBER, ICE_GRID_LEVEL,
ICE_LAST_OUTPUT_TIME, ICE_RES_REFERENCE): the restart falls inside the CFL ramp (cfl-rise-threshold), between two
numbered solutions and on a residual row (res-diter divides N: every run writes a row at its last iteration).

  steady-tecplot     a 1-D stream carrying a density wave (MK, MUSCL, RK2), steady, tecplot ascii: solutions on
                     sol-diter, res-threshold on (a reference that never stops the run), a probe on diter
  steady-vtk         the same, vtk raw
  transient-tecplot  the stream time-accurate: solutions on sol-dtime, a probe on diter and dtime, tecplot ascii
  transient-vtk      the same, vtk raw
  multigrid          the two-block ring of fast/multigrid at 2 levels, steady, tecplot ascii, a probe in block 2:
                     the restart falls in the fine phase, and goes on there instead of running the coarse level again
  res-threshold      the steady stream without the ramp, whose density residual falls at every iteration, and a
                     res-threshold that stops A at an iteration K after the restart (from a run A0 that does not
                     stop): B must stop at K too, which it does only measured against the residual of A's first
                     iteration
  wedge              the one-degree wedge of Axis/AG (axis row, side faces 200), its AG cloud set moving along x:
                     the case whose w and P13, P23 a mirror freezes, steady, tecplot ascii; one side node plane
                     moved by a relative 1e-12, so that a restart which does not freeze them fails
  past-steady-tec  B goes on past its solution of iteration N to N + N/2 and stops (killed, or at the time
                     limit of a queue); its numbered solution of iteration N is put back as the one a restart reads
                     and the restart runs N more: the rows of N + 1 .. N + N/2 must go first (the residual rows and
                     the iteration rows of the probe), steady stream, solutions on sol-diter
  past-transient-vtk the same, time-accurate, vtk raw: the probe's time rows after the solution go too

ICE_RESTART_RANKS=<n> runs every run on n MPI ranks (an MPI build; the ring's blocks, and its probe, then go to
different ranks).

Without the counters every leg fails: the restart ramps the CFL again, numbers its solutions from 1, writes its
time-based solutions on its own clock, counts the residual rows of a time-accurate run from 1, runs the coarse level
again, measures res-threshold against its own first residual, and from a vtk solution drops no row.
"""
import importlib.util
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import Case, ICE_BIN, GREEN, RED, RESET                    # noqa: E402

HERE = Path(__file__).resolve().parent
TEST = HERE.parents[1]
RANKS = int(os.environ.get('ICE_RESTART_RANKS', '1'))
# The serial and the MPI registration run this same directory: separate scratch trees for ctest -j
WORK = HERE / ('work' if RANKS == 1 else 'work-r%d' % RANKS)
N = 12                                   # iterations of each part of B; res-diter, diter divide it
VINJ, TINJ, RINJ, RHO_AL = 10.0, 300.0, 1.0e-5, 1000.0
MASS = RHO_AL * 4.0 / 3.0 * math.pi * RINJ**3


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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


def stream(work, steady, fmt, ramp=True):
    """The stream of the restart test: 32 cells, a density wave entering at 10 m/s."""
    case = Case(work, nx=32, Lx=0.32)
    rho = lambda x: 1.0 + 0.5 * math.sin(2.0 * math.pi * x / 0.32)                  # noqa: E731
    case.particles(rho, VINJ, 0.0, 0.0, TINJ, lambda x: rho(x) / MASS)
    case.boundaries(inlets=[{'gp': VINJ, 'v': VINJ, 'T': TINJ, 'r': RINJ}])
    case.ini(iters=2 * N, cfl=0.8, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=RHO_AL, dt_max=1.0,
             numerics={'time-accurate': not steady, 'cfl-rise-threshold': 3 * N // 2 if ramp else 0})
    io = {'sol-format': fmt, 'sol-overwrite': 'false', 'res-diter': 3}
    # the cell in the middle of the stream
    probe = {'variables': 'rho_p1 u_p1', 'index_position': '1 16 1 1', 'diter': 4}
    if steady:
        io['sol-diter'] = 5
        ini_set(case.dir, 'ICE-Parameters', res_threshold='1e-300')
    else:
        # time-based solutions at iterations 9 and 13, the restart at 12: the next one falls on the clock of A
        io['sol-dtime'] = 1.8e-3
        probe['dtime'] = 1.0e-3
    ini_set(case.dir, 'ICE-IO', **io)
    ini_set(case.dir, 'ICE-Probes', probe1='centre')
    ini_set(case.dir, 'centre', **probe)
    return case.dir


def ring(work):
    """The two-block ring of fast/multigrid at 2 levels: N coarse iterations, then 2N fine ones."""
    d = load('mg_run', HERE.parent / 'multigrid' / 'run.py').build(work, 2)
    ini_set(d, 'ICE-Parameters', iter_threshold=2 * N)
    ini_set(d, 'ICE-Multigrid', level1_iter=2 * N, level2_iter=N)
    ini_set(d, 'ICE-Numerics', cfl_rise_threshold=3 * N // 2)
    ini_set(d, 'ICE-IO', sol_overwrite='false', sol_diter=5, res_diter=3)
    ini_set(d, 'ICE-Probes', probe1='block2')
    ini_set(d, 'block2', variables='rho_p1 T_p1', index_position='2 1 1 1', diter=4)
    return Path(d)


def wedge(work):
    """Axis/AG's case (its make_case.py and input.ini) with u = 2 sin(2 pi x / Lx) instead of the cloud at rest, and
    its +delthe/2 node plane moved by a relative 1e-12. The setup still freezes w, P13 and P23 (the move is far
    inside mesh_fit), but where they are integrated the move seeds them with round-off, as the nodes of a mesh
    written by a mesher do, and their mode across the layer grows: a restart that does not freeze them fails. On the
    exact wedge the two side faces cancel to the bit and w stays 0 even when integrated, so the leg could not tell."""
    import numpy as np
    mc = load('axis_ag', TEST / 'Axis' / 'AG' / 'make_case.py')
    d = Path(work)
    if d.exists():
        shutil.rmtree(str(d))
    for sub in ('INPUT', 'MESH', 'OUTPUT'):
        (d / sub).mkdir(parents=True)
    Rp, rho_m, rho, P, T = 1e-5, 2000.0, 1.0, 1.0, 300.0
    X, Y, Z = mc.wedge(nx=40, nr=20, Lx=0.1, R=0.05, deg=1.0, r0=1e-8)
    Z[:, :, 1] *= 1.0 + 1e-12                        # the +delthe/2 node plane, off by a relative 1e-12
    xc = mc.cc(X, Y, Z)[0]
    one = np.ones(xc.shape)
    n = rho / (4.0 / 3.0 * math.pi * Rp**3 * rho_m)
    mc.write_tec(str(d / 'INPUT/part-ic.tec'),
                 ['x', 'y', 'z', 'rho', 'u', 'v', 'w', 'P11', 'P12', 'P13', 'P22', 'P23', 'P33', 'T', 'n'], X, Y, Z,
                 [rho * one, 2.0 * np.sin(2.0 * math.pi * xc / 0.1), 0 * one, 0 * one, P * one, 0 * one, 0 * one,
                  P * one, 0 * one, P * one, T * one, n * one])
    mc.write_bc(str(d / 'INPUT/part-bc.txt'), *xc.shape, {1: 300, 2: 300, 3: 300, 4: 300, 5: 200, 6: 200})
    (d / 'INPUT/part-phase.txt').write_text('condensed-dispersed phase\nA 1\n')
    mc.write_tec(str(d / 'MESH/mesh.tec'), ['x', 'y', 'z'], X, Y, Z, [])
    shutil.copy(str(TEST / 'Axis' / 'AG' / 'input.ini'), str(d / 'input.ini'))
    ini_set(d, 'ICE-Parameters', iter_threshold=2 * N)
    ini_set(d, 'ICE-Numerics', cfl_rise_threshold=3 * N // 2)
    ini_set(d, 'ICE-IO', sol_overwrite='false', sol_diter=5, res_diter=3, shell_diter=N)
    return d


def run(case_dir):
    """Run the solver in case_dir, on RANKS ranks. Raises if it fails, pointing at the log."""
    env = dict(os.environ, OMP_NUM_THREADS='1', KMP_STACKSIZE='100M', OMP_STACKSIZE='100M')
    cmd = [str(ICE_BIN)] if RANKS == 1 else ['mpirun', '-n', str(RANKS), str(ICE_BIN)]
    with open(str(Path(case_dir) / 'log'), 'w') as out:
        rc = subprocess.call(cmd, cwd=str(case_dir), stdout=out, stderr=subprocess.STDOUT, env=env)
    if rc != 0:
        raise RuntimeError('ICE failed (exit %d) in %s -- see the log there' % (rc, case_dir))


def rows(case_dir):
    """[(iteration, density residual)] of the residual history."""
    out = []
    for line in (Path(case_dir) / 'OUTPUT' / 'part-residual-history.dat').read_text().split('\n'):
        f = line.split()
        if len(f) >= 3:
            out.append((int(f[0]), float(f[2])))
    return out


def tree(d):
    """{relative path: bytes} of every file under d/OUTPUT."""
    out = d / 'OUTPUT'
    return {str(p.relative_to(out)): p.read_bytes() for p in sorted(out.rglob('*')) if p.is_file()}


def report(name, ok, detail):
    print('  %-18s %s  %s' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET), detail))
    return 0 if ok else 1


def leg(name, make, stop=None):
    """A against B; with stop, the iteration at which both must have stopped (the last row of each history)."""
    a, b = make(WORK / name / 'A'), make(WORK / name / 'B')
    ini_set(b, 'ICE-Parameters', iter_threshold=N)
    try:
        run(a)
        run(b)
        ini_set(b, 'ICE-Parameters', newrun='false')
        run(b)
    except RuntimeError as e:
        return report(name, False, str(e))
    return compare(name, a, b, stop)


def went_past_leg(name, steady, fmt):
    """B runs N + N/2 iterations, its numbered solution of iteration N (number 2: solutions every N/2) is copied
    over the final one, and the restart runs N more."""
    ext = '.tec' if fmt.startswith('tecplot') else '.vtm'     # a .vtm names its own block files

    def make(w):
        d = stream(w, steady, fmt)
        ini_set(d, 'ICE-IO', sol_diter=N // 2, sol_dtime='1e30')
        return d
    a, b = make(WORK / name / 'A'), make(WORK / name / 'B')
    ini_set(b, 'ICE-Parameters', iter_threshold=N + N // 2)
    try:
        run(a)
        run(b)
        shutil.copy(str(b / 'OUTPUT' / ('part-field2' + ext)), str(b / 'OUTPUT' / ('part-field' + ext)))
        ini_set(b, 'ICE-Parameters', iter_threshold=N, newrun='false')
        run(b)
    except RuntimeError as e:
        return report(name, False, str(e))
    return compare(name, a, b)


def compare(name, a, b, stop=None):
    """The OUTPUT of A against that of B, file by file."""
    ta, tb = tree(a), tree(b)
    only_a, only_b = sorted(set(ta) - set(tb)), sorted(set(tb) - set(ta))
    differ = [f for f in sorted(set(ta) & set(tb)) if ta[f] != tb[f]]
    ok = not (only_a or only_b or differ)
    detail = '%d files, the same bytes' % len(ta) if ok else \
        'only in A %s, only in B %s, differ %s' % (only_a[:4], only_b[:4], differ[:6])
    if stop is not None:
        ka, kb = rows(a)[-1][0], rows(b)[-1][0]
        ok = ok and ka == kb == stop
        detail = 'stopped at %d (A) and %d (B), threshold met at %d; %s' % (ka, kb, stop, detail)
    return report(name, ok, detail)


def res_threshold_leg(name):
    """res-threshold put between the residual ratio of an iteration K after the restart and the lowest ratio before
    it, from the density residual of every iteration of a run A0 that does not stop; A then stops at K."""
    a0 = stream(WORK / name / 'A0', True, 'tecplot ascii', ramp=False)
    ini_set(a0, 'ICE-Parameters', res_threshold=0)
    ini_set(a0, 'ICE-IO', res_diter=1)
    try:
        run(a0)
    except RuntimeError as e:
        return report(name, False, str(e))
    hist = rows(a0)
    first = next(r for _, r in hist if r > 0.0)
    ratio = [(k, r / first) for k, r in hist]
    for i, (k, q) in enumerate(ratio):
        low = min(p for _, p in ratio[:i]) if i else float('inf')
        if k > N + N // 2 and q < low:
            threshold = 0.5 * (q + low)
            break
    else:
        return report(name, False, 'no iteration after %d whose residual ratio is a new low in A0 (%s)'
                      % (N + N // 2, ratio[-3:]))

    def make(w):
        d = stream(w, True, 'tecplot ascii', ramp=False)
        ini_set(d, 'ICE-Parameters', res_threshold=repr(threshold))
        return d
    return leg(name, make, stop=k)


def main():
    fails = 0
    fails += leg('steady-tecplot', lambda w: stream(w, True, 'tecplot ascii'))
    fails += leg('steady-vtk', lambda w: stream(w, True, 'vtk raw'))
    fails += leg('transient-tecplot', lambda w: stream(w, False, 'tecplot ascii'))
    fails += leg('transient-vtk', lambda w: stream(w, False, 'vtk raw'))
    fails += leg('multigrid', ring)
    fails += res_threshold_leg('res-threshold')
    fails += leg('wedge', wedge)
    fails += went_past_leg('past-steady-tec', True, 'tecplot ascii')
    fails += went_past_leg('past-transient-vtk', False, 'vtk raw')
    print('[fast] %s: a run of 2N iterations and one of N restarted for N more write the same files (%d rank%s)'
          % ('PASS' if fails == 0 else 'FAIL', RANKS, '' if RANKS == 1 else 's'))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
