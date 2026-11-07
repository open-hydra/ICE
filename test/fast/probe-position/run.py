"""Fast test: a probe given by position takes the cell whose centre is nearest; a probe takes the names of the
solution files, the pressures included, as many as its key holds.

Two probes in each case, one by position and one by index in the cell the position names; 8 iterations, a row every
2. The position probe's values must be finite, its rows those of the index probe, and its last row the values of
that cell in the final solution (to 1e-13: a list-directed row need not carry every digit).

  stream   the 1-D stream of the restart tests (32 cells, one layer in y and z), the probe at x = 0.155 on the lower
           faces (y = z = 0), as a 1-D case is usually given: cell 16
  ring     the two-block ring of fast/multigrid, the probe at the centre x of cell 16 of block 2 and y = 0.06, z = 0:
           cell (16, 4) of block 2, found in the second block

RED while a position took the nearest node, the j = 0 and k = 0 planes included, and the node's indices for those of
a cell: the stream's probe read a ghost cell (NaN), the ring's the neighbouring one.

  names-IG     an IG stream of 32 cells, two probes by index on cell 16: every name of the solution file with _p1, and
               every one with _1 (P_p1, P_1 the pressure)
  names-AG     the same for AG (P11_p1 ... P33_p1, P11_1 ... P33_1 the tensor)
  names-AG-24  one AG probe with the 24 names of both spellings

Each pressure component differs from the others (IG P = 2 + x; AG P11 = 3 + x, P12 = 0.3 + 0.1x, P13 = 0.2 + 0.1x,
P22 = 2.5 + x, P23 = 0.1 + 0.1x, P33 = 2 + x, positive definite), so a name read from another column shows; the shipped
Axis cases are uniform in P (1) and in every off-diagonal (0). 8 steady iterations, a row every 2; the last row of each
file must be the named columns of cell 16 in the solution (1e-13). RED while the probes knew no pressure (set-up stop:
unknown variable "P_p1", "P11_p1") and while the names went into an array of 16 that parse (ORION) fills unchecked.
"""
import importlib.util
import math
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import Case, ICE_BIN, GREEN, RED, RESET, read_solution     # noqa: E402

HERE = Path(__file__).resolve().parent
WORK = HERE / 'work'
VARS = ['rho_p1', 'u_p1', 'T_p1']
PRESSURES = {'IG': ['P'], 'AG': ['P11', 'P12', 'P13', 'P22', 'P23', 'P33']}
IC_NAMES = {'IG': ['rp1', 'up1', 'vp1', 'wp1', 'Pp1', 'Tp1', 'np1'],
            'AG': ['rp1', 'up1', 'vp1', 'wp1', 'P11', 'P12', 'P13', 'P22', 'P23', 'P33', 'Tp1', 'np1']}


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


def probes(d, position, index):
    ini_set(d, 'ICE-Parameters', iter_threshold=8)
    ini_set(d, 'ICE-Probes', probe1='bypos', probe2='byidx')
    ini_set(d, 'bypos', variables=' '.join(VARS), position=position, diter=2)
    ini_set(d, 'byidx', variables=' '.join(VARS), index_position=index, diter=2)


def stream(work):
    mass = 1000.0 * 4.0 / 3.0 * math.pi * 1.0e-15
    case = Case(work, nx=32, Lx=0.32)
    rho = lambda x: 1.0 + 0.5 * math.sin(2.0 * math.pi * x / 0.32)                  # noqa: E731
    case.particles(rho, 10.0, 0.0, 0.0, 300.0, lambda x: rho(x) / mass)
    case.boundaries(inlets=[{'gp': 10.0, 'v': 10.0, 'T': 300.0, 'r': 1.0e-5}])
    case.ini(cfl=0.8, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=1000.0, numerics={'time-accurate': False})
    probes(case.dir, '0.155 0.0 0.0', '1 16 1 1')
    return case.dir


def ring(work):
    spec = importlib.util.spec_from_file_location('mg_run', str(HERE.parent / 'multigrid' / 'run.py'))
    mg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mg)
    d = mg.build(work, 1)
    probes(d, '0.7421875 0.06 0.0', '2 16 4 1')
    return Path(d)


def rows(path, nvar=len(VARS)):
    """The rows of a probe file, an iteration or a time and then the nvar variables; written list-directed, a row may
    take more than one line, so the values are counted."""
    vals = [float(v) for v in Path(path).read_text().split()]
    n = 1 + nvar
    return [vals[k:k + n] for k in range(0, len(vals), n)]


def names_stream(work, closure, lists):
    """The 1-D stream of 32 cells with one family of the closure, its pressure components distinct, and one probe by
    index on cell 16 per name list."""
    mass = 1000.0 * 4.0 / 3.0 * math.pi * 1.0e-15
    case = Case(work, nx=32, Lx=0.32)
    rho = lambda x: 1.0 + 0.5 * math.sin(2.0 * math.pi * x / 0.32)                  # noqa: E731
    if closure == 'IG':
        pressure = [lambda x: 2.0 + x]
    else:
        pressure = [lambda x: 3.0 + x, lambda x: 0.3 + 0.1 * x, lambda x: 0.2 + 0.1 * x,
                    lambda x: 2.5 + x, lambda x: 0.1 + 0.1 * x, lambda x: 2.0 + x]
    case.write_ic(IC_NAMES[closure], [rho, 10.0, 0.0, 0.0] + pressure + [300.0, lambda x: rho(x) / mass])
    case.boundaries('extrapolation')
    case.ini(cfl=0.8, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=1000.0, numerics={'time-accurate': False},
             closures=(closure,))
    d = case.dir
    ini_set(d, 'ICE-Parameters', iter_threshold=8)
    ini_set(d, 'ICE-Probes', **{'probe%d' % (m + 1): 'names%d' % (m + 1) for m in range(len(lists))})
    for m, names in enumerate(lists):
        ini_set(d, 'names%d' % (m + 1), variables=' '.join(names), index_position='1 16 1 1', diter=2)
    return d


def names_leg(name, closure, lists):
    """Run the stream; each probe's last row must be its named columns of cell 16 in the solution (a name _1 is the
    solution's _p1)."""
    d = names_stream(WORK / name, closure, lists)
    env = dict(os.environ, OMP_NUM_THREADS='1', KMP_STACKSIZE='100M', OMP_STACKSIZE='100M')
    with open(str(d / 'log'), 'w') as out:
        rc = subprocess.call([str(ICE_BIN)], cwd=str(d), stdout=out, stderr=subprocess.STDOUT, env=env)
    if rc != 0:
        print('  %-11s %s  ICE failed (exit %d), see %s/log' % (name, RED + 'FAIL' + RESET, rc, d))
        return 1
    sol = read_solution(d / 'OUTPUT/part-field.tec')
    ok, notes = True, []
    for m, names in enumerate(lists):
        got = rows(d / ('OUTPUT/names%d.txt' % (m + 1)), len(names))
        want = [sol['var'][sol['names'].index(v if '_p' in v else v.replace('_', '_p'))][15] for v in names]
        last = bool(got) and len(got[-1]) == 1 + len(names) and \
            all(a == b or abs(a - b) <= 1.0e-13 * abs(b) for a, b in zip(got[-1][1:], want))
        ok = ok and last
        notes.append('names%d: %d names, %d rows, last row = cell 16 of the solution %s'
                     % (m + 1, len(names), len(got), last))
        if not last:
            notes.append('%s against %s' % (got[-1][1:] if got else None, want))
    print('  %-11s %s  %s' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET), '; '.join(notes)))
    return 0 if ok else 1


def leg(name, make, zone, cell):
    d = make(WORK / name)
    env = dict(os.environ, OMP_NUM_THREADS='1', KMP_STACKSIZE='100M', OMP_STACKSIZE='100M')
    with open(str(d / 'log'), 'w') as out:
        rc = subprocess.call([str(ICE_BIN)], cwd=str(d), stdout=out, stderr=subprocess.STDOUT, env=env)
    if rc != 0:
        print('  %-7s %s  ICE failed (exit %d), see %s/log' % (name, RED + 'FAIL' + RESET, rc, d))
        return 1
    pos, idx = rows(d / 'OUTPUT/bypos.txt'), rows(d / 'OUTPUT/byidx.txt')
    sol = read_solution(d / 'OUTPUT/part-field.tec', zone=zone)
    want = [sol['var'][sol['names'].index(v)][cell] for v in VARS]
    finite = bool(pos) and all(math.isfinite(v) for r in pos for v in r)
    same = pos == idx
    last = bool(pos) and len(pos[-1]) == 1 + len(VARS) and \
        all(abs(a - b) <= 1.0e-13 * abs(b) for a, b in zip(pos[-1][1:], want))
    ok = finite and same and last
    print('  %-7s %s  %d rows; finite %s, the index probe\'s rows %s, last row = the cell in the solution %s '
          '(%s against %s)' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET), len(pos), finite,
                               same, last, pos[-1][1:] if pos else None, want))
    return 0 if ok else 1


def main():
    fails = 0
    fails += leg('stream', stream, 0, 15)              # cell 16 of 32: 0-based 15
    fails += leg('ring', ring, 1, 3 * 32 + 15)         # block 2, cell (16, 4) of 32 x 8
    for closure in ('IG', 'AG'):
        names = ['rho', 'u', 'v', 'w'] + PRESSURES[closure] + ['T', 'n']
        lists = [[v + '_p1' for v in names], [v + '_1' for v in names]]
        fails += names_leg('names-' + closure, closure, lists)
    fails += names_leg('names-AG-24', 'AG', [sum(lists, [])])        # the AG lists of the loop's last pass
    print('[fast] %s: a probe given by position takes the cell whose centre is nearest; a probe takes the names of the '
          'solution files' % ('PASS' if fails == 0 else 'FAIL'))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
