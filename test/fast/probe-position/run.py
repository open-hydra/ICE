"""Fast test: a probe given by position takes the cell whose centre is nearest.

Two probes in each case, one by position and one by index in the cell the position names; 8 iterations, a row every
2. The position probe's values must be finite, its rows those of the index probe, and its last row the values of
that cell in the final solution (to 1e-13: a list-directed row need not carry every digit).

  stream   the 1-D stream of the restart tests (32 cells, one layer in y and z), the probe at x = 0.155 on the lower
           faces (y = z = 0), as a 1-D case is usually given: cell 16
  ring     the two-block ring of fast/multigrid, the probe at the centre x of cell 16 of block 2 and y = 0.06, z = 0:
           cell (16, 4) of block 2, found in the second block

RED while a position took the nearest node, the j = 0 and k = 0 planes included, and the node's indices for those of
a cell: the stream's probe read a ghost cell (NaN), the ring's the neighbouring one.
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


def rows(path):
    """The rows of a probe file, an iteration or a time and then the variables; written list-directed, a row may take
    more than one line, so the values are counted."""
    vals = [float(v) for v in Path(path).read_text().split()]
    n = 1 + len(VARS)
    return [vals[k:k + n] for k in range(0, len(vals), n)]


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
    print('[fast] %s: a probe given by position takes the cell whose centre is nearest'
          % ('PASS' if fails == 0 else 'FAIL'))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
