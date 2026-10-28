"""Fast test: a restart goes on from the time of the solution it reads.

  tecplot-steady  a 1-D stream carrying a density wave (MK, MUSCL, RK2), steady: run A of N = 20 iterations
                  writes part-field.tec, restart B (newrun = false) of n = 10 iterations reads it. B's rows of
                  part-residual-history.dat are the iterations N + 1 .. N + n, and the solution time of its
                  part-field.tec is -(N + n): a steady solution holds minus its global iteration count. RED when
                  the restart takes the time with its sign for the count: B goes on from iteration -N.
  vtk-steady      the same with sol-format = vtk raw: B reads OUTPUT/part-field.vtm and the block files it lists
                  in OUTPUT/vtk/; rows N + 1 .. N + n, the TIME of B's solution -(N + n).
  vtk-transient   time-accurate, vtk raw: B goes on in time from the TIME of A's solution. B counts its
                  iterations from 1 again, its first row's time is above A's last one by less than 1.5 of A's
                  last step, its times increase, and the TIME of its solution is its last row's time.
"""
import math
import re
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import Case, run_ice, GREEN, RED, RESET                    # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
N, NB = 20, 10
VINJ, TINJ, RINJ, RHO_AL = 10.0, 300.0, 1.0e-5, 1000.0
MASS = RHO_AL * 4.0 / 3.0 * math.pi * RINJ**3


def set_keys(case, **keys):
    """Overwrite or add keys of the case's input.ini (Case.ini writes iter-threshold, sol-format, res-diter)."""
    path = case.dir / 'input.ini'
    lines = path.read_text().split('\n')
    done = set()
    for i, line in enumerate(lines):
        key = line.split('=')[0].strip()
        if key in keys:
            lines[i] = '%s = %s' % (key, keys[key])
            done.add(key)
    for key in keys:
        if key not in done:
            k = lines.index('[ICE-Parameters]')
            lines.insert(k + 1, '%s = %s' % (key, keys[key]))
    path.write_text('\n'.join(lines))


def history(case):
    """[(iteration, time)] of part-residual-history.dat, one row per iteration (res-diter = 1)."""
    rows = []
    for line in (case.dir / 'OUTPUT' / 'part-residual-history.dat').read_text().split('\n'):
        f = line.split()
        if len(f) >= 3:
            rows.append((int(f[0]), float(f[1])))
    return rows


def tec_time(case):
    """The SOLUTIONTIME of the first zone header of OUTPUT/part-field.tec, None when it has none."""
    m = re.search(r'SOLUTIONTIME\s*=\s*([-+0-9.EDed]+)', (case.dir / 'OUTPUT' / 'part-field.tec').read_text())
    return float(m.group(1)) if m else None


def vtk_time(case):
    """The TIME field data of the first block file that OUTPUT/part-field.vtm lists (raw appended data: a 4-byte
    size, then the value), None when the file holds none."""
    vtm = case.dir / 'OUTPUT' / 'part-field.vtm'
    block = re.search(r'file="([^"]+\.vts)"', vtm.read_text()).group(1)
    data = (vtm.parent / block).read_bytes()
    head = data.find(b'<AppendedData')
    m = re.search(rb'Name="TIME"[^>]*offset="(\d+)"', data[:head])
    if m is None:
        return None
    start = data.find(b'_', head) + 1 + int(m.group(1))
    return struct.unpack('<d', data[start + 4:start + 12])[0]


def report(name, ok, detail):
    print('  %-15s %s  %s' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET), detail))
    return 0 if ok else 1


def stream(work, steady, fmt):
    """Run A (N iterations) and restart B (NB iterations) of the stream; the rows of A and those B added."""
    case = Case(work, nx=32, Lx=0.32)
    rho = lambda x: 1.0 + 0.5 * math.sin(2.0 * math.pi * x / 0.32)                  # noqa: E731
    case.particles(rho, VINJ, 0.0, 0.0, TINJ, lambda x: rho(x) / MASS)
    case.boundaries(inlets=[{'gp': VINJ, 'v': VINJ, 'T': TINJ, 'r': RINJ}])
    case.ini(iters=N, cfl=0.8, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=RHO_AL,
             numerics={'time-accurate': not steady})
    set_keys(case, **{'res-diter': 1, 'sol-format': fmt})
    try:
        run_ice(case.dir)
    except RuntimeError:
        return case, None, None, None, None
    rows_a = history(case)
    time_a = tec_time(case) if fmt.startswith('tecplot') else vtk_time(case)
    set_keys(case, **{'iter-threshold': NB, 'newrun': 'false'})
    try:
        run_ice(case.dir)
    except RuntimeError:
        return case, rows_a, time_a, None, None
    rows = history(case)
    time_b = tec_time(case) if fmt.startswith('tecplot') else vtk_time(case)
    return case, rows_a, time_a, rows[len(rows_a):], time_b


def steady_leg(name, fmt):
    case, rows_a, time_a, rows_b, time_b = stream(WORK / name, True, fmt)
    if rows_b is None:
        return report(name, False, 'run A or the restart failed (see %s/log)' % case.dir)
    its = [r[0] for r in rows_b]
    ok = its == list(range(N + 1, N + NB + 1)) and time_a == -N and time_b == -(N + NB)
    return report(name, ok, 'A wrote %d rows, time %r; B rows %s..%s (%d), time %r (want %d..%d, %r)'
                  % (len(rows_a), time_a, its[0] if its else None, its[-1] if its else None, len(its), time_b,
                     N + 1, N + NB, -(N + NB)))


def transient_leg(name, fmt):
    case, rows_a, time_a, rows_b, time_b = stream(WORK / name, False, fmt)
    if rows_b is None:
        return report(name, False, 'run A or the restart failed (see %s/log)' % case.dir)
    its = [r[0] for r in rows_b]
    times = [r[1] for r in rows_b]
    last_step = rows_a[-1][1] - rows_a[-2][1]
    ok = (its == list(range(1, NB + 1)) and time_a is not None and times and times[0] > rows_a[-1][1]
          and times[0] - rows_a[-1][1] <= 1.5 * last_step and all(b > a for a, b in zip(times, times[1:]))
          and time_b is not None and abs(time_b - times[-1]) <= 1e-9 * abs(times[-1]))
    return report(name, ok, 'A ends at %.10e (TIME %r); B rows %s..%s, times %.10e .. %.10e, TIME %r'
                  % (rows_a[-1][1], time_a, its[0] if its else None, its[-1] if its else None,
                     times[0] if times else float('nan'), times[-1] if times else float('nan'), time_b))


def main():
    fails = 0
    fails += steady_leg('tecplot-steady', 'tecplot ascii')
    fails += steady_leg('vtk-steady', 'vtk raw')
    fails += transient_leg('vtk-transient', 'vtk raw')
    print('[fast] %s: a restart goes on from the time of its solution' % ('PASS' if fails == 0 else 'FAIL'))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
