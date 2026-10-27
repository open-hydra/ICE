"""Fast test: input.ini edited during a run, with an error in it, under ini-error = warn and stop.

The case is a 1-D MK stream carrying a density wave, time-accurate, its step set by the CFL number (dt-max out of
the way); input.ini is read again every 50 iterations and the residual history is written at every one. Once the
first row of the history is there (set-up is over), the harness replaces input.ini, in one rename, by a copy whose
cfl is 0.4 instead of 0.8 (a valid change), whose sol-diter = seven does not read as an integer, and whose ini-error
names the other mode.

  warn      ini-error left out (warn): the readings after the edit are not taken, so the run goes on to its last
            iteration with exit code 0, warns, writes no [ERROR], and its residual history and field are those of the
            run without the edit, byte for byte. Taking the valid cfl of the edited file would change the step from
            that reading on; taking its ini-error = stop would stop the run at the next reading.
  stop      ini-error = stop in the file from the start: the first reading after the edit stops the run with the
            error of sol-diter and a non-zero exit code; the edited file's ini-error = warn does not apply to its own
            reading.
  empty     ini-error left out, input.ini replaced by comment lines only: not taken either (the file holds no key),
            the run goes on to its end with the field and history of the run without the edit.
"""
import math
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import Case, ICE_BIN, run_ice, GREEN, RED, RESET            # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
ITERS, INI_DITER = 8000, 50
VINJ, TINJ, RINJ, RHO_AL = 10.0, 300.0, 1.0e-5, 1000.0
MASS = RHO_AL * 4.0 / 3.0 * math.pi * RINJ**3
SOL_ERROR = '[ERROR] input.ini, [ICE-IO] sol-diter = seven: the value is not an integer'


def set_key(text, section, key, value):
    """text of an input.ini with key = value in [section]: the line of the key replaced, or added after the header."""
    lines, current, done = text.split('\n'), None, False
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith('['):
            current = s.strip('[]')
        elif current == section and s.split('=')[0].strip() == key:
            lines[i], done = '%s = %s' % (key, value), True
    if not done:
        k = lines.index('[%s]' % section)
        lines.insert(k + 1, '%s = %s' % (key, value))
    return '\n'.join(lines)


def build(work, ini_error=None):
    """The stream, with ini-diter, a history row per iteration and, when given, ini-error."""
    case = Case(work, nx=100, Lx=0.32)
    rho = lambda x: 1.0 + 0.5 * math.sin(2.0 * math.pi * x / 0.32)                  # noqa: E731
    case.particles(rho, VINJ, 0.0, 0.0, TINJ, lambda x: rho(x) / MASS)
    case.boundaries(inlets=[{'gp': VINJ, 'v': VINJ, 'T': TINJ, 'r': RINJ}])
    case.ini(iters=ITERS, cfl=0.8, rk='RK2', drag='NoDrag', heat='NoHeat', rho_al=RHO_AL, dt_max=1.0)
    path = case.dir / 'input.ini'
    text = set_key(path.read_text(), 'ICE-IO', 'res-diter', 1)
    text = set_key(text, 'ICE-IO', 'ini-diter', INI_DITER)
    if ini_error is not None:
        text = set_key(text, 'ICE-IO', 'ini-error', ini_error)
    path.write_text(text)
    return case


def edited(case, ini_error):
    """The edit: a valid cfl, a sol-diter that does not read, and the other mode."""
    text = (case.dir / 'input.ini').read_text()
    text = set_key(text, 'ICE-Numerics', 'cfl', 0.4)
    text = set_key(text, 'ICE-IO', 'sol-diter', 'seven')
    return set_key(text, 'ICE-IO', 'ini-error', ini_error)


def history(case):
    """[(iteration, row text)] of the residual history."""
    rows = []
    path = case.dir / 'OUTPUT' / 'part-residual-history.dat'
    if path.exists():
        for line in path.read_text().split('\n'):
            f = line.split()
            if f and f[0].isdigit():
                rows.append((int(f[0]), line))
    return rows


def run_with_edit(case, text):
    """Run the case; once the history has its first row, input.ini becomes text (one rename). (exit, edit iteration)."""
    env = dict(os.environ, OMP_NUM_THREADS='1', KMP_STACKSIZE='100M', OMP_STACKSIZE='100M')
    with open(str(case.dir / 'log'), 'w') as out, open(str(case.dir / 'err'), 'w') as err:
        proc = subprocess.Popen([str(ICE_BIN)], cwd=str(case.dir), stdout=out, stderr=err, env=env)
        at, start = None, time.time()
        while proc.poll() is None and time.time() - start < 300:
            rows = history(case)
            if rows:
                tmp = case.dir / 'input.ini.edit'
                tmp.write_text(text)
                os.replace(str(tmp), str(case.dir / 'input.ini'))
                at = rows[-1][0]
                break
            time.sleep(0.01)
        try:
            rc = proc.wait(timeout=300)
        except subprocess.TimeoutExpired:
            proc.kill()
            rc = proc.wait()
    return rc, at


def same_outputs(a, b):
    return all((a.dir / 'OUTPUT' / f).read_bytes() == (b.dir / 'OUTPUT' / f).read_bytes()
               if (a.dir / 'OUTPUT' / f).exists() and (b.dir / 'OUTPUT' / f).exists() else False
               for f in ('part-field.tec', 'part-residual-history.dat'))


def report(name, ok, detail):
    print('  %-7s %s  %s' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET), detail))
    return 0 if ok else 1


def not_taken(case, ref, name, text, warning):
    rc, at = run_with_edit(case, text)
    log = (case.dir / 'log').read_text(errors='replace')
    warned = log.count('is not taken: ')
    if at is None:
        return report(name, False, 'exit %d before the edit (no history row)' % rc)
    if rc == 0 and warned == 0:
        return report(name, False, 'no reading after the edit at iteration %d: harness late?' % at)
    last = history(case)[-1][0] if history(case) else 0
    ok = rc == 0 and warning in log and '[ERROR]' not in log and last == ITERS and same_outputs(case, ref)
    return report(name, ok, 'edit at iteration %d, exit %d, %d readings not taken, last iteration %d, outputs %s'
                  % (at, rc, warned, last, 'as without the edit' if same_outputs(case, ref) else 'DIFFER'))


def stop(case, text):
    rc, at = run_with_edit(case, text)
    log = (case.dir / 'log').read_text(errors='replace')
    if at is None:
        return report('stop', False, 'exit %d before the edit (no history row)' % rc)
    last = history(case)[-1][0] if history(case) else 0
    if rc == 0:
        return report('stop', False, 'exit 0 at iteration %d (edit at %d): harness late?' % (last, at))
    ok = SOL_ERROR in log and at <= last < ITERS
    return report('stop', ok, 'edit at iteration %d, exit %d at iteration %d, %s'
                  % (at, rc, last, 'error given' if SOL_ERROR in log else 'error NOT given'))


if __name__ == '__main__':
    print('input.ini edited during the run with an error, ini-error = warn | stop')
    ref = build(WORK / 'reference')
    run_ice(ref.dir)
    case = build(WORK / 'warn')
    rc = not_taken(case, ref, 'warn', edited(case, 'stop'),
                   'is not taken: ' + SOL_ERROR[len('[ERROR] '):])
    case = build(WORK / 'stop', ini_error='stop')
    rc += stop(case, edited(case, 'warn'))
    case = build(WORK / 'empty')
    rc += not_taken(case, ref, 'empty', '; input.ini emptied during the run\n',
                    'is not taken: input.ini holds no key')
    sys.exit(1 if rc else 0)
