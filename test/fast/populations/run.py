"""Fast test: a multi-population BC file must reach the right particle family.

ATLAS writes one copy of the whole boundary table per (material, population) pair
of the dispersed phase, block by block, each copy carrying that population's
injection data. ICE binds copy c to family c. Getting that wrong has no symptom
other than the wrong size class being injected, so it is pinned three ways:

  proportional  three families, three copies whose mass fluxes are 1, 2 and 4
                times a reference, and an initial loading scaled the same way.
                The equations are linear in the loading at fixed velocity and
                radius, so family p must come out as 2^(p-1) times family 1, to
                round-off, with its velocity and temperature fields identical to
                family 1's -- those do not carry the loading at all.
  shared        one copy and three families: the single record speaks for all of
                them, which is what every case written before populations existed
                relies on. The three families must stay identical.
  mismatch      two copies and three families: neither shared nor one-per-family,
                so the setup must stop instead of guessing.
  residuals     the residual history of the first two runs: after the five norms
                of the domain, the same five of each family, in turn. Family p of
                the proportional run has 2^(p-1) times family 1's density residual
                and family 1's velocity and temperature ones; the shared run's
                three families have the same; and in every row each norm of the
                domain is the root of the sum of the families' squares.
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import Case, GREEN, RED, RESET, ICE_BIN                  # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
NVAR = 6                      # MK: rho, u, v, w, T, n
RHO, U, V, W, T, N = range(NVAR)

GP0, VINJ, TINJ, RINJ = 2.0, 10.0, 300.0, 1.0e-5
RHO0, N0 = 1.0e-3, 1.0e-3 / (1000.0 * 4.0 / 3.0 * 3.141592653589793 * RINJ**3)


def build(work, copies, families, scaled):
    """A duct with an inlet on face 1, `copies` copies of the table, `families` families."""
    case = Case(work, nx=16, ny=4, Lx=0.16, Ly=0.04)
    scales = [2.0**c for c in range(families)] if scaled else None
    case.particles(RHO0, 0.0, 0.0, 0.0, TINJ, N0, families=families, scales=scales)
    case.boundaries(inlets=[{'gp': GP0 * 2.0**c, 'v': VINJ, 'T': TINJ, 'r': RINJ}
                            for c in range(copies)])
    case.ini(iters=30, cfl=0.8, rk='RK2', closures=['MK'] * families)
    return case


def family(sol, p):
    """{name: values} of family p (1-based) in a multi-family solution."""
    base = (p - 1) * NVAR
    return [sol['var'][base + v] for v in range(NVAR)]


def report(name, ok, detail):
    print('  %-14s %s  %s' % (name, (GREEN + 'PASS' + RESET) if ok else (RED + 'FAIL' + RESET),
                              detail))
    return 0 if ok else 1


TOL = 1.0e-12          # one scaling and one division of round-off, no more


def relative(a, b):
    """max |a - b| over the field, as a fraction of the field's own magnitude."""
    scale = max(abs(x) for x in b) or 1.0
    return max(abs(x - y) for x, y in zip(a, b)) / scale


def proportional():
    """Copy c drives family c: family p must be 2^(p-1) times family 1."""
    sol = build(WORK / 'proportional', copies=3, families=3, scaled=True).run()
    f1 = family(sol, 1)
    worst, worst_same = 0.0, 0.0
    for p in (2, 3):
        fp = family(sol, p)
        s = 2.0 ** (p - 1)
        for v in (RHO, N):
            worst = max(worst, relative(fp[v], [s * x for x in f1[v]]))
        for v in (U, V, W, T):
            worst_same = max(worst_same, relative(fp[v], f1[v]))
    injected = max(f1[RHO])
    if injected <= RHO0:
        return report('proportional', False,
                      'nothing was injected (max rho %g, initial %g)' % (injected, RHO0))
    return report('proportional', worst <= TOL and worst_same <= TOL,
                  'loading off by %.1e, shared fields by %.1e (tolerance %.0e)'
                  % (worst, worst_same, TOL))


def shared():
    """One copy, three families: the record speaks for all of them."""
    sol = build(WORK / 'shared', copies=1, families=3, scaled=False).run()
    f1 = family(sol, 1)
    worst = 0.0
    for p in (2, 3):
        fp = family(sol, p)
        for v in range(NVAR):
            worst = max(worst, max(abs(a - b) for a, b in zip(fp[v], f1[v])))
    return report('shared', worst == 0.0, 'families differ by %g' % worst)


def mismatch():
    """Two copies, three families: the setup must stop and say why."""
    case = build(WORK / 'mismatch', copies=2, families=3, scaled=False)
    out = subprocess.run([str(ICE_BIN)], cwd=str(case.dir),
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    text = out.stdout.decode('utf-8', 'replace')
    ok = out.returncode != 0 and 'copies of the boundary table' in text
    return report('mismatch', ok, 'exit %d, %s' % (
        out.returncode, 'reported' if 'copies of the boundary table' in text else 'no message'))


NRES = 5               # the norms of a residual row: rho, u, v, w, T
TOL_ROW = 1.0e-9       # the history prints 10 significant digits (E20.10)


def history(run, families):
    """The rows of a run's residual history, each split as (domain, [family 1, ...]) after the
    iteration and the time; None, and the offending row, if a row has not 2 + 5 (1 + families) values."""
    rows = []
    for line in (WORK / run / 'OUTPUT' / 'part-residual-history.dat').read_text().split('\n'):
        if not line.strip():
            continue
        v = [float(x) for x in line.split()]
        if len(v) != 2 + NRES * (1 + families):
            return None, line
        rows.append((v[2:2 + NRES], [v[2 + NRES * p:2 + NRES * (p + 1)] for p in range(1, families + 1)]))
    return rows, None


def residuals():
    """The family columns of the residual history: their scaling, and the domain's from them."""
    worst_scale, worst_sum, nrow = 0.0, 0.0, 0
    for run, scaled in (('proportional', True), ('shared', False)):
        rows, bad = history(run, 3)
        if rows is None:
            return report('residuals', False, '%s: a row without 2 + 5 x 4 values: %r' % (run, bad[:80]))
        for dom, fam in rows:
            nrow += 1
            for c in range(NRES):
                s = [2.0**p if (scaled and c == 0) else 1.0 for p in range(3)]
                for p in (1, 2):
                    ref = s[p] * fam[0][c]
                    worst_scale = max(worst_scale, abs(fam[p][c] - ref) / (abs(ref) or 1.0))
                root = sum(f[c]**2 for f in fam)**0.5
                worst_sum = max(worst_sum, abs(dom[c] - root) / (abs(dom[c]) or 1.0))
    if nrow == 0:
        return report('residuals', False, 'no residual row written')
    return report('residuals', worst_scale <= TOL_ROW and worst_sum <= TOL_ROW,
                  '%d rows: families off their scaling by %.1e, domain off the root of the sum by %.1e '
                  '(tolerance %.0e)' % (nrow, worst_scale, worst_sum, TOL_ROW))


if __name__ == '__main__':
    print('Population binding')
    rc = proportional() + shared() + mismatch() + residuals()
    sys.exit(1 if rc else 0)
