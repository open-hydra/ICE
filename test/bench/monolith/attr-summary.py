#!/usr/bin/env python3
"""Summary and verdict of an ice-attr.slurm job (plan openmp-loose-ends, L1).

    attr-summary.py <ice-attr-<id>.out> [commit order ...]

Per binary and case set (R24 = the 24-block cut, R1 = the single block, both 1 x 80,
one thread group): the median and the MAD (median absolute deviation) of steady_iter
over the passes, and the medians of the ICE Detail regions of each run's last window,
grouped as the commits touched them -- per-cell (dt copy zero source resid update
diag; e69ebcc), ghost/bound (ghost1 ghost2 bound; fddd1cb), flux (055f6d4).

The region times come from the IOK line, i.e. each run's LAST window only: one sample per run.
For a verdict on one region, read the run logs instead and average the steady half of the
`ICE Detail` windows (job 304520: a +2.1 ms last-window rise of compute_source was noise there).

Verdict per consecutive pair of commits: a step is REAL when the gap of the medians
exceeds the larger of 3 x the larger MAD of the two and 1 % of the earlier median,
otherwise NOISE. Python 3.6 (monolith's).
"""
import re
import sys
from collections import OrderedDict, defaultdict

GROUPS = OrderedDict([
    ('percell', ('dt', 'copy', 'zero', 'source', 'resid', 'update', 'diag')),
    ('ghostbound', ('ghost1', 'ghost2', 'bound')),
    ('flux', ('flux',)),
])


def median(v):
    v = sorted(v)
    n = len(v)
    if n == 0:
        return float('nan')
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def mad(v):
    m = median(v)
    return median([abs(x - m) for x in v])


def detail(s):
    p = [x for x in s.split(';') if x]
    return {p[i]: float(p[i + 1]) for i in range(0, len(p) - 1, 2)}


def main():
    path = sys.argv[1]
    order = sys.argv[2:]
    steady = defaultdict(list)            # (commit, set) -> [s]
    phases = defaultdict(lambda: defaultdict(list))
    fails = []
    for line in open(path, errors='replace'):
        if line.startswith('IFAIL'):
            fails.append(line.strip())
        if not line.startswith('IOK '):
            continue
        f = dict(re.findall(r'(\w+)=(\S+)', line))
        m = re.match(r'a-(\w+)-(R24|R1)-p(\d+)$', f['tag'])
        if not m:
            continue
        c, cs = m.group(1), m.group(2)
        if c not in order:
            order.append(c)
        steady[(c, cs)].append(float(f['steady_iter']))
        if 'detail' in f:
            d = detail(f['detail'])
            for g, keys in GROUPS.items():
                phases[(c, cs)][g].append(sum(d.get(k, 0.0) for k in keys))
    for x in fails:
        print('FAILED RUN', x)
    for cs in ('R24', 'R1'):
        print('\n## %s, 1 x 80, ICE_THREAD_GROUPS=1' % cs)
        print('| commit | n | steady median s | MAD s | min | max | per-cell | ghost/bound | flux |')
        print('|---|---|---|---|---|---|---|---|---|')
        for c in order:
            v = steady.get((c, cs), [])
            if not v:
                continue
            ph = phases[(c, cs)]
            print('| %s | %d | %.4f | %.4f | %.4f | %.4f | %s | %s | %s |' % (
                c, len(v), median(v), mad(v), min(v), max(v),
                *['%.4f' % median(ph[g]) if ph[g] else '-' for g in GROUPS]))
        print('\n| step | gap of medians | threshold | verdict | per-cell gap | ghost/bound gap | flux gap |')
        print('|---|---|---|---|---|---|---|')
        prev = None
        for c in order:
            if (c, cs) not in steady:
                continue
            if prev is not None:
                a, b = steady[(prev, cs)], steady[(c, cs)]
                gap = median(b) - median(a)
                thr = max(3.0 * max(mad(a), mad(b)), 0.01 * median(a))
                pg = []
                for g in GROUPS:
                    pa, pb = phases[(prev, cs)][g], phases[(c, cs)][g]
                    pg.append('%+.4f' % (median(pb) - median(pa)) if pa and pb else '-')
                print('| %s -> %s | %+.4f (%+.1f %%) | %.4f | %s | %s | %s | %s |' % (
                    prev, c, gap, 100.0 * gap / median(a), thr,
                    'REAL' if abs(gap) > thr else 'noise', *pg))
            prev = c
        have = [c for c in order if (c, cs) in steady]
        if len(have) > 2:
            a, b = steady[(have[0], cs)], steady[(have[-1], cs)]
            gap = median(b) - median(a)
            thr = max(3.0 * max(mad(a), mad(b)), 0.01 * median(a))
            print('| %s -> %s (whole) | %+.4f (%+.1f %%) | %.4f | %s | | | |' % (
                have[0], have[-1], gap, 100.0 * gap / median(a), thr, 'REAL' if abs(gap) > thr else 'noise'))


if __name__ == '__main__':
    main()
