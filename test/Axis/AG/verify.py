#!/usr/bin/env python3
"""Axis/AG: a uniform AG cloud at rest (rho = 1, P = I) in a closed one-degree wedge about x. The
wedge side faces (code 200) must carry the pressure that balances the difference between the outer
and the inner radial face of each cell; without it every cell feels a net force P*delta*dx*dr
toward the axis. A cloud at rest stays at rest, and discretely to roundoff: a closed cell has
sum(A n) = 0 and a uniform state reconstructs to itself. PASS needs a readable, finite output of
the full run and, over every cell, max(|u|,|v|,|w|)/c, max|rho - 1| and max|Pii - 1| (P11, P22,
P33) below 1e-10, with c = sqrt((P11 + P22 + P33)/rho) the closure's sound speed. Exit 0 on pass,
1 on fail."""
import os, re, sys
import numpy as np
np.seterr(all='ignore')                                  # a broken run may carry P < 0 or rho = 0
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_case import read_field

NAME = 'Axis AG'
HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-10
PCOL = {'P11': 4, 'P22': 7, 'P33': 9}                    # diagonal pressure columns of the output (0-based)

def fail(why):
    print('%s  -->  FAIL (%s)' % (NAME, why)); sys.exit(1)

def numbers(path):
    L = open(path).read().split('\n'); zi = next(i for i, l in enumerate(L) if 'ZONE' in l)
    return np.array([float(t) for l in L[zi+1:] for t in l.split()])

out = os.path.join(HERE, 'OUTPUT', 'part-field.tec')
try:
    I, J, K, st, v = read_field(out)
    allv = numbers(out)
except (OSError, StopIteration, ValueError, AttributeError) as e:
    fail('no readable OUTPUT/part-field.tec: %s' % e)
try:
    log = open(os.path.join(HERE, 'logfile')).read()
except OSError:
    log = ''
niter = int(re.search(r'^\s*iter-threshold\s*=\s*(\d+)', open(os.path.join(HERE, 'input.ini')).read(), re.M).group(1))

bad = []
if not np.isfinite(allv).all():
    bad.append('%d non-finite values' % int((~np.isfinite(allv)).sum()))
if 'Time of operation' not in log:
    bad.append('no "Time of operation" in logfile')
if st is None or abs(st) != niter:
    bad.append('SOLUTIONTIME %s, not the %d iterations of the run' % (st, niter))
rho, vel = v[0], np.abs(np.stack(v[1:4])).max(axis=0)
P = [v[q] for q in PCOL.values()]
c = np.sqrt(3.0*sum(P)/len(P)/rho)
mach = vel/c
im = np.unravel_index(np.nanargmax(mach), mach.shape) if np.isfinite(mach).any() else (0, 0, 0)
m, drho = float(np.max(mach)), float(np.abs(rho - 1.0).max())
dP = max(float(np.abs(p - 1.0).max()) for p in P)
if not m < TOL:
    bad.append('max(|u|,|v|,|w|)/c %.3e at i=%d j=%d' % (m, im[0] + 1, im[1] + 1))
if not drho < TOL:
    bad.append('max|rho - 1| %.3e' % drho)
if not dP < TOL:
    bad.append('max|%s - 1| %.3e' % ('|'.join(PCOL), dP))
if bad:
    fail('; '.join(bad))
print('%s  -->  PASS (%d cells: max(|u|,|v|,|w|)/c %.3e, max|rho - 1| %.3e, max|%s - 1| %.3e)'
      % (NAME, rho.size, m, drho, '|'.join(PCOL), dP))
