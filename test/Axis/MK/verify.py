#!/usr/bin/env python3
"""Axis/MK: a one-degree wedge about x whose side faces carry code 200 and whose axis face
carries 300. A pressureless cloud moving along x with a density that depends on r only is
stationary at any order, so the axis row (j = 1) must keep its initial density: a regression
pin, not the discriminator (that is a build trapping invalid operations with a signalling-NaN
ghost sentinel, in which the Jameson detector reads the side-face ghosts). PASS also needs a
readable, finite output of the full run and the boundary census in the log. Exit 0 on pass, 1
on fail."""
import os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_case import read_field

NAME = 'Axis MK'
HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-9
CENSUS = [('Symmetry', 80), ('Axisymmetry', 1600)]      # faces 3+4: 2 x 40; faces 5+6: 2 x 40 x 20

def fail(why):
    print('%s  -->  FAIL (%s)' % (NAME, why)); sys.exit(1)

def numbers(path):
    L = open(path).read().split('\n'); zi = next(i for i, l in enumerate(L) if 'ZONE' in l)
    return np.array([float(t) for l in L[zi+1:] for t in l.split()])

out = os.path.join(HERE, 'OUTPUT', 'part-field.tec')
try:
    I, J, K, st, v = read_field(out)
    ic = read_field(os.path.join(HERE, 'INPUT', 'part-ic.tec'))[4]
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
for name, count in CENSUS:
    if not re.search(r'^\s*%s\s+%d\s*$' % (name, count), log, re.M):
        bad.append('census line "%s %d" missing' % (name, count))
d_axis = float(np.abs(v[0][:, 0, :] - ic[0][:, 0, :]).max())
d_all = float(np.abs(v[0] - ic[0]).max())
if not d_axis < TOL:
    bad.append('axis row max|rho - rho_ic| %.3e >= %.0e' % (d_axis, TOL))
if bad:
    fail('; '.join(bad))
print('%s  -->  PASS (axis row max|rho - rho_ic| %.3e, whole field %.3e; census %s)'
      % (NAME, d_axis, d_all, ', '.join('%s %d' % c for c in CENSUS)))
