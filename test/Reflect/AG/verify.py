#!/usr/bin/env python3
"""Reflect/AG: a sheared AG blob moving along a code-200 symmetry plane. The mirrored tensor H P H
makes the plane's x-momentum flux P12 vanish, and the background's fluxes through the 400 faces
cancel, so the total x-momentum sum(rho u V) must be conserved: PASS iff its relative change
from the initial condition is < 1e-9. Exit 0 on pass, 1 on fail."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_case import read_field, NX, NY, LX, LY, H

HERE = os.path.dirname(os.path.abspath(__file__))
NAME, T_END, TOL = 'Reflect/AG', 0.1, 1e-9

def fail(why):
    print('%s  -->  FAIL (%s)' % (NAME, why)); sys.exit(1)

try:
    I, J, K, t, v = read_field(os.path.join(HERE, 'OUTPUT', 'part-field.tec'))
    _, _, _, _, v0 = read_field(os.path.join(HERE, 'INPUT', 'part-ic.tec'))
except (OSError, StopIteration, ValueError) as e:
    fail('no readable OUTPUT/part-field.tec: %s' % e)
if len(v) != 12 or len(v0) != 12:
    fail('expected the 12 AG variables, read %d and %d' % (len(v), len(v0)))
if not all(np.isfinite(a).all() for a in v):
    fail('non-finite values in the solution')
if t is None or t < T_END*(1.0 - 1e-12):
    fail('the run stopped at t = %s, before time-threshold = %g' % (t, T_END))

V = (LX/NX)*(LY/NY)*H
M0, M1 = float((v0[0]*v0[1]).sum()*V), float((v[0]*v[1]).sum()*V)
drift = abs(M1 - M0)/M0
print('%s: t = %.6e, M_x = %.15e at t = 0, %.15e at t, relative change %.3e (tolerance %.0e)'
      % (NAME, t, M0, M1, drift, TOL))
if not drift < TOL:
    fail('x-momentum drifts by %.3e: the symmetry plane exerts a shear' % drift)
print('%s  -->  PASS (x-momentum conserved to %.2e)' % (NAME, drift))
