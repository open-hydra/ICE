#!/usr/bin/env python3
"""Thermal/IG: a pressurised blob expanding from rest in a closed box at a uniform T = 1. With the
pressure work 2.5 P u_n in the energy flux the thermal energy rho c_s T moves with the mass, so T
must stay uniform: PASS iff max|T - 1| < 1e-9 over the cells with rho > 1e-3. sum(E V) is printed
as a sanity number only, not a gate: the scheme is conservative with any flux, so the sum moves
only through floors and boundaries. Exit 0 on pass, 1 on fail."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_case import read_field, NX, NY, LX, LY, H, T0

HERE = os.path.dirname(os.path.abspath(__file__))
NAME, CS, T_END, TOL = 'Thermal/IG', 1.0, 0.1, 1e-9

def fail(why):
    print('%s  -->  FAIL (%s)' % (NAME, why)); sys.exit(1)

# A table would replace specific-heat = 1 and shrink the signature under the RED threshold
if os.path.exists(os.path.join(HERE, 'INPUT', 'part-properties.dat')):
    fail('INPUT/part-properties.dat present: a table replaces specific-heat = 1 and hides the defect')
try:
    I, J, K, t, v = read_field(os.path.join(HERE, 'OUTPUT', 'part-field.tec'))
    _, _, _, _, v0 = read_field(os.path.join(HERE, 'INPUT', 'part-ic.tec'))
except (OSError, StopIteration, ValueError) as e:
    fail('no readable OUTPUT/part-field.tec: %s' % e)
if len(v) != 7 or len(v0) != 7:
    fail('expected the 7 IG variables, read %d and %d' % (len(v), len(v0)))
if not all(np.isfinite(a).all() for a in v):
    fail('non-finite values in the solution')
if t is None or t < T_END*(1.0 - 1e-12):
    fail('the run stopped at t = %s, before time-threshold = %g' % (t, T_END))

def energy(q):                                   # E = rho c_s T + (rho |u|^2 + 3 P)/2
    return q[0]*CS*q[5] + 0.5*(q[0]*(q[1]**2 + q[2]**2 + q[3]**2) + 3.0*q[4])

V = (LX/NX)*(LY/NY)*H
rho, T = v[0], v[5]
live = rho > 1e-3
dT = float(np.abs(T[live] - T0).max())
umax = float(np.sqrt(v[1]**2 + v[2]**2 + v[3]**2).max())
E0, E1 = float(energy(v0).sum()*V), float(energy(v).sum()*V)
print('%s: t = %.6e, %d of %d cells with rho > 1e-3, max|u| = %.3e' % (NAME, t, live.sum(), rho.size, umax))
print('   sum(E V) = %.15e at t = 0, %.15e at t (relative change %.2e; sanity only)' % (E0, E1, abs(E1 - E0)/E0))
print('   max|T - 1| = %.3e (tolerance %.0e)' % (dT, TOL))
if umax < 0.1:
    fail('the blob did not expand, max|u| = %.2e: nothing was tested' % umax)
if not dT < TOL:
    fail('max|T - 1| = %.3e: the thermal energy is not advected with the mass' % dT)
print('%s  -->  PASS (max|T - 1| = %.2e over %d cells)' % (NAME, dT, live.sum()))
