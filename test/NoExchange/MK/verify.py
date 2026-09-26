#!/usr/bin/env python3
"""NoExchange/MK: with drag = NoDrag, heat = NoHeat and emiss = 0 the particles must keep their
initial state EXACTLY (at rest, 300 K, unchanged density) although the gas moves at 10 m/s and is
100 K hotter. Exit 0 on pass, 1 on fail."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_case import read_field

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    I, J, K, st, v = read_field(os.path.join(HERE, 'OUTPUT', 'part-field.tec'))
except (OSError, StopIteration, ValueError) as e:
    print('NoExchange/MK  -->  FAIL (no readable OUTPUT/part-field.tec: %s)' % e); sys.exit(1)
rho, u, vv, w, T = v[0], v[1], v[2], v[3], v[4]
checks = [('rho', rho, 1e-3), ('u', u, 0.0), ('v', vv, 0.0), ('w', w, 0.0), ('T', T, 300.0)]
bad = [(n, float(abs(a - x).max())) for n, a, x in checks if (a != x).any()]
if bad:
    print('NoExchange/MK  -->  FAIL ' + ', '.join('%s max|dev| %.3e' % b for b in bad)); sys.exit(1)
print('NoExchange/MK  -->  PASS (%d cells: rho, u, v, w, T bit-for-bit at their initial values)' % rho.size)
