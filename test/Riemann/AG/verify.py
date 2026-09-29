#!/usr/bin/env python3
"""Riemann/AG: an anisotropic Sod tube (P11/P22 = 1e7) through explicit Euler and first order, against
the exact gamma = 3 Riemann solution the 1-D AG system reduces to (Toro, Riemann Solvers, ch. 4).
PASS iff every value is finite, min(P11) > 1e3 (no floor hits: the floor is 1e-25, the states
1e4..1e5), the L1 density error sum|rho - rho_ex|/sum(rho_ex) <= 0.05 (first-order Rusanov on 400
cells: about 0.02), max|P22/rho - 0.01| < 1e-9 (P22/rho is passively advected), and the total
variation of rho and of P11 is at most 1.05 times the exact one. The last check is what sees an
under-estimated Rusanov speed: the global step shrinks as the star region speeds up, so the
under-dissipated scheme does not blow up here, it oscillates (total variation 1.6x and 3.6x the
exact one) while its L1 error stays as small as the correct scheme's. Exit 0 on pass, 1 on fail."""
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_case import read_field, NX, LX, X0, RHO_L, P_L, RHO_R, P_R, PT_RATIO

HERE = os.path.dirname(os.path.abspath(__file__))
NAME, GAMMA, T_END, L1_TOL, P11_MIN, PT_TOL, TV_TOL = 'Riemann/AG', 3.0, 2.7e-4, 0.05, 1e3, 1e-9, 1.05

def fail(why):
    print('%s  -->  FAIL (%s)' % (NAME, why)); sys.exit(1)

def exact_riemann(WL, WR, gam):
    """Exact Riemann solution of the Euler equations, W = (rho, u, p): Newton on
    f_L(p) + f_R(p) + u_R - u_L = 0 from p0 = (p_L + p_R)/2, then Toro's sampling.
    Returns (p*, u*, sample) with sample(S) = (rho, u, p) at S = (x - x0)/t."""
    rL, uL, pL = WL; rR, uR, pR = WR
    cL, cR = math.sqrt(gam*pL/rL), math.sqrt(gam*pR/rR)
    g1, g2, gm = (gam - 1.0)/(2.0*gam), (gam + 1.0)/(2.0*gam), (gam - 1.0)/(gam + 1.0)
    def fK(p, rK, pK, cK):                       # f_K and its derivative
        if p > pK:                               # shock
            A, B = 2.0/((gam + 1.0)*rK), gm*pK
            q = math.sqrt(A/(p + B))
            return (p - pK)*q, q*(1.0 - 0.5*(p - pK)/(B + p))
        return 2.0*cK/(gam - 1.0)*((p/pK)**g1 - 1.0), (p/pK)**(-g2)/(rK*cK)   # rarefaction
    p = 0.5*(pL + pR)
    for _ in range(100):
        fl, dl = fK(p, rL, pL, cL); fr, dr = fK(p, rR, pR, cR)
        pn = max(p - (fl + fr + uR - uL)/(dl + dr), 1e-12*p)
        done = abs(pn - p) <= 1e-15*0.5*(pn + p)
        p = pn
        if done:
            break
    ps = p
    us = 0.5*(uL + uR) + 0.5*(fK(ps, rR, pR, cR)[0] - fK(ps, rL, pL, cL)[0])
    def sample(S):
        if S <= us:                              # left of the contact
            if ps > pL:
                if S <= uL - cL*math.sqrt(g2*ps/pL + g1): return rL, uL, pL
                return rL*(ps/pL + gm)/(gm*ps/pL + 1.0), us, ps
            if S <= uL - cL: return rL, uL, pL
            if S > us - cL*(ps/pL)**g1: return rL*(ps/pL)**(1.0/gam), us, ps
            b = 2.0/(gam + 1.0) + gm/cL*(uL - S)
            return rL*b**(2.0/(gam - 1.0)), 2.0/(gam + 1.0)*(cL + 0.5*(gam - 1.0)*uL + S), pL*b**(2.0*gam/(gam - 1.0))
        if ps > pR:                              # right of the contact
            if S >= uR + cR*math.sqrt(g2*ps/pR + g1): return rR, uR, pR
            return rR*(ps/pR + gm)/(gm*ps/pR + 1.0), us, ps
        if S >= uR + cR: return rR, uR, pR
        if S <= us + cR*(ps/pR)**g1: return rR*(ps/pR)**(1.0/gam), us, ps
        b = 2.0/(gam + 1.0) - gm/cR*(uR - S)
        return rR*b**(2.0/(gam - 1.0)), 2.0/(gam + 1.0)*(-cR + 0.5*(gam - 1.0)*uR + S), pR*b**(2.0*gam/(gam - 1.0))
    return ps, us, sample

# The oracle first: Toro's Sod test (Table 4.3, gamma = 1.4) has p* = 0.30313, u* = 0.92745
ps, us, _ = exact_riemann((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 1.4)
if abs(ps - 0.30313) > 1e-5 or abs(us - 0.92745) > 1e-5:
    fail('exact solver self-test: p* = %.6f, u* = %.6f against Toro 0.30313, 0.92745' % (ps, us))

try:
    I, J, K, t, v = read_field(os.path.join(HERE, 'OUTPUT', 'part-field.tec'))
except (OSError, StopIteration, ValueError) as e:
    fail('no readable OUTPUT/part-field.tec: %s' % e)
if len(v) != 12:
    fail('expected the 12 AG variables, read %d' % len(v))
if not all(np.isfinite(a).all() for a in v):
    fail('non-finite values in the solution')
if t is None or t < T_END*(1.0 - 1e-12):
    fail('the run stopped at t = %s, before time-threshold = %g' % (t, T_END))

rho, P11, P22 = v[0][:, :, 0], v[4][:, :, 0], v[7][:, :, 0]
ps, us, sample = exact_riemann((RHO_L, 0.0, P_L), (RHO_R, 0.0, P_R), GAMMA)
xc = (np.arange(NX) + 0.5)*LX/NX
ex = np.array([sample((x - X0)/t) for x in xc])
rho_ex, p_ex = ex[:, 0], ex[:, 2]
l1 = float(np.abs(rho.mean(axis=1) - rho_ex).sum()/rho_ex.sum())
tv = lambda a: float(np.abs(np.diff(a)).sum())
tv_rho = tv(rho.mean(axis=1))/tv(rho_ex)
tv_p11 = tv(P11.mean(axis=1))/tv(p_ex)
p11min = float(P11.min())
dpt = float(np.abs(P22/rho - PT_RATIO).max())
print('%s: t = %.6e, exact gamma = 3 solution p* = %.6e, u* = %.6e' % (NAME, t, ps, us))
print('   min(P11) = %.6e (> %.0e), L1(rho) = %.4e (<= %.2f), max|P22/rho - %.2f| = %.3e (< %.0e)'
      % (p11min, P11_MIN, l1, L1_TOL, PT_RATIO, dpt, PT_TOL))
print('   total variation / exact: rho %.4f, P11 %.4f (<= %.2f)' % (tv_rho, tv_p11, TV_TOL))
bad = []
if not p11min > P11_MIN:
    bad.append('min(P11) = %.3e: floor hits' % p11min)
if not l1 <= L1_TOL:
    bad.append('L1(rho) = %.3e' % l1)
if not dpt < PT_TOL:
    bad.append('max|P22/rho - 0.01| = %.3e' % dpt)
if not (tv_rho <= TV_TOL and tv_p11 <= TV_TOL):
    bad.append('total variation %.3f (rho), %.3f (P11) times the exact: oscillations' % (tv_rho, tv_p11))
if bad:
    fail(', '.join(bad))
print('%s  -->  PASS (L1(rho) = %.3e against the exact gamma = 3 solution)' % (NAME, l1))
