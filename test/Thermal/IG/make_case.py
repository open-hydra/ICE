#!/usr/bin/env python3
import numpy as np
def write_tec(path, names, X, Y, Z, cellvars):
    """ICE/ORION BLOCK tecplot: x,y,z nodal (i fastest), then cell vars, one value per line. X: (I,J,K) nodes."""
    I, J, K = X.shape; nv = 3 + len(cellvars)
    with open(path, 'w') as f:
        f.write(' VARIABLES =' + ' '.join('"%s"' % n for n in names) + '\n')
        f.write(' ZONE  T = B1-CD, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED)\n' % (I, J, K, nv))
        for a in (X, Y, Z): f.write('\n'.join('%.15E' % v for v in a.ravel(order='F')) + '\n')
        for c in cellvars:  f.write('\n'.join('%.15E' % v for v in c.ravel(order='F')) + '\n')   # c: (I-1,J-1,K-1)
def write_bc(path, ni, nj, nk, codes, payload=lambda f, i, j, k: ''):
    """ATLAS 6-column bc.txt: 'b i j k f code' per boundary cell, payload line after 401/402/403 records."""
    with open(path, 'w') as f:
        def rec(i, j, k, face):
            f.write('%8d%8d%8d%8d%8d%8d\n' % (1, i, j, k, face, codes[face]))
            if codes[face] in (401, 402, 403): f.write('   ' + payload(face, i, j, k) + '\n')
        for k in range(1, nk+1):
            for j in range(1, nj+1): rec(1, j, k, 1); rec(ni, j, k, 2)
        for k in range(1, nk+1):
            for i in range(1, ni+1): rec(i, 1, k, 3); rec(i, nj, k, 4)
        for j in range(1, nj+1):
            for i in range(1, ni+1): rec(i, j, 1, 5); rec(i, j, nk, 6)
def slab(nx, ny, Lx, Ly, h=1e-3):            # planar slab, K=2, z = ±h/2 (Doisneau layout)
    x = np.linspace(0, Lx, nx+1); y = np.linspace(0, Ly, ny+1); z = np.array([-h/2, h/2])
    return np.meshgrid(x, y, z, indexing='ij')
def wedge(nx, nr, Lx, R, deg=1.0, r0=0.0):   # axisymmetric wedge about x, K=2; r0=1e-8 puts the axis row where ATLAS puts it
    x = np.linspace(0, Lx, nx+1); r = np.linspace(r0, R, nr+1); th = np.radians([-deg/2, deg/2])
    X, Rr, Th = np.meshgrid(x, r, th, indexing='ij'); return X, Rr*np.cos(Th), Rr*np.sin(Th)
def cc(X, Y, Z):                              # cell centres from the 8 nodes
    m = lambda A: 0.125*(A[:-1,:-1,:-1]+A[1:,:-1,:-1]+A[:-1,1:,:-1]+A[1:,1:,:-1]+A[:-1,:-1,1:]+A[1:,:-1,1:]+A[:-1,1:,1:]+A[1:,1:,1:])
    return m(X), m(Y), m(Z)
def read_field(path):
    """returns (I,J,K, solutiontime, [cell arrays (I-1,J-1,K-1) per variable]) of an ICE part-field.tec"""
    import re
    L = open(path).read().split('\n'); zi = next(i for i, l in enumerate(L) if 'ZONE' in l)
    g = lambda k: int(re.search(r'\b%s\s*=\s*(\d+)' % k, L[zi]).group(1)); I, J, K = g('I'), g('J'), g('K')
    st = re.search(r'SOLUTIONTIME\s*=\s*([-+.\dE]+)', L[zi]); st = float(st.group(1)) if st else None
    v = np.array([float(t) for l in L[zi+1:] for t in l.split()]); Nn = I*J*K; Nc = (I-1)*(J-1)*max(K-1, 1)
    nv = (len(v) - 3*Nn)//Nc
    return I, J, K, st, [v[3*Nn+q*Nc:3*Nn+(q+1)*Nc].reshape((I-1, J-1, max(K-1, 1)), order='F') for q in range(nv)]

# Thermal/IG: a pressurised Gaussian blob expanding from rest in a closed box at a uniform
# temperature T = 1. specific-heat = 1 makes the thermal energy rho c_s T = rho comparable with
# the pressure energy 1.5 P = 1.5 rho, so a missing pressure work in the energy flux moves T by
# O(theta/(c_s T)) = O(1); with the correct flux T stays uniform to round-off. INPUT/ must hold
# NO part-properties.dat: a table replaces specific-heat and would shrink the signature below
# the gate's own threshold (verify.py refuses to judge if one is present).
NX, NY, LX, LY, H = 100, 100, 1.0, 1.0, 0.01     # H >= dx: a thinner slab would bind dt in z
SIGMA, THETA, T0, RP, RHO_M = 0.05, 1.0, 1.0, 1e-6, 2000.0
if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    X, Y, Z = slab(nx=NX, ny=NY, Lx=LX, Ly=LY, h=H)
    xc, yc, zc = cc(X, Y, Z)
    rho = 0.01 + np.exp(-((xc - 0.5)**2 + (yc - 0.5)**2)/(2*SIGMA**2))
    n = rho/(4.0/3.0*np.pi*RP**3*RHO_M)                # consistent with RP and density = 2000
    zero = np.zeros_like(rho)
    write_tec('INPUT/part-ic.tec', ['x','y','z','rp1','up1','vp1','wp1','Pp1','Tp1','np1'], X, Y, Z,
              [rho, zero, zero, zero, rho*THETA, T0 + zero, n])
    write_bc('INPUT/part-bc.txt', NX, NY, 1, {1: 300, 2: 300, 3: 300, 4: 300, 5: 0, 6: 0})
    with open('INPUT/part-phase.txt', 'w') as f:
        f.write('condensed-dispersed phase\nA 1\n')
    if os.path.exists('INPUT/part-properties.dat'):
        os.remove('INPUT/part-properties.dat')
    write_tec('MESH/mesh.tec', ['x','y','z'], X, Y, Z, [])
