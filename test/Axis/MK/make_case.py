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

# One-degree wedge about the x axis, 40 x 20 x 1 cells, the axis row at r = 1e-8 as ATLAS
# places it: face 3 is a regular face of tiny area and carries 300, the wedge side faces
# 5/6 carry 200. A pressureless cloud moves along x at 20 m/s with a density that depends
# on r only, fed by a 402 inlet with the same profile: the field is stationary at any order.
# The Jameson detector reads the side-face ghosts of every cell, which only a 200 branch fills.
if __name__ == '__main__':
    import math, os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    R, Rp, rho_m, u0, T0 = 0.05, 1e-5, 2000.0, 20.0, 300.0
    X, Y, Z = wedge(nx=40, nr=20, Lx=0.1, R=R, deg=1.0, r0=1e-8)
    shp = tuple(s - 1 for s in X.shape)
    one = np.ones(shp)
    xc, yc, zc = cc(X, Y, Z)
    rho = 1.0 + 4.0*(np.sqrt(yc**2 + zc**2)/R)**2
    n = rho / (4.0/3.0*math.pi*Rp**3*rho_m)              # consistent with Rp: the solver recomputes Rp from rho and n
    write_tec('INPUT/part-ic.tec', ['x','y','z','rho','u','v','w','T','n'], X, Y, Z,
              [rho, u0*one, 0*one, 0*one, T0*one, n])
    # 402 payload: mass flux rho*u of the row, speed, direction along the face normal, T, radius.
    # Full precision: the inlet ghost density is massflux/speed, and 7 digits would move the rows by 1e-7
    inlet = lambda f, i, j, k: '%.16E %.16E normal, normal, %.16E %.16E' % (rho[0, j-1, k-1]*u0, u0, T0, Rp)
    write_bc('INPUT/part-bc.txt', shp[0], shp[1], shp[2],
             {1: 402, 2: 400, 3: 300, 4: 300, 5: 200, 6: 200}, inlet)
    with open('INPUT/part-phase.txt', 'w') as f:
        f.write('condensed-dispersed phase\nA 1\n')
    write_tec('MESH/mesh.tec', ['x','y','z'], X, Y, Z, [])
