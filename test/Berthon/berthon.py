"""Comparison of a Berthon Riemann problem against its exact solution.

The three cases are one-dimensional Riemann problems for the anisotropic
Gaussian closure. Their exact solutions are not ICE's own output: they are the
analytical wave patterns tabulated in `Results/<case>/<case>_exact.tec` as
piecewise-linear profiles of

    density, u, v, P11, P22, and det(P) = P11 P22 - P12^2,

the last being the quantity the Gaussian closure transports along the contact.

The solutions contain shocks and contacts, so a pointwise comparison would only
measure how many cells the scheme smears them over. The check is therefore the
L1 error over the domain, normalised by the variation of the exact profile, at
a tolerance that a monotone second-order scheme on this mesh meets with room to
spare but a wrong wave pattern does not.
"""
import math
import re
import sys
from pathlib import Path

#: Exact-solution zone name -> index of the matching ICE primitive variable, or
#: a callable evaluated on the twelve AG primitives of one cell.
QUANTITIES = (
    ('Density',      'rho',    lambda q: q[0]),
    ('X Velocity',   'u',      lambda q: q[1]),
    ('Y Velocity',   'v',      lambda q: q[2]),
    ('XX Pressure',  'P11',    lambda q: q[4]),
    ('YY Pressure',  'P22',    lambda q: q[7]),
    ('Det Pressure', 'det(P)', lambda q: q[4]*q[7] - q[5]*q[5]),
)


def read_exact(path):
    """Read the piecewise-linear exact profiles, one Tecplot zone each."""
    lines = Path(path).read_text().splitlines()
    zones, i = {}, 0
    while i < len(lines):
        m = re.search(r'ZONE\s+T\s*=\s*"([^"]+)"', lines[i])
        if not m:
            i += 1
            continue
        n = int(re.search(r'\bI\s*=\s*(\d+)', lines[i + 1]).group(1))
        i += 3                                   # past the I=/DATAPACKING lines
        vals = []
        while len(vals) < 2 * n:
            vals += [float(t) for t in lines[i].split()]
            i += 1
        zones[m.group(1)] = (vals[:n], vals[n:])
    return zones


def read_field(path):
    """Read ICE's cell-centred AG solution along x, on a 1D (nx,1,1) block."""
    lines = Path(path).read_text().splitlines()
    head = next(k for k, l in enumerate(lines) if 'ZONE' in l)
    I, J, K = (int(re.search(r'\b%s\s*=\s*(\d+)' % a, lines[head]).group(1))
               for a in 'IJK')
    nums = []
    for line in lines[head + 1:]:
        nums += [float(t) for t in line.split()]
    nnode = I * J * K
    ncell = (I - 1) * (J - 1) * max(1, K - 1)
    x = nums[:nnode]
    var = [nums[3*nnode + v*ncell: 3*nnode + (v+1)*ncell] for v in range(12)]
    xc = [0.5 * (x[i] + x[i + 1]) for i in range(I - 1)]
    return xc, var


def sample(xe, ve, x):
    """Value of a piecewise-linear profile at x; repeated abscissae are jumps."""
    if x <= xe[0]:
        return ve[0]
    if x >= xe[-1]:
        return ve[-1]
    for i in range(len(xe) - 1):
        if xe[i] <= x <= xe[i + 1]:
            dx = xe[i + 1] - xe[i]
            if dx <= 0.0:
                continue
            return ve[i] + (ve[i + 1] - ve[i]) * (x - xe[i]) / dx
    return ve[-1]


def verify(case, tol):
    """Compare OUTPUT/part-field.tec against the stored exact solution."""
    here = Path(__file__).resolve().parent
    exact = read_exact(here / 'Results' / case / ('%s_exact.tec' % case))
    xc, var = read_field(Path('OUTPUT/part-field.tec'))

    green, red, reset = '\033[92m', '\033[91m', '\033[0m'
    print('Berthon %s  (%d cells)' % (case, len(xc)))
    print('   %-14s %-12s %s' % ('quantity', 'L1 error', 'tolerance'))

    worst = 0.0
    for zone, label, get in QUANTITIES:
        xe, ve = exact[zone]
        ref = [sample(xe, ve, x) for x in xc]
        num = [get([var[v][c] for v in range(12)]) for c in range(len(xc))]
        span = max(ve) - min(ve)
        err = sum(abs(a - b) for a, b in zip(num, ref)) / (len(xc) * span)
        worst = max(worst, err)
        print('   %-14s %-12.3e %.3e  %s' %
              (label, err, tol, 'ok' if err <= tol else 'FAILED'))

    passed = worst <= tol
    print('Berthon %s  -->  %s' % (case, (green + 'PASS' + reset) if passed
                                   else (red + 'FAIL' + reset)))
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(verify(sys.argv[1], float(sys.argv[2])))
