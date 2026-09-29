"""Shared machinery for ICE's code-verification cases.

These cases differ from the ones under test/Doisneau: they compare ICE against a
solution known independently of ICE (a closed form, or an accurate integration of
the ODE the solver is supposed to integrate), not against a stored output. They
build their own input files in a scratch directory, so there is no case data to
keep in the repository and no reference to regenerate.

The physical setup is always spatially uniform (or uniformly advected), so the
answer per cell is the same everywhere and the mesh only has to be big enough to
exercise the flux loops.
"""
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ICE_BIN = Path(os.environ.get('ICE_BIN', str(ROOT / 'bin' / 'ICE')))

GREEN, RED, RESET = '\033[92m', '\033[91m', '\033[0m'

TOLL = 1.e-20          # the guard ICE adds to its denominators (Lib_Drag, Mod_Sources)


# ---------------------------------------------------------------------------
#  Case construction
# ---------------------------------------------------------------------------

def linspace(a, b, n):
    """n equally spaced values from a to b inclusive."""
    if n == 1:
        return [a]
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def _fmt(values):
    return '\n'.join(' %.15E' % v for v in values)


def _nodal(xn, yn, zn, component):
    """Expand a node coordinate over the (K, J, I) ordering Tecplot BLOCK uses."""
    out = []
    for z in zn:
        for y in yn:
            for x in xn:
                out.append((x, y, z)[component])
    return out


def write_tec(path, names, xn, yn, zn, cellvars, zone='B1'):
    """Structured single-zone Tecplot BLOCK file: nodal x,y,z then cell-centred data."""
    nvar = 3 + len(names)
    head = [' VARIABLES ="x" "y" "z" ' + ' '.join('"%s"' % v for v in names),
            ' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, '
            'VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED)' % (zone, len(xn), len(yn), len(zn), nvar)]
    body = [_fmt(_nodal(xn, yn, zn, c)) for c in (0, 1, 2)]
    body += [_fmt(v) for v in cellvars]
    Path(path).write_text('\n'.join(head + body) + '\n')


def write_bc(path, nx, ny, nz, mode='extrapolation', inlets=None):
    """ATLAS-format BC file for a single block.

    ICE reads exactly 2*(ny*nz + nx*nz + nx*ny) records per copy of the boundary
    table, one per boundary face cell, so every face has to be present - including
    the degenerate k faces, which get the null code 0. 'periodic-x' closes faces 1
    and 2 onto each other with code 201.

    `inlets` makes face 1 a code-402 inlet and turns the file into a
    multi-population one: one (gp, v, T, r) dict per copy, written in the order
    ATLAS writes its (material, population) copies, so copy c belongs to family c.
    """
    copies = inlets if inlets else [None]
    table = []

    def header(i, j, k, f, code):
        return '%8d%8d%8d%8d%8d%8d' % (1, i, j, k, f, code)

    def connection(i, j, k, f):
        # The ghost cells of face f take the cells behind the opposite face
        src, fs = (nx, 2) if f == 1 else (1, 1)
        return ['%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (1, src, j, k, fs, 1, 0, 0, 1)]

    def inlet(spec):
        # ATLAS writes nine fields; ICE consumes the first six. `normal,` in the two
        # angle columns is its way of saying "along the face normal".
        return ['%14.5E%14.5E%16s%16s%14.5E%14.5E%14.5E %s%14.5E'
                % (spec['gp'], spec['v'], 'normal,', 'normal,', spec['T'], spec['r'],
                   0.0, 'Dirac', 0.0)]

    side = 201 if mode == 'periodic-x' else 400
    for f, (i0, j0, k0) in ((1, (1, None, None)), (2, (nx, None, None))):
        for k in range(1, nz + 1):
            for j in range(1, ny + 1):
                if f == 1 and inlets:
                    table.append((header(i0, j, k, f, 402), 'inlet'))
                    continue
                table.append((header(i0, j, k, f, side), None))
                if side == 201:
                    table.append((connection(i0, j, k, f), None))
    for f, j0 in ((3, 1), (4, ny)):
        for k in range(1, nz + 1):
            for i in range(1, nx + 1):
                table.append((header(i, j0, k, f, 400), None))
    for f, k0 in ((5, 1), (6, nz)):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                table.append((header(i, j, k0, f, 0), None))

    rows = []
    for spec in copies:
        for entry, kind in table:
            rows += entry if isinstance(entry, list) else [entry]
            if kind == 'inlet':
                rows += inlet(spec)

    Path(path).write_text('\n'.join(rows) + '\n')


def write_ini(path, sections):
    """input.ini from {section: {key: value}}."""
    out = []
    for section, keys in sections.items():
        out.append('[%s]' % section)
        for key, value in keys.items():
            if isinstance(value, bool):
                value = 'true' if value else 'false'
            out.append('%-22s = %s' % (key, value))
        out.append('')
    Path(path).write_text('\n'.join(out))


class Case(object):
    """A generated ICE case in a scratch directory."""

    def __init__(self, work, nx=8, ny=1, Lx=1.0, Ly=None, Lz=None, x0=0.0, y0=0.0):
        self.dir = Path(work)
        if self.dir.exists():
            shutil.rmtree(str(self.dir))
        (self.dir / 'INPUT').mkdir(parents=True)
        (self.dir / 'OUTPUT').mkdir()
        self.nx, self.ny, self.nz = nx, ny, 1
        self.Lx = Lx
        self.Ly = Ly if Ly is not None else Lx / nx
        self.Lz = Lz if Lz is not None else self.Ly
        self.xn = linspace(x0, x0 + Lx, nx + 1)
        self.yn = linspace(y0, y0 + self.Ly, ny + 1)
        self.zn = [0.0, self.Lz]
        self.xc = [0.5 * (self.xn[i] + self.xn[i + 1]) for i in range(nx)]
        self.yc = [0.5 * (self.yn[j] + self.yn[j + 1]) for j in range(ny)]

    @property
    def centres(self):
        """(x, y) of every cell centre, in the order the solution arrays use."""
        return [(x, y) for y in self.yc for x in self.xc]

    def cells(self, value):
        """Cell-centred field from a constant, a function of x, or a function of (x, y).

        The 2-D form is what the vortex case needs; a 1-D function is still broadcast
        across j, so every case written before y existed keeps working unchanged.
        """
        if callable(value):
            code = getattr(value, '__code__', None)
            if code is not None and code.co_argcount == 2:
                return [value(x, y) for y in self.yc for x in self.xc]
            return [value(x) for _ in range(self.ny) for x in self.xc]
        return [value] * (self.nx * self.ny)

    def particles(self, rho, u, v, w, T, n, families=1, scales=None):
        """Write the condensed-phase initial condition (MK: rho, u, v, w, T, n).

        `families` repeats the same field for each [ICE-Family*] of the run: ICE reads
        the families' variables one after the other, so the file simply carries as many
        six-variable groups as there are families. `scales` multiplies the loading of
        each one, which is how a case can be made exactly proportional between families.
        """
        names, values = [], []
        scale = scales if scales else [1.0] * families
        for p in range(1, families + 1):
            names += ['%s%d' % (v, p) for v in ('rp', 'up', 'vp', 'wp', 'Tp', 'np')]
            fields = [self.cells(f) for f in (rho, u, v, w, T, n)]
            # rho and n carry the loading, so scaling both leaves the radius alone
            for iv in (0, 5):
                fields[iv] = [scale[p-1] * x for x in fields[iv]]
            values += fields
        write_tec(self.dir / 'INPUT/part-ic.tec', names,
                  self.xn, self.yn, self.zn, values, zone='B1-CD')

    def gas(self, rho, u, v, w, T, R, gam, k, mu):
        """Write the frozen carrier phase. Its presence switches on one-way coupling."""
        write_tec(self.dir / 'INPUT/gas.tec',
                  ['rho', 'u', 'v', 'w', 'T', 'R', 'gam', 'k', 'mu', 'dt'],
                  self.xn, self.yn, self.zn,
                  [self.cells(f) for f in (rho, u, v, w, T, R, gam, k, mu, 0.0)], zone='B1-GAS')

    def boundaries(self, mode='extrapolation', inlets=None):
        write_bc(self.dir / 'INPUT/part-bc.txt', self.nx, self.ny, self.nz, mode, inlets)

    def ini(self, t_end=None, iters=1000000000, cfl=0.8, rk='RK2', drag='Stokes',
            heat='Stokes', reconstruction='MUSCL', limiter='vanleer', rho_al=1000.0,
            cs=900.0, dt_max=None, physics=None, closures=('MK',)):
        numerics = {'time-scheme': rk, 'cfl': cfl, 'time-accurate': True,
                    'space-reconstruction': reconstruction, 'flux-limiter': limiter}
        if dt_max is not None:
            numerics['dt-max'] = dt_max
        # `physics` adds to or overrides the [ICE-Physics] defaults, which is how the
        # evaporation case reaches the vapour keys without every other case carrying them
        phys = {'drag': drag, 'heat-transfer': heat,
                'density': rho_al, 'specific-heat': cs, 'emissivity': 0.0}
        phys.update(physics or {})
        sections = {'ICE-Parameters': {'iter-threshold': iters,
                                       'time-threshold': t_end if t_end is not None else 1e30,
                                       'res-threshold': 0.0},
                    'ICE-Numerics': numerics}
        for p, closure in enumerate(closures, start=1):
            sections['ICE-Family%d' % p] = {'closure': closure}
        write_ini(self.dir / 'input.ini', {
            **sections,
            'ICE-Physics': phys,
            'ICE-IO': {'ic-format': 'tecplot ascii', 'sol-format': 'tecplot ascii',
                       'shell-diter': 1000000000, 'sol-diter': 1000000000,
                       'res-diter': 1000000000},
        })

    def run(self, threads=1):
        run_ice(self.dir, threads)
        return read_solution(self.dir / 'OUTPUT/part-field.tec')


def run_ice(case_dir, threads=1):
    """Run the solver in case_dir. Raises if it fails, pointing at the log."""
    if not ICE_BIN.exists():
        raise SystemExit('[verification] no ICE binary at %s - build first' % ICE_BIN)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), KMP_STACKSIZE='100M',
               OMP_STACKSIZE='100M')
    with open(str(case_dir / 'log'), 'w') as out, open(str(case_dir / 'err'), 'w') as err:
        rc = subprocess.call([str(ICE_BIN)], cwd=str(case_dir), stdout=out, stderr=err, env=env)
    if rc != 0:
        raise RuntimeError('ICE failed (exit %d) in %s - see log and err there' % (rc, case_dir))


def read_solution(path):
    """First zone of a Tecplot BLOCK output: {'time', 'nx', 'ny', 'var'[v][cell]}."""
    text = Path(path).read_text()
    zones = re.split(r'^\s*ZONE', text, flags=re.IGNORECASE | re.MULTILINE)[1:]
    header, body = zones[0].split('\n', 1)
    I, J, K = (int(re.search(r'\b%s\s*=\s*(\d+)' % a, header, re.IGNORECASE).group(1)) for a in 'IJK')
    stime = re.search(r'SOLUTIONTIME\s*=\s*([-+0-9.EDed]+)', header, re.IGNORECASE)
    values = [float(t) for t in body.split() if re.fullmatch(r'[-+0-9.Ee]+', t)]
    nn = I * J * K
    nc = (I - 1) * (J - 1) * max(1, K - 1)
    nvar = (len(values) - 3 * nn) // nc
    out = {'time': float(stime.group(1)) if stime else None,
           'nx': I - 1, 'ny': J - 1,
           'var': [values[3 * nn + v * nc: 3 * nn + (v + 1) * nc] for v in range(nvar)]}
    return out


# ---------------------------------------------------------------------------
#  Reference physics - an independent implementation of what ICE integrates
# ---------------------------------------------------------------------------

def drag_coefficient(law, Re, Ma, G, Tr):
    """Cd for each correlation ICE offers (Lib_Drag.f90)."""
    if law == 'Newton':
        return 0.45
    if law == 'Stokes':
        return 24.0 / (Re + TOLL)
    if law == 'Schlichting':
        return 24.0 / (Re + TOLL) * (1.0 + 3.0 * Re / 16.0)
    if law == 'Schiller-Naumann':
        return 24.0 / (Re + TOLL) * (1.0 + 0.15 * Re ** 0.687)
    if law == 'Wen-Yu':
        return 24.0 / (Re + TOLL) * (1.0 + 0.15 * Re ** 0.687) if Re <= 1000.0 else 0.43
    if law == 'Putnam':
        return 24.0 / (Re + TOLL) * (1.0 + Re ** (2.0 / 3.0) / 6.0) if Re < 1000.0 else 0.4392
    if law == 'Clift-Gauvin':
        return 24.0 / (Re + TOLL) * (1.0 + 0.15 * Re ** 0.687
                                     + 0.0175 * Re / (1.0 + 4.25e4 * Re ** (-1.16)))
    if law == 'Morsi-Alexander':
        for hi, a in ((0.1, (0.0, 24.0, 0.0)), (1.0, (3.69, 22.73, 0.0903)),
                      (10.0, (1.222, 29.1667, -3.8889)), (100.0, (0.6167, 46.5, -116.67)),
                      (1000.0, (0.3644, 98.33, -2778.0)), (5000.0, (0.357, 148.62, -4.75e4)),
                      (10000.0, (0.46, -490.546, 57.87e4))):
            if Re <= hi:
                a1, a2, a3 = a
                break
        else:
            a1, a2, a3 = 0.5191, -1662.5, 5.4167e6
        return a1 + a2 / (Re + TOLL) + a3 / (Re * Re + TOLL)
    if law == 'Carlson-Hoglund':
        Cd0 = drag_coefficient('Wen-Yu', Re, Ma, G, Tr)
        return Cd0 * (1.0 + math.exp(-0.427 / (Ma ** 4.63 + TOLL) - 3.0 / (Re ** 0.88 + TOLL))) \
                   / (1.0 + Ma / (Re + TOLL) * (3.82 + 1.28 * math.exp(-1.25 * Re / (Ma + TOLL))))
    if law == 'Henderson':
        if Ma <= 1.0:
            return _henderson_sub(Re, Ma, G, Tr)
        if Ma >= 1.75:
            return _henderson_sup(Re, Ma, G, Tr)
        Cd1 = _henderson_sub(Re, 1.00, G, Tr)
        Cd2 = _henderson_sup(Re, 1.75, G, Tr)
        return Cd1 + 4.0 / 3.0 * (Ma - 1.0) * (Cd2 - Cd1)
    if law in ('Crowe', 'Hermsen'):
        if law == 'Crowe':
            gfun = 10.0 ** (1.25 * (1.0 + math.tanh(0.77 * math.log10(Re + TOLL) - 1.92)))
            hfun = 2.3 + 1.7 * Tr ** 0.5 - 2.3 * math.tanh(1.17 * math.log10(Ma + TOLL))
        else:
            gfun = (1.0 + Re * (12.278 + 0.548 * Re)) / (1.0 + 11.278 * Re)
            hfun = 5.6 / (1.0 + Ma) + 1.7 * Tr ** 0.5
        Cd0 = drag_coefficient('Wen-Yu', Re, Ma, G, Tr)
        return 2.0 + (Cd0 - 2.0) * _exp(-3.07 * G ** 0.5 * Ma / (Re + TOLL) * gfun) \
                   + hfun / (Ma * G ** 0.5 + TOLL) * _exp(-Re / (2 * Ma + TOLL))
    raise ValueError('unknown drag law ' + law)


def _exp(x):
    """exp() that underflows to zero the way Fortran's does."""
    return math.exp(x) if x > -700.0 else 0.0


def _henderson_sub(Re, Ma, G, Tr):
    s = Ma * (0.5 * G) ** 0.5
    return 24.0 * (Re + s * (4.33 + (3.65 - 1.53 * Tr) / (1.0 + 0.353 * Tr)
                             * _exp(-0.247 * Re / (s + TOLL))) + TOLL) ** (-1.0) \
         + _exp(-0.5 * Ma / (Re ** 0.5 + TOLL)) \
           * ((4.5 + 0.38 * (0.03 * Re + 0.48 * Re ** 0.5)) / (1.0 + 0.03 * Re + 0.48 * Re ** 0.5)
              + 0.1 * Ma ** 2 + 0.2 * Ma ** 8) \
         + 0.6 * s * (1.0 - _exp(-Ma / (Re + TOLL)))


def _henderson_sup(Re, Ma, G, Tr):
    return (0.9 + 0.34 / (Ma ** 2 + TOLL)
            + 1.86 * (Ma / (Re + TOLL)) ** 0.5 * (2.0 + 2.0 / (Ma ** 2 * 0.5 * G + TOLL)
                                                  + 1.058 / (Ma * (0.5 * G) ** 0.5 + TOLL) * Tr ** 0.5)
            - 1.0 / (Ma ** 4 * G ** 2 * 0.25 + TOLL)) / (1.0 + 1.86 * (Ma / (Re + TOLL)) ** 0.5)


def nusselt(law, Re, Pr, Ma):
    """Nu for each correlation ICE offers (Lib_Heat.f90)."""
    if law == 'Stokes':
        return 2.0
    if law == 'JAXA1':
        return 2.5 * Re ** 0.15 + 0.04 * Re
    if law == 'JAXA2':
        return 2.0 + 0.37 * Re ** 0.6 * Pr ** (1.0 / 3.0)
    if law == 'JAXA3':
        return 1.0 / (1.0 / (2.0 + 0.645 * Re ** 0.5 * Pr ** (1.0 / 3.0)) + 3.42 * Ma / (Re * Pr))
    if law == 'Chang':
        return 2.0 + 0.459 * Re ** 0.55 * Pr ** (1.0 / 3.0)
    if law == 'Ranz-Marshall':
        return 2.0 + 0.6 * Re ** 0.5 * Pr ** (1.0 / 3.0)
    if law == 'Kavanau-Drake':
        Nu = 2.0 + 0.459 * Re ** 0.55 * Pr ** 0.33
        return Nu / (1.0 + 3.42 * Ma / (Re * Pr + TOLL) * Nu)
    raise ValueError('unknown heat law ' + law)


class Physics(object):
    """The uniform gas and particle cloud shared by the relaxation cases.

    ICE has no volume fraction: the state is the bulk density rho_p = alpha * rho_al
    and the number density n, from which it recovers the particle radius. The values
    here are chosen so that radius comes back as dp/2.
    """

    def __init__(self, rho_g=1.2, ug=10.0, vg=0.0, Tg=300.0, R=287.05, gam=1.4,
                 mu=1.8e-5, kg=0.026, rho_al=1000.0, cs=900.0, dp=1.0e-4, alpha=1.0e-3):
        self.rho_g, self.ug, self.vg, self.Tg = rho_g, ug, vg, Tg
        self.R, self.gam, self.mu, self.kg = R, gam, mu, kg
        self.rho_al, self.cs, self.dp, self.alpha = rho_al, cs, dp, alpha
        self.rp = 0.5 * dp
        self.rho_p = alpha * rho_al
        self.n = self.rho_p / (rho_al * (4.0 / 3.0) * math.pi * self.rp ** 3)
        self.Pr = gam / (gam - 1.0) * R * mu / kg
        self.sound = math.sqrt(gam * R * Tg)

    @property
    def tau_stokes(self):
        """rho_al * dp^2 / (18 mu): the Stokes relaxation time."""
        return self.rho_al * self.dp ** 2 / (18.0 * self.mu)

    @property
    def tau_thermal(self):
        """rho_al * cs * dp^2 / (12 kg): thermal relaxation at Nu = 2."""
        return self.rho_al * self.cs * self.dp ** 2 / (12.0 * self.kg)

    def tau_drag(self, law, slip, Tp):
        """ICE's relaxation time: 8 rho_al Rp / (3 rho_g Cd |du|), as in Mod_Sources."""
        Re = 2.0 * self.rho_g * self.rp * abs(slip) / self.mu
        Ma = abs(slip) / self.sound
        Cd = drag_coefficient(law, Re, Ma, self.gam, Tp / self.Tg)
        return 8.0 * self.rho_al * self.rp / (3.0 * self.rho_g * Cd * abs(slip) + 1.0e-20)

    def heat_rate(self, law, slip, Tp):
        """dT/dt from ICE's convective exchange term, 2 Nu k pi Rp n (Tg - Tp) / (rho_p cs)."""
        Re = 2.0 * self.rho_g * self.rp * abs(slip) / self.mu
        Ma = abs(slip) / self.sound
        Nu = nusselt(law, Re, self.Pr, Ma)
        return 2.0 * Nu * self.kg * math.pi * self.rp * self.n * (self.Tg - Tp) \
            / (self.rho_p * self.cs)

    def reference(self, t_end, drag='Stokes', heat='Stokes', up0=0.0, Tp0=None, steps=5000):
        """Integrate the (u, T) relaxation ODE with RK4 - the reference solution."""
        Tp0 = self.Tg if Tp0 is None else Tp0
        h = t_end / steps
        u, T = up0, Tp0

        def rhs(u, T):
            slip = self.ug - u
            return (slip / self.tau_drag(drag, slip, T), self.heat_rate(heat, slip, T))

        for _ in range(steps):
            k1 = rhs(u, T)
            k2 = rhs(u + 0.5 * h * k1[0], T + 0.5 * h * k1[1])
            k3 = rhs(u + 0.5 * h * k2[0], T + 0.5 * h * k2[1])
            k4 = rhs(u + h * k3[0], T + h * k3[1])
            u += h / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            T += h / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        return u, T


# ---------------------------------------------------------------------------
#  Reporting
# ---------------------------------------------------------------------------

def observed_order(errors, sizes):
    """Order between consecutive refinements: log(e1/e2) / log(h1/h2)."""
    out = []
    for i in range(len(errors) - 1):
        if errors[i + 1] <= 0.0:
            out.append(float('inf'))
        else:
            out.append(math.log(errors[i] / errors[i + 1]) / math.log(float(sizes[i + 1]) / sizes[i]))
    return out


class Report(object):
    """Collects checks and turns them into one PASS/FAIL line and an exit code."""

    def __init__(self, title):
        self.title = title
        self.failures = 0

    def check(self, ok, text):
        print('   %s %s' % ('OK  ' if ok else 'FAIL', text))
        if not ok:
            self.failures += 1
        return ok

    def close(self, work=None):
        result = '%sPASS%s' % (GREEN, RESET) if not self.failures else '%sFAIL%s' % (RED, RESET)
        print('%-36s -->  %s' % (self.title, result))
        if self.failures == 0 and work is not None and os.environ.get('ICE_KEEP_WORK') is None:
            shutil.rmtree(str(work), ignore_errors=True)
        sys.exit(1 if self.failures else 0)
