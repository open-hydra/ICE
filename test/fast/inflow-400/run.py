"""Fast test: the run says where particles enter through extrapolation faces (400).

A 400 face copies its boundary cell into its ghost, so particles moving into the domain enter with the state of that
cell. The run prints a line at the first stage in which a cell that carries the phase (density above rho_empty) moves
into the domain faster than round-off (above 1e-12 of the signal speed for IG and AG, of its speed for MK, and above
the floor's velocity 1e-25 m/s) through the 400 faces of a face of a block (the
iteration, the faces, the largest |w|, w the velocity along the normal out of the domain; |w|/sqrt(P/rho) for IG,
|w|/sqrt(Pnn/rho) for AG with Pnn the pressure along that normal), one more each time that largest value at a stage is
more than 10 times that of the last line, and at the end of the run one line per such boundary with the stages it let
particles in. On a 1-D duct (16 cells, 400 on faces 1 to 4, 0 on 5 and 6), 20 time-accurate RK2 steps:

  mk        MK moving along -x at 10 m/s: it enters through face 2 only, at iteration 1, one face, |w| = 10 m/s,
            in every one of the 40 stages; the cloud stays uniform, so the flux is that of the cell.
  ig, ag    the same with IG (P = 25, rho = 1: |w|/sqrt(P/rho) = 2) and AG (P11 = 100, P22 = P33 = 1: the normal of
            face 2 is x, so |w|/sqrt(Pnn/rho) = 1).
  none      a 402 inlet on face 1 and the cloud moving along +x: nothing enters through a 400 face, no line.
  empty     MK in cells 5 to 12 moving along -x, the others empty: no line. An empty cell takes the density floor as
            its density and its velocity, which points into the domain through faces 1 and 3, and the cells next to
            it carry round-off from it in their velocity, far below the floor.
  sparse    MK at 1e-8 kg/m3, below ICE's empty-cell density rho_empty, moving into the domain through face 2: no
            line, the cells do not carry the phase.
  inclined  empty with face 1 inclined: no line. Its empty cell moves into the domain along the inclined normal, and
            the inclined cell, once the cloud is in it, carries round-off across the duct (|w| about 1e-17 of its
            10 m/s), below 1e-12 of its speed.
  rest      MK at rest in cells 5 to 12, the others empty: no line. The round-off the floor leaves in the cells at
            rest points into the domain, but is slower than the floor's velocity.
  still     IG and AG at rest up to round-off, u = -1e-15 m/s (2e-16 and 1e-16 of the signal speed sqrt(Pnn/rho)),
            beside face 2: no line. All of their velocity points into the domain, so only the bar against the
            signal speed takes them out.
  growth    a cloud at rest that a gas moving along -x drags in through face 2 (dt-max far below the drag time): the
            first line at iteration 1 and another with |w| more than 10 times the first; 39 stages (not the first).
  threads   mk on 2 threads: the same lines.
  levels    mk on 16 x 2 cells with 2 grid levels (5 steps on level 2, 10 on level 1, steady): a line per level,
            with the level, at iterations 1 and 6, and 10 and 20 stages.
"""
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'verification'))
from common import ICE_BIN, Report, inlet_payload, linspace, write_ini, write_tec  # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
NX, LX, STEPS, U = 16, 0.16, 20, 10.0
RHO_AL, CS, RADIUS = 1000.0, 900.0, 1e-5
MASS = RHO_AL * 4.0 / 3.0 * math.pi * RADIUS**3

FIRST = re.compile(r' 400 boundary of block (\d+), face (\d+)(?:, family \d+)?(?:, grid level (\d+))?: particles entering '
                   r'at iteration (-?\d+), faces(?: above 10 times its last line)?: (\d+), largest (\S+): (\S+)')
TOTAL = re.compile(r' 400 boundary of block (\d+), face (\d+)(?:, family \d+)?(?:, grid level (\d+))?: particles entered in '
                   r'(\d+) stages, first at iteration (-?\d+) \(faces: (\d+), largest (\S+): (\S+?)(?: m/s)?\), largest '
                   r'over the run: (\S+)')


def bc(path, nx, ny, codes, inlet=None):
    """One block of nx x ny x 1 cells, face f carrying codes[f]; a 402 record carries inlet as its payload."""
    rows = []
    for f, ii, jj in ((1, [1], range(1, ny + 1)), (2, [nx], range(1, ny + 1)),
                      (3, range(1, nx + 1), [1]), (4, range(1, nx + 1), [ny]),
                      (5, range(1, nx + 1), range(1, ny + 1)), (6, range(1, nx + 1), range(1, ny + 1))):
        for j in jj:
            for i in ii:
                rows.append('%8d%8d%8d%8d%8d%8d' % (1, i, j, 1, f, codes[f]))
                if codes[f] == 402:
                    rows.append('   ' + inlet)
    Path(path).write_text('\n'.join(rows) + '\n')


def write_mesh(path, names, xn, yn, zn, cells, tilt):
    """write_tec with the nodes of i = 0 moved by -tilt*y along x: face 1 inclined about z."""
    lines = [' VARIABLES ="x" "y" "z" ' + ' '.join('"%s"' % v for v in names),
             ' ZONE  T = B1-CD, I=%d, J=%d, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED)'
             % (len(xn), len(yn), 3 + len(names))]
    nodes = [(x - (tilt * y if i == 0 else 0.0), y, z) for z in zn for y in yn for i, x in enumerate(xn)]
    for c in range(3):
        lines += [' %.15E' % node[c] for node in nodes]
    for cell in cells:
        lines += [' %.15E' % a for a in cell]
    Path(path).write_text('\n'.join(lines) + '\n')


def build(name, closure='MK', u=-U, cells=None, codes=None, inlet=None, ny=1, levels=1, tilt=0.0, gas_u=None,
          rho0=None):
    """The duct, a uniform cloud of density rho0 moving along x at u (cells: the 1-based cells that hold it, every
    one by default)."""
    work = WORK / name
    if work.exists():
        shutil.rmtree(str(work))
    (work / 'INPUT').mkdir(parents=True)
    (work / 'OUTPUT').mkdir()
    xn, yn, zn = linspace(0.0, LX, NX + 1), linspace(0.0, ny * LX / NX, ny + 1), [0.0, LX / NX]
    full = [1.0 if cells is None or i + 1 in cells else 0.0 for _ in range(ny) for i in range(NX)]
    if rho0 is None:
        rho0 = 1e-3 if closure == 'MK' else 1.0
    rho = [rho0 * a for a in full]
    vel = [u * a for a in full]
    zero, temp, n = [0.0] * len(full), [300.0] * len(full), [r / MASS for r in rho]
    if closure == 'MK':
        names, fields = ['rp1', 'up1', 'vp1', 'wp1', 'Tp1', 'np1'], [rho, vel, zero, zero, temp, n]
    elif closure == 'IG':
        names, fields = ['rho', 'u', 'v', 'w', 'P', 'T', 'n'], [rho, vel, zero, zero, [25.0 * a for a in full], temp, n]
    else:
        p11, p22 = [100.0 * a for a in full], list(full)
        names = ['rho', 'u', 'v', 'w', 'P11', 'P12', 'P13', 'P22', 'P23', 'P33', 'T', 'n']
        fields = [rho, vel, zero, zero, p11, zero, zero, p22, zero, p22, temp, n]
    write_mesh(work / 'INPUT/part-ic.tec', names, xn, yn, zn, fields, tilt)
    if gas_u is not None:
        c = len(full)
        write_tec(work / 'INPUT/gas.tec', ['rho', 'u', 'v', 'w', 'T', 'R', 'gam', 'k', 'mu', 'dt'], xn, yn, zn,
                  [[1.2] * c, [gas_u] * c, zero, zero, temp, [287.05] * c, [1.4] * c, [0.026] * c, [1.8e-5] * c, zero],
                  zone='B1-GAS')
    codes = codes or {1: 400, 2: 400, 3: 400, 4: 400, 5: 0, 6: 0}
    for m in range(1, levels + 1):
        r = 2 ** (m - 1)
        bc(work / 'INPUT' / ('part-bc.txt' if m == 1 else 'part-bc%d.txt' % m), NX // r, max(ny // r, 1), codes, inlet)
    (work / 'INPUT/part-phase.txt').write_text('condensed-dispersed phase\nA 1\n')
    numerics = {'time-scheme': 'RK2', 'cfl': 0.5, 'time-accurate': levels == 1,
                'space-reconstruction': 'MUSCL', 'flux-limiter': 'vanleer'}
    if gas_u is not None:
        numerics['dt-max'] = 1e-6
    sections = {
        'ICE-Parameters': {'iter-threshold': STEPS if levels == 1 else 10, 'time-threshold': 1e30, 'res-threshold': 0.0},
        'ICE-Numerics': numerics,
        'ICE-Family1': {'closure': closure},
        'ICE-Physics': {'drag': 'Stokes' if gas_u is not None else 'NoDrag', 'heat-transfer': 'NoHeat',
                        'density': RHO_AL, 'specific-heat': CS, 'emissivity': 0.0},
        'ICE-IO': {'ic-format': 'tecplot ascii', 'sol-format': 'tecplot ascii', 'shell-diter': 1000000000,
                   'sol-diter': 1000000000, 'res-diter': 1000000000},
    }
    if levels > 1:
        sections['ICE-Multigrid'] = {'levels': levels, 'level1-iter': 10, 'level2-iter': 5}
    write_ini(work / 'input.ini', sections)
    return work


def run(work, threads=1):
    """The lines of the 400 report in ICE's log: (first-stage and growth lines, end-of-run lines, all of them)."""
    if not ICE_BIN.exists():
        raise SystemExit('[fast] no ICE binary at %s -- build first' % ICE_BIN)
    with open(str(work / 'log'), 'w') as out:
        rc = subprocess.call([str(ICE_BIN)], cwd=str(work), stdout=out, stderr=subprocess.STDOUT,
                             env=dict(os.environ, OMP_NUM_THREADS=str(threads), KMP_STACKSIZE='100M'))
    log = (work / 'log').read_text(errors='replace')
    if rc != 0 or 'Time of operation' not in log:
        raise SystemExit('[fast] ICE failed (exit %d) in %s -- see the log there' % (rc, work))
    lines = [l for l in log.split('\n') if l.startswith(' 400 boundary of ')]
    return ([FIRST.match(l).groups() for l in lines if FIRST.match(l)],
            [TOTAL.match(l).groups() for l in lines if TOTAL.match(l)], lines)


def near(text, value):
    return abs(float(text) - value) <= 1e-4 * value


def one_boundary(r, name, got, label, x, stages):
    """The report of a duct whose cloud enters through face 2 of block 1 only, from iteration 1 in every stage."""
    first, total, lines = got
    ok = len(lines) == 2 and len(first) == 1 and len(total) == 1
    if ok:
        b, f, lev, it, nf, lab, val = first[0]
        ok = (b, f, it, nf, lab) == ('1', '2', '1', '1', label) and near(val, x)
        b, f, lev, ns, it, nf, lab, val, run_max = total[0]
        ok = ok and (b, f, ns, it, nf, lab) == ('1', '2', str(stages), '1', '1', label) and near(val, x) \
            and near(run_max, x)
    r.check(ok, '%-7s one line at iteration 1 (face 2, 1 face, %s %g) and %d stages at the end: %s'
            % (name, label, x, stages, lines if lines else 'no line'))
    return lines


def main():
    r = Report('Inflow through 400 faces')
    mk = one_boundary(r, 'mk', run(build('mk')), '|w|', U, 2 * STEPS)
    one_boundary(r, 'ig', run(build('ig', 'IG')), '|w|/sqrt(P/rho)', U / 5.0, 2 * STEPS)
    one_boundary(r, 'ag', run(build('ag', 'AG')), '|w|/sqrt(Pnn/rho)', U / 10.0, 2 * STEPS)

    inlet = inlet_payload(1.0e-2, U, 300.0, RADIUS)
    lines = run(build('none', u=U, codes={1: 402, 2: 400, 3: 400, 4: 400, 5: 0, 6: 0}, inlet=inlet))[2]
    r.check(not lines, 'none    a 402 inlet and particles leaving through 400 faces: no line %s' % (lines or ''))
    lines = run(build('empty', cells=range(5, 13)))[2]
    r.check(not lines, 'empty   empty cells at the 400 faces (density and velocity at the floor): no line %s'
            % (lines or ''))
    lines = run(build('sparse', rho0=1e-8))[2]
    r.check(not lines, 'sparse  a cloud below rho_empty moving into the domain: no line %s' % (lines or ''))
    lines = run(build('inclined', cells=range(5, 13), tilt=0.5))[2]
    r.check(not lines, 'inclined face 1 inclined (its empty cell moves in), round-off across the inclined cell: no line '
            '%s' % (lines or ''))
    lines = run(build('rest', u=0.0, cells=range(5, 13)))[2]
    r.check(not lines, 'rest    a cloud at rest beside empty cells: no line %s' % (lines or ''))
    lines = run(build('still-ig', 'IG', u=-1e-15))[2] + run(build('still-ag', 'AG', u=-1e-15))[2]
    r.check(not lines, 'still   IG and AG at rest up to round-off (1e-16 of the signal speed): no line %s'
            % (lines or ''))

    first, total, lines = run(build('growth', u=0.0, gas_u=-U))
    ok = len(first) >= 2 and len(total) == 1 and all(g[:2] == ('1', '2') for g in first + total)
    if ok:
        ok = first[0][3] == '1' and 'above 10 times' not in lines[0] and 'above 10 times its last line' in lines[1] \
            and float(first[1][6]) > 10.0 * float(first[0][6]) and total[0][3] == str(2 * STEPS - 1)
    r.check(ok, 'growth  a first line at iteration 1, a second above 10 times it, %d stages: %s' % (2 * STEPS - 1, lines))

    lines = run(build('threads'), threads=2)[2]
    r.check(lines == mk, 'threads 2 threads: the lines of 1 thread')

    first, total, lines = run(build('levels', ny=2, levels=2))
    want_first = [('1', '2', '2', '1', '1'), ('1', '2', '1', '6', '2')]
    want_total = [('1', '2', '2', '10', '1', '1'), ('1', '2', '1', '20', '6', '2')]
    ok = [g[:5] for g in first] == want_first and [g[:6] for g in total] == want_total
    r.check(ok, 'levels  level 2 from iteration 1 (1 face) in 10 stages, level 1 from 6 (2 faces) in 20: %s' % lines)
    r.close(WORK)


if __name__ == '__main__':
    main()
