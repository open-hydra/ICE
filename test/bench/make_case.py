#!/usr/bin/env python3
"""Write a 3-D multi-block, multi-family ICE case for the parallel gates and legs.

    make_case.py OUTDIR [--preset small|cache|bench3] [options]

The campaign case of test/bench/monolith (192^3, one block per rank, one IG
family) never exercised what a production mesh does to the parallel code: several
blocks per rank, several families, and block corners whose cells carry three
boundary records. This generator writes such a case from a handful of numbers:

  * an aligned box of NX x NY x NZ cells cut into NBX x NBY x NBZ blocks, every
    cut a connection (101), every periodic wrap a connection back across the
    domain (201), the remaining faces symmetry (300) or extrapolation (400);
  * one material (`A <families>`) with one family per closure named in
    --families, in that order; family p carries 2^(p-1) times the loading of
    family 1, so the families are distinguishable and a bug that touches only
    the last one shows;
  * a density that is above rho_empty EVERYWHERE and varies in all three
    directions, rho0 (1 + amp cos 2pi x/Lx cos 2pi y/Ly cos 2pi z/Lz), with a
    uniform velocity that differs per direction: every interface and every
    corner carries changing data from iteration 1, and a transposed halo shows;
  * optionally a frozen Taylor-Green carrier (one-way coupling).

Everything is written in the ASCII formats the verification cases use, through
test/verification/common.py, so Python 3.6 with the standard library is enough.
`--solo K` writes only family K of the list, as family 1, with the same state:
the oracle of the FamilyHalo3D gate (families are independent within a step, so
family K alone must equal family K in company, bit for bit, when dt-max binds
every cell).

case.json records every parameter, the cell and record counts, and the number
of iterations after which a signal has crossed one block width:
N = ceil(1.25 nx_block dx / (Ux dt_max)).
"""
import argparse
import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'verification'))
from common import write_tec_blocks, write_ini, linspace  # noqa: E402

PRESETS = {
    # gates: seconds, three closures, corners with three boundary records
    'small':  dict(cells=(24, 16, 12), blocks=(3, 2, 2), families='MK,IG,AG',
                   bc=('periodic', 'periodic', 'sym-out'), coupled='off'),
    # L3-resident overhead leg
    'cache':  dict(cells=(32, 32, 16), blocks=(1, 1, 1), families='IG',
                   bc=('periodic', 'periodic', 'periodic'), coupled='off'),
    # the multi-family factor of the ghost pipeline, 24 blocks
    'bench3': dict(cells=(96, 64, 48), blocks=(4, 2, 3), families='MK,IG,AG',
                   bc=('periodic', 'periodic', 'periodic'), coupled='taylor-green'),
}

NAMES = {
    'MK': ['rp', 'up', 'vp', 'wp', 'Tp', 'np'],
    'IG': ['rp', 'up', 'vp', 'wp', 'Pp', 'Tp', 'np'],
    'AG': ['rp', 'up', 'vp', 'wp', 'Pp{p}1', 'Pp{p}12', 'Pp{p}13', 'Pp{p}22',
           'Pp{p}23', 'Pp{p}33', 'Tp', 'np'],
}

FACE_CODE = {'sym': 300, 'out': 400}


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('outdir')
    ap.add_argument('--preset', choices=sorted(PRESETS), default='small')
    ap.add_argument('--cells', nargs=3, type=int, metavar=('NX', 'NY', 'NZ'))
    ap.add_argument('--blocks', nargs=3, type=int, metavar=('NBX', 'NBY', 'NBZ'))
    ap.add_argument('--length', nargs=3, type=float, metavar=('LX', 'LY', 'LZ'))
    ap.add_argument('--dx', type=float, default=1.0e-3, help='cell size when --length is absent')
    ap.add_argument('--families', help='comma-separated closures, one family each (MK, IG, AG)')
    ap.add_argument('--solo', type=int, default=0, help='write only family K of the list, as family 1')
    ap.add_argument('--bc-x', choices=['periodic', 'sym', 'sym-out'])
    ap.add_argument('--bc-y', choices=['periodic', 'sym', 'sym-out'])
    ap.add_argument('--bc-z', choices=['periodic', 'sym', 'sym-out'])
    ap.add_argument('--velocity', nargs=3, type=float, default=(10.0, 7.0, 4.0))
    ap.add_argument('--amp', type=float, default=0.3)
    ap.add_argument('--rho0', type=float, default=1.0e-3)
    ap.add_argument('--a-disp', type=float, default=10.0, help='dispersion speed a: P = rho a^2 / 3')
    ap.add_argument('--radius', type=float, default=1.0e-5)
    ap.add_argument('--rho-mat', type=float, default=1000.0)
    ap.add_argument('--T0', type=float, default=300.0)
    ap.add_argument('--coupled', choices=['off', 'taylor-green'])
    ap.add_argument('--U0', type=float, default=1.0, help='Taylor-Green amplitude')
    ap.add_argument('--scheme', default='RK2')
    ap.add_argument('--cfl', type=float, default=0.5)
    ap.add_argument('--dt-max', type=float, default=1.0e-5)
    ap.add_argument('--irs', choices=['off', 'on'], default='off')
    ap.add_argument('--shock-detector', default='none')
    ap.add_argument('--iters', type=int, default=0, help='0: the crossing count from case.json')
    ap.add_argument('--shell-diter', type=int, default=10)
    ap.add_argument('--timers', choices=['off', 'on'], default='off')
    ap.add_argument('--sol-format', default='vtk raw')
    a = ap.parse_args(argv)

    p = PRESETS[a.preset]
    a.cells = tuple(a.cells or p['cells'])
    a.blocks = tuple(a.blocks or p['blocks'])
    a.families = (a.families or p['families']).split(',')
    a.bc = (a.bc_x or p['bc'][0], a.bc_y or p['bc'][1], a.bc_z or p['bc'][2])
    a.coupled = a.coupled or p['coupled']
    if a.length is None:
        a.length = tuple(n * a.dx for n in a.cells)
    for c in a.families:
        if c not in NAMES:
            sys.exit('unknown closure %s' % c)
    for d in range(3):
        if a.cells[d] % a.blocks[d]:
            sys.exit('%d cells do not split into %d blocks in direction %d' % (a.cells[d], a.blocks[d], d + 1))
        if a.cells[d] // a.blocks[d] < 4:
            sys.exit('fewer than 4 cells per block in direction %d: MUSCL and the second ghost need 4' % (d + 1))
    if a.solo and not 1 <= a.solo <= len(a.families):
        sys.exit('--solo %d: the list has %d families' % (a.solo, len(a.families)))
    return a


def block_id(bi, bj, bk, nb):
    """1-based block number, bi fastest: the order ATLAS and MDB use."""
    return 1 + bi + nb[0] * (bj + nb[1] * bk)


def block_nodes(a, bi, bj, bk):
    """Node coordinates of block (bi, bj, bk)."""
    out = []
    for d, b in enumerate((bi, bj, bk)):
        n = a.cells[d] // a.blocks[d]
        h = a.length[d] / a.cells[d]
        out.append(linspace(b * n * h, (b + 1) * n * h, n + 1))
    return out


def centres(nodes):
    return [0.5 * (nodes[i] + nodes[i + 1]) for i in range(len(nodes) - 1)]


def family_fields(a, closure, p, xc, yc, zc):
    """Cell-centred fields of one family in i-fastest order, family p scaled 2^(p-1)."""
    scale = 2.0 ** (p - 1)
    kx, ky, kz = (2.0 * math.pi / L for L in a.length)
    vol = a.rho_mat * 4.0 / 3.0 * math.pi * a.radius ** 3
    P0 = a.rho0 * a.a_disp ** 2 / 3.0
    U, V, W = a.velocity
    rho, u, v, w, T, n = [], [], [], [], [], []
    for z in zc:
        for y in yc:
            for x in xc:
                r = scale * a.rho0 * (1.0 + a.amp * math.cos(kx * x) * math.cos(ky * y) * math.cos(kz * z))
                rho.append(r); u.append(U); v.append(V); w.append(W); T.append(a.T0); n.append(r / vol)
    zeros = [0.0] * len(rho)
    press = [P0] * len(rho)
    if closure == 'MK':
        return [rho, u, v, w, T, n]
    if closure == 'IG':
        return [rho, u, v, w, press, T, n]
    return [rho, u, v, w, press, zeros, zeros, press, zeros, press, T, n]


def family_names(closure, p):
    return [nm.format(p=p) + ('' if '{p}' in nm else str(p)) for nm in NAMES[closure]]


def gas_fields(a, xc, yc, zc):
    """Frozen Taylor-Green carrier: rho u v w T R gam k mu dt."""
    kx, ky, kz = (2.0 * math.pi / L for L in a.length)
    u, v, w = [], [], []
    for z in zc:
        for y in yc:
            for x in xc:
                u.append(a.U0 * math.sin(kx * x) * math.cos(ky * y) * math.cos(kz * z))
                v.append(-a.U0 * math.cos(kx * x) * math.sin(ky * y) * math.cos(kz * z))
                w.append(0.0)
    m = len(u)
    return [[1.2] * m, u, v, w, [300.0] * m, [287.05] * m, [1.4] * m, [0.026] * m, [1.8e-5] * m, [0.0] * m]


def bc_records(a):
    """The ATLAS-format boundary table, one copy shared by every family."""
    nb = a.blocks
    nloc = [a.cells[d] // nb[d] for d in range(3)]
    rows = []

    def header(b, i, j, k, f, code):
        return '%8d%8d%8d%8d%8d%8d' % (b, i, j, k, f, code)

    def payload(bs, i, j, k, fs):
        return '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (bs, i, j, k, fs, 1, 0, 0, 1)

    for bk in range(nb[2]):
        for bj in range(nb[1]):
            for bi in range(nb[0]):
                b = block_id(bi, bj, bk, nb)
                pos = (bi, bj, bk)
                nx, ny, nz = nloc
                for d in range(3):
                    lo, hi = (1, nloc[d])
                    for f, side in ((2 * d + 1, -1), (2 * d + 2, +1)):
                        nbr = list(pos)
                        nbr[d] += side
                        code = None
                        if 0 <= nbr[d] < nb[d]:
                            code = 101
                        elif a.bc[d] == 'periodic':
                            nbr[d] %= nb[d]
                            code = 201
                        elif a.bc[d] == 'sym':
                            code = 300
                        else:                                  # sym-out
                            code = 300 if side < 0 else 400
                        bs = block_id(nbr[0], nbr[1], nbr[2], nb) if code in (101, 201) else 0
                        own = lo if side < 0 else hi          # the boundary cell's index along d
                        src = hi if side < 0 else lo          # the source cell's index along d
                        fs = f + 1 if side < 0 else f - 1     # the neighbour's facing face
                        # the two transverse directions, inner index fastest
                        t1, t2 = [e for e in range(3) if e != d]
                        for k2 in range(1, nloc[t2] + 1):
                            for k1 in range(1, nloc[t1] + 1):
                                ijk = [0, 0, 0]
                                ijk[d], ijk[t1], ijk[t2] = own, k1, k2
                                rows.append(header(b, ijk[0], ijk[1], ijk[2], f, code))
                                if code in (101, 201):
                                    sjk = list(ijk)
                                    sjk[d] = src
                                    rows.append(payload(bs, sjk[0], sjk[1], sjk[2], fs))
    return rows


def build_case(argv):
    a = parse_args(argv)
    out = Path(a.outdir)
    (out / 'INPUT').mkdir(parents=True, exist_ok=True)
    (out / 'OUTPUT').mkdir(exist_ok=True)
    nb = a.blocks

    families = list(enumerate(a.families, start=1))
    if a.solo:
        families = [(a.solo, a.families[a.solo - 1])]

    part_blocks, gas_blocks, names = [], [], []
    for pw, (p, closure) in enumerate(families, start=1):
        names += family_names(closure, pw)
    for bk in range(nb[2]):
        for bj in range(nb[1]):
            for bi in range(nb[0]):
                b = block_id(bi, bj, bk, nb)
                xn, yn, zn = block_nodes(a, bi, bj, bk)
                xc, yc, zc = centres(xn), centres(yn), centres(zn)
                fields = []
                for p, closure in families:
                    fields += family_fields(a, closure, p, xc, yc, zc)
                part_blocks.append((xn, yn, zn, fields, 'B%d-CD' % b))
                if a.coupled != 'off':
                    gas_blocks.append((xn, yn, zn, gas_fields(a, xc, yc, zc), 'B%d-GAS' % b))

    write_tec_blocks(out / 'INPUT/part-ic.tec', names, part_blocks)
    if gas_blocks:
        write_tec_blocks(out / 'INPUT/gas.tec', ['rho', 'u', 'v', 'w', 'T', 'R', 'gam', 'k', 'mu', 'dt'], gas_blocks)
    rows = bc_records(a)
    (out / 'INPUT/part-bc.txt').write_text('\n'.join(rows) + '\n')
    (out / 'INPUT/part-phase.txt').write_text('condensed-dispersed phase\nA %d\n' % len(families))

    # iterations after which a signal has crossed one block width in x
    dx = a.length[0] / a.cells[0]
    nxb = a.cells[0] // nb[0]
    crossing = int(math.ceil(1.25 * nxb * dx / (abs(a.velocity[0]) * a.dt_max)))
    iters = a.iters or crossing
    wave = max(abs(v) for v in a.velocity) + a.a_disp
    dt_cfl = a.cfl * min(a.length[d] / a.cells[d] for d in range(3)) / wave

    numerics = {
        'time-scheme': a.scheme, 'cfl': a.cfl, 'time-accurate': True,
        'dt-max': a.dt_max, 'tau-factor': 0.0,
        'space-reconstruction': 'MUSCL', 'flux-limiter': 'vanleer',
        'shock-detector': a.shock_detector, 'irs': a.irs == 'on',
    }
    sections = {
        'ICE-Parameters': {'iter-threshold': iters, 'time-threshold': 1.0e30, 'res-threshold': 0.0},
        'ICE-Numerics': numerics,
        'ICE-Physics': {'drag': 'Schiller-Naumann', 'heat-transfer': 'Ranz-Marshall',
                        'density': a.rho_mat, 'specific-heat': 4180.0, 'emissivity': 0.0},
        'ICE-IO': {'ic-format': 'tecplot ascii', 'sol-format': a.sol_format,
                   'shell-diter': a.shell_diter, 'sol-diter': 1000000000,
                   'res-diter': 1000000000, 'timers': a.timers == 'on'},
    }
    for pw, (p, closure) in enumerate(families, start=1):
        sections['ICE-Family%d' % pw] = {'closure': closure}
    write_ini(out / 'input.ini', sections)

    ncells = a.cells[0] * a.cells[1] * a.cells[2]
    info = {
        'preset': a.preset, 'cells': a.cells, 'blocks': nb, 'length': a.length,
        'families': [c for _, c in families], 'family_ids': [p for p, _ in families],
        'solo': a.solo, 'bc': a.bc, 'coupled': a.coupled, 'velocity': a.velocity,
        'amp': a.amp, 'rho0': a.rho0, 'a_disp': a.a_disp, 'radius': a.radius, 'T0': a.T0,
        'scheme': a.scheme, 'cfl': a.cfl, 'dt_max': a.dt_max, 'dt_cfl_estimate': dt_cfl,
        'dt_max_binds': a.dt_max < dt_cfl, 'iters': iters, 'crossing_iters': crossing,
        'ncells': ncells, 'nblocks': nb[0] * nb[1] * nb[2],
        'bc_records': sum(1 for r in rows if len(r.split()) == 6),
    }
    (out / 'case.json').write_text(json.dumps(info, indent=1, sort_keys=True) + '\n')
    return info


if __name__ == '__main__':
    info = build_case(sys.argv[1:])
    print('wrote %s: %d cells, %d blocks, families %s, %d boundary records, %d iterations (dt-max binds: %s)'
          % (sys.argv[1], info['ncells'], info['nblocks'], ','.join(info['families']),
             info['bc_records'], info['iters'], info['dt_max_binds']))
