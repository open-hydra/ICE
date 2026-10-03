#!/usr/bin/env python3
"""Fast tests: a 3-D multi-block, multi-family case must give the same answer
whatever the thread count, the rank count, or the company of its families.

    ICE_EQ_MODE=threads|groups|ranks|hybrid|family [ICE_EQ_ITERS=n] [ICE_EQ_MASK_CORNERS=1] [ICE_EQ_EXTRA='--scheme euler'] python3 -B run.py

The gates before this one (openmp-equiv, mpi-equiv) run 2-D single-family cases
in which no populated cell crosses a block interface within their 200 iterations,
and compare a Tecplot ASCII file written with 15 significant digits. This one
runs the `small` preset of test/bench/make_case.py -- 24x16x12 cells in 3x2x2
blocks, three closures, periodic in x and y, symmetry/outlet in z, a density that
is nowhere empty and moves across every interface from iteration 1 -- and
compares the raw-binary VTK payload byte for byte.

  threads   1 thread against 2, 3, 4 and 8.
            RED on the code before the owner-computes boundary loop, and it was
            measured RED on 2026-09-30: a block-corner cell carries three
            boundary records, `compute_bound` accumulates them with `!$OMP
            ATOMIC` in whatever order the dynamic schedule hands the threads,
            and (a+b)+c is not (a+c)+b. A cell with two records is safe, since
            the residual is zero when they land and a+b commutes.
            ICE_EQ_ITERS=1 ICE_EQ_MASK_CORNERS=1 ICE_EQ_EXTRA='--scheme euler'
            compares ONE STAGE with the eight corner cells of every block left
            out: GREEN today, and the proof that nothing else in the step
            depends on the thread count. A second stage already carries a
            corner's difference into its stencil neighbours, so the mask is
            exact only for a single-stage step; the full-length RK2 leg is the
            target of the owner-computes rewrite.
  groups    1 thread against thread groups (ICE_THREAD_GROUPS): 2 groups of 2
            threads, 3 of 2, 4 of 2, and 2 uneven groups of 1 and 2 threads -- the
            12 blocks owned by groups, first touch, per-cell phases, ghost fill,
            boundary fluxes and the flux kernel's rounds all by group. Each run's
            banner must report the groups in use, or the leg FAILs: a binary that
            ignored the variable would only repeat the threads leg.
            Measured RED on 2026-10-03: a group's ghost records left out (NaN
            ghosts) and every group sweeping group 0's blocks in the flux rounds.
  ranks     1 rank against 2, 3, 4 and 5 (12 blocks, so 5 is uneven), and a
            2x1x1 layout on 3 ranks (more ranks than blocks). Needs an MPI build.
            The case carries the three families, so this is also the halo's
            per-family gate: measured RED on 2026-09-30 with the exchange stubbed
            to use family 1's schedule for every family (the `family` leg is
            single-rank and stayed green under that stub).
  hybrid    (2 ranks, 2 threads), (3, 4) and (4, 3) against 1x1. Needs both.
  family    the MK+IG+AG run against each family run alone (`--solo K`): the
            families are independent within a step and dt-max binds every cell,
            so family K in company must equal family K alone, value for value.
            RED if the ghost fill skips a family: measured on 2026-09-30 with
            the per-family fill stubbed to family 1's records for every family
            (NaN within the step). The halo's leg is `ranks`.

Every leg first checks the banner: a serial binary reports no thread count, and
`OMP_NUM_THREADS` is then ignored, so the run would compare a file with itself.
Such a build is a SKIP (exit 77), never a PASS.
"""
import glob
import os
import re
import struct
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'test' / 'bench'))
from make_case import build_case  # noqa: E402

ICE_BIN = Path(os.environ.get('ICE_BIN', str(ROOT / 'bin' / 'ICE')))
WORK = Path(os.environ.get('ICE_SCRATCH', str(HERE / 'work')))
MODE = os.environ.get('ICE_EQ_MODE', 'threads')
ITERS = int(os.environ.get('ICE_EQ_ITERS', '0'))          # 0: the crossing count of the case
MASK = os.environ.get('ICE_EQ_MASK_CORNERS', '0') == '1'
EXTRA = os.environ.get('ICE_EQ_EXTRA', '').split()      # extra make_case.py arguments, e.g. --scheme euler
MPIRUN = os.environ.get('MPIRUN', 'mpirun')
SKIP = 77


def say(ok, what, yes, no):
    print('[equiv3d] %s: %s %s' % ('PASS' if ok else 'FAIL', what, yes if ok else no))
    return ok


def run(case, threads, ranks, extra_env=None):
    """Run ICE in case; returns the log text. Raises on a non-zero exit."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), KMP_STACKSIZE='100M', OMP_STACKSIZE='100M')
    env.update(extra_env or {})
    cmd = [str(ICE_BIN)] if ranks == 1 else [MPIRUN, '-np', str(ranks), str(ICE_BIN)]
    with open(str(case / 'log'), 'w') as out, open(str(case / 'err'), 'w') as err:
        rc = subprocess.call(cmd, cwd=str(case), stdout=out, stderr=err, env=env)
    log = (case / 'log').read_text(errors='replace')
    if rc != 0:
        print(log[-2000:])
        raise SystemExit('[equiv3d] ICE failed (exit %d) in %s' % (rc, case))
    return log


def assert_build(log, threads, ranks):
    """The binary must really be the build the leg needs: SKIP otherwise."""
    if threads > 1:
        m = re.search(r'Number of threads\s+-->\s+(\d+)', log)
        if 'Serial execution' in log or m is None:
            print('[equiv3d] SKIP: %s is not an OpenMP build' % ICE_BIN)
            sys.exit(SKIP)
        if int(m.group(1)) != threads:
            raise SystemExit('[equiv3d] FAIL: asked for %d threads, the banner says %s' % (threads, m.group(1)))
    if ranks > 1:
        m = re.search(r'Number of ranks\s+-->\s+(\d+)', log)
        if m is None:
            print('[equiv3d] SKIP: %s is not an MPI build' % ICE_BIN)
            sys.exit(SKIP)
        if int(m.group(1)) != ranks:
            raise SystemExit('[equiv3d] FAIL: asked for %d ranks, the banner says %s' % (ranks, m.group(1)))


# --- raw-appended VTK -------------------------------------------------------

def read_vts(path):
    """({name: bytes} of the CellData arrays, (ni, nj, nk)) of one raw-appended
    .vts, without decoding: byte equality is the comparison."""
    data = Path(path).read_bytes()
    head_end = data.find(b'<AppendedData')
    header = data[:head_end].decode('ascii', 'replace')
    m = re.search(r'header_type="(\w+)"', header)
    hsize = 8 if (m and m.group(1) == 'UInt64') else 4
    hfmt = '<Q' if hsize == 8 else '<I'
    m = re.search(r'WholeExtent="([+-]?\d+) ([+-]?\d+) ([+-]?\d+) ([+-]?\d+) ([+-]?\d+) ([+-]?\d+)"', header)
    e = [int(x) for x in m.groups()]
    dims = (e[1] - e[0], e[3] - e[2], e[5] - e[4])
    start = data.find(b'_', head_end) + 1
    arrays = {}
    cell = header[header.find('<CellData'):header.find('</CellData>')]
    # ORION writes the names doubly quoted (Name=""rho_p1""): tolerate any run of quotes
    for name, off in re.findall(r'<DataArray[^>]*Name="+([^"]+)"+[^>]*offset="(\d+)"', cell):
        off = start + int(off)
        n = struct.unpack(hfmt, data[off:off + hsize])[0]
        arrays[name] = data[off + hsize:off + hsize + n]
    return arrays, dims


def vts_files(case):
    return sorted(glob.glob(str(case / 'OUTPUT' / 'vtk' / '*.vts')))


def check_payload(case):
    """G0: every array holds exactly 8 bytes per cell -- raw doubles, not text."""
    files = vts_files(case)
    if not files:
        raise SystemExit('[equiv3d] FAIL: no .vts written in %s (sol-format = vtk raw?)' % case)
    for f in files:
        arrays, dims = read_vts(f)
        ncells = dims[0] * dims[1] * dims[2]
        if not arrays:
            raise SystemExit('[equiv3d] FAIL: no CellData arrays parsed in %s' % f)
        for name, buf in arrays.items():
            if len(buf) != 8 * ncells:
                raise SystemExit('[equiv3d] FAIL: %s %s holds %d bytes for %d cells' % (f, name, len(buf), ncells))
    return len(files)


def corner_cells(dims):
    """Flat indices (i fastest) of the eight cells with three boundary faces."""
    ni, nj, nk = dims
    return {(i - 1) + ni * ((j - 1) + nj * (k - 1))
            for i in (1, ni) for j in (1, nj) for k in (1, nk)}


def differing_cells(a, b):
    """Per block file: (name, dims, set of differing cell indices) over both runs."""
    fa, fb = vts_files(a), vts_files(b)
    if len(fa) != len(fb) or not fa:
        raise SystemExit('[equiv3d] FAIL: %d and %d block files in %s and %s' % (len(fa), len(fb), a, b))
    out = []
    for x, y in zip(fa, fb):
        ax, dx = read_vts(x)
        ay, dy = read_vts(y)
        if dx != dy or set(ax) != set(ay):
            raise SystemExit('[equiv3d] FAIL: %s and %s do not describe the same block' % (x, y))
        diff = set()
        for name in ax:
            p, q = ax[name], ay[name]
            if p != q:
                diff |= {c for c in range(len(p) // 8) if p[8 * c:8 * c + 8] != q[8 * c:8 * c + 8]}
        out.append((Path(x).name, dx, diff))
    return out


def compare(a, b, what):
    """Byte comparison; with the corner mask, cells with three boundary faces
    are left out and the verdict says how many cells differed where."""
    blocks = differing_cells(a, b)
    ndiff = sum(len(d) for _, _, d in blocks)
    ncorner = sum(len(d & corner_cells(dims)) for _, dims, d in blocks)
    if MASK:
        rest = ndiff - ncorner
        return say(rest == 0, what,
                   'bit-identical outside the %d masked corner cells (%d of them differ)'
                   % (8 * len(blocks), ncorner),
                   '%d cells differ outside the corners (%d corner cells too)' % (rest, ncorner))
    return say(ndiff == 0, what, 'bit-identical',
               '%d cells differ, %d of them block corners with three boundary records' % (ndiff, ncorner))


def family_arrays(case, k):
    """{basename: [bytes per block]} of family k's arrays, in block order."""
    out = {}
    for f in vts_files(case):
        arrays, _ = read_vts(f)
        for name, buf in arrays.items():
            if name.endswith('_p%d' % k):
                out.setdefault(name[:-len('_p%d' % k)], []).append(buf)
    return out


# --- legs -------------------------------------------------------------------

def make(name, extra=()):
    case = WORK / name
    args = [str(case), '--preset', 'small'] + list(extra)
    if ITERS:
        args += ['--iters', str(ITERS)]
    build_case(args + EXTRA)
    return case


def leg_threads():
    ref = make('t1')
    run(ref, 1, 1)
    print('[equiv3d] %d blocks written, raw payload checked; %s'
          % (check_payload(ref), 'corner cells masked' if MASK else 'every cell compared'))
    ok = True
    for t in (2, 3, 4, 8):
        case = make('t%d' % t)
        assert_build(run(case, t, 1), t, 1)
        ok &= compare(ref, case, '1 thread vs %d threads:' % t)
    return ok


def leg_groups():
    ref = make('g1')
    run(ref, 1, 1, {'ICE_THREAD_GROUPS': '1'})
    print('[equiv3d] %d blocks written, raw payload checked; every cell compared' % check_payload(ref))
    ok = True
    for groups, t in ((2, 4), (3, 6), (4, 8), (2, 3)):
        case = make('g%d-t%d' % (groups, t))
        log = run(case, t, 1, {'ICE_THREAD_GROUPS': str(groups)})
        assert_build(log, t, 1)
        m = re.search(r'thread groups (\d+) \(asked', log)
        if m is None or int(m.group(1)) != groups:
            raise SystemExit('[equiv3d] FAIL: asked for %d thread groups, the banner says %s'
                             % (groups, m.group(1) if m else 'nothing'))
        ok &= compare(ref, case, '1 thread vs %d threads in %d groups:' % (t, groups))
    return ok


def leg_ranks():
    ref = make('r1')
    run(ref, 1, 1)
    check_payload(ref)
    ok = True
    for r in (2, 3, 4, 5):
        case = make('r%d' % r)
        assert_build(run(case, 1, r), 1, r)
        ok &= compare(ref, case, '1 rank vs %d ranks:' % r)
    # more ranks than blocks: two blocks on three ranks, one of them idle
    ref2 = make('r1-2blk', ['--blocks', '2', '1', '1'])
    run(ref2, 1, 1)
    case = make('r3-2blk', ['--blocks', '2', '1', '1'])
    assert_build(run(case, 1, 3), 1, 3)
    ok &= compare(ref2, case, '2 blocks on 1 rank vs 3 ranks (one idle):')
    return ok


def leg_hybrid():
    ref = make('h1x1')
    run(ref, 1, 1)
    check_payload(ref)
    ok = True
    for r, t in ((2, 2), (3, 4), (4, 3)):
        case = make('h%dx%d' % (r, t))
        assert_build(run(case, t, r), t, r)
        ok &= compare(ref, case, '1x1 vs %d ranks x %d threads:' % (r, t))
    return ok


def leg_family():
    together = make('fam-all')
    run(together, 1, 1)
    check_payload(together)
    ok = True
    for k in (1, 2, 3):
        alone = make('fam-solo%d' % k, ['--solo', str(k)])
        run(alone, 1, 1)
        got, want = family_arrays(together, k), family_arrays(alone, 1)
        same = set(got) == set(want) and all(got[n] == want[n] for n in want)
        ok &= say(same, 'family %d alone vs in company (MK+IG+AG):' % k,
                  'equal value for value', 'differ')
    return ok


def main():
    if not ICE_BIN.exists():
        raise SystemExit('[equiv3d] no ICE binary at %s -- build first' % ICE_BIN)
    WORK.mkdir(parents=True, exist_ok=True)
    legs = {'threads': leg_threads, 'groups': leg_groups, 'ranks': leg_ranks, 'hybrid': leg_hybrid, 'family': leg_family}
    if MODE not in legs:
        raise SystemExit('[equiv3d] unknown ICE_EQ_MODE %s' % MODE)
    ok = legs[MODE]()
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
