"""O. Solidification in a steady stream - IGLOO's solid-box as an Eulerian stream, and melting in a heated one.

Case N's material in a 1-D steady stream (box 0.15 m, nx cells, ny = 1): a code-402 inlet on face 1 (29.5 kg/(s m^2),
10 m/s, radius 15 um), extrapolation elsewhere, drag = NoDrag, heat-transfer = Stokes, local time stepping, the IC the
inlet state everywhere (every cell approaches its steady state from above) unless stated. At steady state a first-order
upwind solution is a recursion from the inlet, a = (A/m)(dx/u): the liquid e = (e_up + a T_g)/(1 + a/c_l) while its
temperature stays above T_nuc, otherwise the equilibrium state (the solid, the plateau e_up + a (T_g - T_m), or the
liquid); chi is the inflow's between T_nuc and T_m, 1 below, 0 above. It is ICE's discrete solution to rounding.

  O0  every cell liquid (T-melt 200 K): the solidifying run's file holds the plain run's blocks byte for byte, and
      two blocks of zeros for f and chi. MK/Saurel and IG/Rusanov, MUSCL.
  O1  first order, nx = 120, against the recursion per cell, with its energy balance.
  O2  MUSCL at nx = 120, 240, 480 against the closed form: the front within a cell, T = T_m on the plateau, the
      liquid error falling with dx - at first order: ICE holds the inlet state at the ghost-cell centre, half a cell
      upstream of the face, which shifts the profile by dx/2, (T_0 - T_g) dx/(2 L_l) = 10.6 K at nx = 120 - and the
      exit within the solid's first-order bound. Every nx must be a fixed point: the family reconstructs
      T_liq = e/c_l, continuous through the recalescence, so the front does not depend on its own cell's phase
      (reconstructing T instead, nx = 480 cycles with period 2, the front cell's chi 0.6/1.0 and an f ripple of 2.9e-3).
  O3  IG/Rusanov, MUSCL nx = 240: the front of O2 within a cell and chi < 1/2 in every cell upstream of it. The
      Rusanov leak is c/(2u + c) = 5e-5 at the inlet's P = 1e-6, and with T_liq reconstructed nothing else reaches
      upstream (reconstructing T, the cell upstream of the front held 0.2).
  O3b IG/Rusanov above the no-flip bound: O1's stream scaled to u = 0.1 m/s (box 1.5e-3 m), where the inlet's
      P = 1e-6 gives c = sqrt(3P/rho_p) = 2.5u and 4u; first order against the driver's own Rusanov model.
  O4  a solid injected at 1500 K into a 3000 K gas, h-fus = 4e5: solid, melting at T_m, liquid - the melting law;
      first order, nx = 120 and 240, MK/Saurel and IG/Rusanov, from an empty domain.
  O5  a steady run with the gas at 300 K, then a restart in place with the gas at 1000 K: the front moves back to
      the 1000 K run's cell (no lock), per cell equal to a fresh 1000 K run.
  O6  O2's stream at nx = 120 on two mesh blocks joined by a connection after cell 40 (on the freezing plateau, so
      the partner cells carry f and chi), at first order with RK2 and at MUSCL with RK3: the one-block run, the
      two-block run and, with ICE_MPI_RANKS = n > 1, the two-block run on n MPI ranks are fixed points. At first order
      the two-block run equals the one-block run value for value: the connection hands each block its partner's
      cells, a zero slope keeps the cell lengths out of the faces, and a cell's residual sums the same two fluxes. With
      MUSCL the connection face is reconstructed with compute_bound's lengths, which differ from the interior's in the
      last bits (the pre-existing dl_g2 ~ dl_g1 of the ledger), so they agree to 1e-13 in T (relative) and 1e-14 in f
      and chi (measured 8.6e-15 and 1.0e-15). The MPI run equals the serial two-block run value for value.

Every steady run is compared with a twin at 40 nx + 1 iterations - an odd gap, not a multiple of 3, so a cycle of
period 2 or 3 shows - in T, f and chi, cell by cell, to 1e-12 (else it is not a fixed point).
"""
import importlib.util
import math
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import (Case, Report, inlet_payload, read_solution, run_ice, write_tec_blocks,  # noqa: E402
                    ICE_BIN)

_spec = importlib.util.spec_from_file_location('N_solid', str(HERE.parent / 'N-solid-cooling' / 'run.py'))
N = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(N)

WORK = HERE / 'work'
LX, U_IN, T_IN, GP = 0.15, 10.0, 2400.0, 29.5
RHO_IN = GP / U_IN
RHO_E = 1.0e-20                            # the content of an empty domain


# ---------------------------------------------------------------------------
#  The discrete oracles
# ---------------------------------------------------------------------------

def jump(Tl, h):
    f0 = N.CL * (N.TM - Tl) / h
    if f0 < 1.0:
        return N.TM, f0
    return N.TM - (N.CL * (N.TM - Tl) - h) / N.CS, 1.0


def state(e, chi, h):
    """solid_state: (T, f, chi') of energy e with the transported chi."""
    Tl = e / N.CL
    if Tl >= N.TM:
        return Tl, 0.0, 0.0
    if Tl <= N.TN:
        T, f = jump(Tl, h)
        return T, f, 1.0
    if chi >= 0.5:
        T, f = jump(Tl, h)
        return T, f, chi
    return Tl, 0.0, chi


def energy(T, f, h):
    return N.CL * T - f * (h - (N.CL - N.CS) * (N.TM - T))


def recursion(nx, Tg, T0, h=N.HF, chi0=0.0, lx=LX, u=U_IN):
    """First-order upwind steady state, cell by cell from the inlet: [(T, f, chi, e)]."""
    a = (N.A / N.M) * (lx / nx) / u
    es = N.CL * N.TM - h
    e_up = energy(T0, 1.0 if chi0 >= 0.5 else 0.0, h)
    chi_up, out = chi0, []
    for _ in range(nx):
        if chi_up >= 0.5:
            e = (e_up + a * (Tg - N.TM + es / N.CS)) / (1.0 + a / N.CS)
            if e >= es:
                e = e_up + a * (Tg - N.TM)
                if e >= N.CL * N.TM:
                    e = (e_up + a * Tg) / (1.0 + a / N.CL)
        else:
            e = (e_up + a * Tg) / (1.0 + a / N.CL)
            if e / N.CL <= N.TN:
                e = (e_up + a * (Tg - N.TM + es / N.CS)) / (1.0 + a / N.CS)
                if e >= es:
                    e = e_up + a * (Tg - N.TM)
        T, f, chi = state(e, chi_up, h)
        out.append((T, f, chi, e))
        e_up, chi_up = e, chi
    return out


def windows(nx, Tg, T0, h=N.HF):
    """Cells whose liquid inflow admits both a liquid and a nucleated steady state (plan section 1.4 (e))."""
    a = (N.A / N.M) * (LX / nx) / U_IN
    en = N.CL * N.TN
    rec = recursion(nx, Tg, T0, h)
    ups = [N.CL * T0] + [r[3] for r in rec[:-1]]
    return [i + 1 for i in range(nx) if (i == 0 or rec[i - 1][2] < 0.5)
            and en + a * (N.TN - Tg) < ups[i] <= en + a * (N.TM - Tg)]


def rusanov_model(nx, Tg, T0, ratio, h=N.HF, cfl=0.4, tol=1.0e-9, itmax=400000):
    """The steady first-order Rusanov solution of e and chi (uniform rho and u, signal speed c = ratio u)."""
    sig = cfl / (1.0 + ratio)
    k = sig * (LX / nx / U_IN) * N.A / N.M
    F = lambda l, r: 0.5 * ((2.0 + ratio) * l - ratio * r)                       # noqa: E731
    e_in, q_in = N.CL * T0, 0.0
    e, q = [e_in] * nx, [0.0] * nx
    T = [0.0] * nx
    for i in range(nx):
        T[i], _, q[i] = state(e[i], q[i], h)
    for _ in range(itmax):
        ne, nq = [0.0] * nx, [0.0] * nx
        for i in range(nx):
            el, ql = (e_in, q_in) if i == 0 else (e[i - 1], q[i - 1])
            er, qr = (e[i], q[i]) if i == nx - 1 else (e[i + 1], q[i + 1])
            ne[i] = e[i] + sig * (F(el, e[i]) - F(e[i], er)) + k * (Tg - T[i])
            nq[i] = q[i] + sig * (F(ql, q[i]) - F(q[i], qr))
        dm = max(abs(ne[i] - e[i]) for i in range(nx))
        out = [state(ne[i], nq[i], h) for i in range(nx)]
        T = [o[0] for o in out]
        e, q = ne, [o[2] for o in out]
        if dm < tol:
            break
    return [(out[i][0], out[i][1], out[i][2], e[i]) for i in range(nx)]


def closed_form(x, Tg=600.0, T0=T_IN, h=N.HF):
    """The cooling stream at x = u t (N's closed form in time)."""
    return N.cooling(x / U_IN, T0, Tg, h)[0]


# ---------------------------------------------------------------------------
#  Runs
# ---------------------------------------------------------------------------

def stream(label, nx, Tg, phase=N.PHASE, closure='MK', order='first-order', T0=T_IN, ic=None, lx=LX, u=U_IN,
           gp=GP, iters=None, cfl=0.8, rk='RK2'):
    """The steady stream; ic = (rho, T) of the initial content (the inlet state by default)."""
    case = Case(WORK / label.lower(), nx=nx, Lx=lx)
    rho0, Tic = ic if ic is not None else (gp / u, T0)
    base = N.LAYOUT[closure]
    # the pseudo-pressure scales with the content, so an empty IC keeps the inflow's small signal speed
    Pic = N.PSEUDO * rho0 / (gp / u)
    fields = [rho0, u, 0.0, 0.0] + dict(IG=[Pic], AG=[Pic, 0.0, 0.0, Pic, 0.0, Pic]).get(closure, []) \
        + [Tic, rho0 / (N.RHO_M * 4.0 / 3.0 * math.pi * N.RP ** 3)]
    case.write_ic(['%s1' % s for s in base], fields)
    case.gas(rho=N.RHO_G, u=u, v=0.0, w=0.0, T=Tg, R=N.R_G, gam=N.GAM, k=N.KG, mu=N.MU)
    case.boundaries('extrapolation', inlets=[inlet_payload(gp, u, T0, N.RP)])
    case.phase(phase)
    limiter = 'vanleer' if order == 'MUSCL' else 'none'
    case.ini(iters=iters or 20 * nx, cfl=cfl, rk=rk, drag='NoDrag', heat='Stokes', rho_al=N.RHO_M, cs=N.CL,
             reconstruction=order, limiter=limiter, numerics={'time-accurate': False}, closures=(closure,))
    return case


def steady(rep, label, builder, **kw):
    """A steady run at 20·nx iterations and its twin at 40·nx + 1, an odd gap and not a multiple of 3, so that a cycle of
    period 2 or 3 shows; T (relative), f and chi compared cell by cell to 1e-12. The second run's solution, or None (FAIL
    recorded) if a run failed or the two differ."""
    nx = kw['nx']
    base = kw.pop('iters', None) or 20 * nx
    one = N.attempt(rep, label, builder(label, iters=base, **kw))
    if one is None:
        return None
    two = N.attempt(rep, label + ' (twin)', builder(label + '-twin', iters=2 * base + 1, **kw))
    if two is None:
        return None
    closure = kw.get('closure', 'MK')
    a, b = cells(one, closure), cells(two, closure)
    dT = max(abs(x[0] - y[0]) / abs(y[0]) for x, y in zip(a, b))
    df, dchi = (max(abs(x[k] - y[k]) for x, y in zip(a, b)) for k in (1, 2))
    rep.check(max(dT, df, dchi) <= 1.0e-12, '%s: a fixed point (T, f, chi move %.1e, %.1e, %.1e from %d to %d iterations)'
              % (label, dT, df, dchi, base, 2 * base + 1))
    return two


def stream2(label, nx, Tg, split, ranks=1, **kw):
    """stream() on two mesh blocks joined by a connection after cell `split` (block 1: cells 1..split), run on `ranks`
    MPI ranks; the case's run() returns both zones as one solution, cells in x order."""
    case = stream(label, nx, Tg, **kw)
    assert case.ny == 1 and 0 < split < nx
    xs = (case.xn[:split + 1], case.xn[split:])
    for name, zone in (('part-ic.tec', 'CD'), ('gas.tec', 'GAS')):
        one = read_solution(case.dir / 'INPUT' / name)
        write_tec_blocks(case.dir / 'INPUT' / name, one['names'],
                         [(xs[0], case.yn, case.zn, [v[:split] for v in one['var']], 'B1-' + zone),
                          (xs[1], case.yn, case.zn, [v[split:] for v in one['var']], 'B2-' + zone)])
    inlet = (case.dir / 'INPUT/part-bc.txt').read_text().splitlines()[1]
    rows = []
    for b, n in ((1, split), (2, nx - split)):
        def rec(i, f, code, b=b):
            return '%8d%8d%8d%8d%8d%8d' % (b, i, 1, 1, f, code)
        if b == 1:
            rows += [rec(1, 1, 402), inlet, rec(n, 2, 101), '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (2, 1, 1, 1, 1, 1, 0, 0, 1)]
        else:
            rows += [rec(1, 1, 101), '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (1, split, 1, 1, 2, 1, 0, 0, 1), rec(n, 2, 400)]
        rows += [rec(i, f, 400) for f in (3, 4) for i in range(1, n + 1)]
        rows += [rec(i, f, 0) for f in (5, 6) for i in range(1, n + 1)]
    (case.dir / 'INPUT/part-bc.txt').write_text('\n'.join(rows) + '\n')
    case.run = lambda threads=1: run_blocks(case, ranks)
    return case


def run_blocks(case, ranks):
    """Run a two-block case on `ranks` MPI ranks (serially when 1); both zones of its solution as one."""
    if ranks > 1:
        env = dict(os.environ, OMP_NUM_THREADS='1', KMP_STACKSIZE='100M')
        with open(str(case.dir / 'log'), 'w') as out, open(str(case.dir / 'err'), 'w') as err:
            rc = subprocess.call(['mpirun', '-np', str(ranks), str(ICE_BIN)], cwd=str(case.dir), stdout=out,
                                 stderr=err, env=env)
        if rc != 0:
            raise RuntimeError('ICE failed (exit %d) in %s - see log and err there' % (rc, case.dir))
    else:
        run_ice(case.dir)
    zones = [read_solution(case.dir / 'OUTPUT/part-field.tec', zone=z) for z in (0, 1)]
    return {'time': zones[0]['time'], 'nx': zones[0]['nx'] + zones[1]['nx'], 'ny': 1, 'names': zones[0]['names'],
            'var': [a + b for a, b in zip(zones[0]['var'], zones[1]['var'])]}


def cells(sol, closure='MK'):
    c = N.columns(closure)
    return [(sol['var'][c['T']][i], sol['var'][c['f']][i], sol['var'][c['chi']][i]) for i in range(sol['nx'])]


def first_nucleated(cs):
    for i, (T, f, chi) in enumerate(cs):
        if f > 0.0:
            return i + 1
    return None


def against_recursion(rep, label, cs, rec, tolT, tolf=1.0e-9):
    wT = max(abs(c[0] - r[0]) / r[0] for c, r in zip(cs, rec))
    wf = max(abs(c[1] - r[1]) for c, r in zip(cs, rec))
    wchi = max(abs(c[2] - r[2]) for c, r in zip(cs, rec))
    print('       worst per-cell difference: T %.2e relative, f %.2e, chi %.2e' % (wT, wf, wchi))
    rep.check(wT <= tolT and wf <= tolf and wchi <= tolf,
              '%s: every cell on the recursion (T to %.0e, f and chi to %.0e)' % (label, tolT, tolf))


def o0(rep):
    for closure in ('MK', 'IG'):
        label = 'O0-%s' % closure
        print('   %s every cell liquid (T-melt 200 K), MUSCL nx = 120, against the same case without tokens' % label)
        runs = {}
        for tag, phase in (('plain', 'A 1'),
                           ('solid', 'A 1 solidification=on T-melt=200 T-nuc=100 h-fus=1.07e6 cp-solid=600')):
            case = stream('%s-%s' % (label, tag), 120, 600.0, phase=phase, closure=closure, order='MUSCL')
            if N.attempt(rep, '%s %s' % (label, tag), case) is None:
                return
            runs[tag] = (case.dir / 'OUTPUT/part-field.tec').read_text()
        nb, nc = len(N.LAYOUT[closure]), 120
        def tokens(text):
            body = text.split('ZONE', 1)[1].split('\n', 1)[1]
            return body.split()
        tp, ts = tokens(runs['plain']), tokens(runs['solid'])
        nn = 3 * 121 * 2 * 2
        same = tp[:nn + nb * nc] == ts[:nn + nb * nc]
        zeros = len(ts) == nn + (nb + 2) * nc and all(float(v) == 0.0 for v in ts[nn + nb * nc:])
        rep.check(same and zeros, '%s: the plain blocks byte for byte, f and chi all zero' % label)


def o1(rep):
    print('   O1  first order, nx = 120, against the recursion')
    win = windows(120, 600.0, T_IN)
    rep.check(not win, 'O1 setup: no window cell at nx = 120 (%s)' % (win or 'none'))
    sol = steady(rep, 'O1', stream, nx=120, Tg=600.0)
    if sol is None:
        return None
    cs, rec = cells(sol), recursion(120, 600.0, T_IN)
    fr, frr = first_nucleated(cs), first_nucleated([(r[0], r[1], r[2]) for r in rec])
    print('       front cell %s (recursion %s), exit %.4f K (recursion %.4f K)' % (fr, frr, cs[-1][0], rec[-1][0]))
    against_recursion(rep, 'O1', cs, rec, 1.0e-9)
    a = (N.A / N.M) * (LX / 120) / U_IN
    e_out = energy(cs[-1][0], cs[-1][1], N.HF)
    lost = sum(a * (T - 600.0) for T, f, chi in cs)
    err = abs((N.CL * T_IN - e_out) - lost) / lost
    rep.check(err <= 1.0e-10, 'O1: e_in - e_out = the heat given to the gas, to %.1e' % err)
    return fr


def o2(rep):
    Ll = U_IN * N.TAU_L
    xn = Ll * math.log((T_IN - 600.0) / (N.TN - 600.0))
    f0 = N.CL * (N.TM - N.TN) / N.HF
    xs = xn + U_IN * (1.0 - f0) * N.M * N.HF / (N.A * (N.TM - 600.0))
    Ls = U_IN * N.TAU_S
    errs, fronts = [], {}
    for nx in (120, 240, 480):
        print('   O2  MUSCL, nx = %d' % nx)
        sol = steady(rep, 'O2-%d' % nx, stream, nx=nx, Tg=600.0, order='MUSCL')
        if sol is None:
            return fronts
        dx = LX / nx
        cs = cells(sol)
        fr = first_nucleated(cs)
        fronts[nx] = fr
        home = int(xn // dx) + 1
        rep.check(fr is not None and abs(fr - home) <= 1, 'O2-%d: front cell %s, x_n = %.6f m is in cell %d' % (nx, fr, xn, home))
        plateau = [i for i in range(nx) if i * dx > xn + dx and (i + 1) * dx < xs - dx]
        rep.check(all(cs[i][0] == N.TM and cs[i][2] == 1.0 for i in plateau),
                  'O2-%d: T = T_m exactly and chi = 1 in the %d plateau cells' % (nx, len(plateau)))
        liquid = [i for i in range(nx) if (i + 1) * dx < xn - 2 * dx]
        errs.append(max(abs(cs[i][0] - closed_form((i + 0.5) * dx)) for i in liquid))
        xc = (nx - 0.5) * dx
        Tcf = closed_form(xc)
        bound = (Tcf - 600.0) / Ls * 2 * dx
        rep.check(abs(cs[-1][0] - Tcf) <= bound, 'O2-%d: exit %.4f K, closed form %.4f K, within %.1f K'
                  % (nx, cs[-1][0], Tcf, bound))
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    print('       liquid errors %s, ratios %s' % (', '.join('%.3e' % e for e in errs), ', '.join('%.2f' % r for r in ratios)))
    rep.check(all(r >= 1.8 for r in ratios), 'O2: the liquid error falls at first order (at least 1.8x per refinement)')
    return fronts


def o3(rep, front_mk):
    print('   O3  IG/Rusanov, MUSCL nx = 240')
    sol = steady(rep, 'O3', stream, nx=240, Tg=600.0, closure='IG', order='MUSCL')
    if sol is None:
        return
    cs = cells(sol, 'IG')
    fr = first_nucleated(cs)
    up = max([c[2] for c in cs[:fr - 1]] or [0.0]) if fr else float('nan')
    print('       front cell %s (MK %s), max chi upstream %.3e' % (fr, front_mk, up))
    rep.check(fr is not None and front_mk is not None and abs(fr - front_mk) <= 1 and up < 0.5,
              'O3: front within a cell of MK\'s, chi < 1/2 upstream of it')


def o3b(rep):
    for ratio, gp, front in ((2.5, 4.8e-6, 30), (4.0, 1.875e-6, 29)):
        rho = gp / 0.1
        c = math.sqrt(3.0 * N.PSEUDO / rho)
        print('   O3b IG/Rusanov first order, u = 0.1 m/s, rho_p = %.4e kg/m^3: c = %.4f m/s = %.2f u' % (rho, c, c / 0.1))
        sol = steady(rep, 'O3b-%g' % ratio, stream, nx=120, Tg=600.0, closure='IG', lx=LX * 0.01, u=0.1, gp=gp,
                     iters=60 * 120)
        if sol is None:
            continue
        cs = cells(sol, 'IG')
        model = rusanov_model(120, 600.0, T_IN, c / 0.1)
        fr, frm = first_nucleated(cs), first_nucleated([(m[0], m[1], m[2]) for m in model])
        up = max([x[2] for x in cs[:fr - 1]] or [0.0]) if fr else float('nan')
        print('       front cell %s (model %s, upwind 31), exit %.4f K (model %.4f K), max chi upstream %.4f'
              % (fr, frm, cs[-1][0], model[-1][0], up))
        against_recursion(rep, 'O3b c = %gu' % ratio, cs, model, 1.0e-6, 1.0e-6)
        rep.check(fr == front and up < 0.5 and 31 - fr <= 2,
                  'O3b c = %gu: front at cell %s (expected %d, at most two cells above the upwind 31)' % (ratio, fr, front))


def o4(rep):
    phase = 'A 1 solidification=on h-fus=4e5 cp-solid=600'
    for closure in ('MK', 'IG'):
        for nx in (120, 240):
            label = 'O4-%s-%d' % (closure, nx)
            print('   %s a solid at 1500 K into a 3000 K gas, h-fus = 4e5' % label)
            sol = steady(rep, label, stream, nx=nx, Tg=3000.0, phase=phase, closure=closure, T0=1500.0,
                         ic=(RHO_E, 1500.0))
            if sol is None:
                continue
            cs = cells(sol, closure)
            rec = recursion(nx, 3000.0, 1500.0, h=4.0e5, chi0=1.0)
            mush = next((i + 1 for i, c in enumerate(cs) if c[1] < 1.0), None)
            liq = next((i + 1 for i, c in enumerate(cs) if c[2] < 0.5), None)
            print('       first mush cell %s, first liquid cell %s, exit %.4f K (recursion %.4f K)'
                  % (mush, liq, cs[-1][0], rec[-1][0]))
            against_recursion(rep, label, cs, rec, 1.0e-9 if closure == 'MK' else 1.0e-4,
                              1.0e-9 if closure == 'MK' else 1.0e-4)
            rep.check(all(c[0] <= N.TM for c in cs if c[1] > 0.0), '%s: no cell with f > 0 above T_m' % label)


def o5(rep):
    print('   O5  gas at 300 K, then a restart in place at 1000 K (cfl 0.4)')
    for Tg in (300.0, 1000.0):
        win = windows(120, Tg, T_IN)
        rep.check(not win, 'O5 setup: no window cell at %g K (%s)' % (Tg, win or 'none'))
    case = stream('O5', 120, 300.0, cfl=0.4)
    sol1 = N.attempt(rep, 'O5 run 1', case)
    if sol1 is None:
        return
    print('       run 1: front cell %s, exit %.4f K' % (first_nucleated(cells(sol1)), cells(sol1)[-1][0]))
    N.restart_in_place(case, {'newrun': False},
                       gas=dict(rho=N.RHO_G, u=U_IN, v=0.0, w=0.0, T=1000.0, R=N.R_G, gam=N.GAM, k=N.KG, mu=N.MU))
    try:
        run_ice(case.dir)
    except RuntimeError as exc:
        rep.check(False, 'O5 restart: %s' % exc)
        return
    cs = cells(read_solution(case.dir / 'OUTPUT/part-field.tec'))
    ctrl = steady(rep, 'O5-control', stream, nx=120, Tg=1000.0, cfl=0.4)
    rec = recursion(120, 1000.0, T_IN)
    print('       restart: front cell %s, exit %.4f K (fresh 1000 K recursion %s, %.4f K)'
          % (first_nucleated(cs), cs[-1][0], first_nucleated([(r[0], r[1], r[2]) for r in rec]), rec[-1][0]))
    against_recursion(rep, 'O5 restart', cs, rec, 1.0e-9)
    if ctrl is not None:
        cc = cells(ctrl)
        w = max(abs(a[0] - b[0]) / b[0] for a, b in zip(cs, cc))
        rep.check(w <= 1.0e-9, 'O5: the restart equals the fresh 1000 K run to %.1e' % w)


def o6(rep):
    ranks = int(os.environ.get('ICE_MPI_RANKS', '1'))
    for order, rk in (('first-order', 'RK2'), ('MUSCL', 'RK3')):
        label = 'O6-%s' % ('first' if order == 'first-order' else 'muscl')
        print("   %s O2's stream, nx = 120, %s with %s: one block, two blocks joined after cell 40%s"
              % (label, order, rk, (', and those on %d MPI ranks' % ranks) if ranks > 1 else ''))
        one = steady(rep, label + '-1', stream, nx=120, Tg=600.0, order=order, rk=rk)
        runs = [('two blocks', steady(rep, label + '-2', stream2, nx=120, Tg=600.0, order=order, rk=rk, split=40))]
        if ranks > 1:
            sol = steady(rep, label + '-mpi', stream2, nx=120, Tg=600.0, order=order, rk=rk, split=40, ranks=ranks)
            if sol is not None:
                witness = (WORK / (label + '-mpi-twin').lower() / 'log').read_text(errors='replace')
                rep.check('Number of ranks   --> %4d' % ranks in witness and 'MPI partition: 2 blocks over %d ranks'
                          % ranks in witness, '%s: the solver reported %d MPI ranks and the two blocks spread over them '
                          '(a serial binary under mpirun would not)' % (label, ranks))
            runs.append(('two blocks on %d ranks' % ranks, sol))
        if one is None:
            continue
        a = cells(one)
        print('       one block: front cell %s, exit %.4f K' % (first_nucleated(a), a[-1][0]))
        tol = (0.0, 0.0) if order == 'first-order' else (1.0e-13, 1.0e-14)
        ref = [('the one-block run', a, tol)]
        for n, (what, sol) in enumerate(runs):
            if sol is None:
                continue
            b = cells(sol)
            if n == 0:
                ref.append(('the serial two-block run', b, (0.0, 0.0)))
            other, c, (tT, tf) = ref[0] if n == 0 else ref[-1]
            dT = max(abs(x[0] - y[0]) / y[0] for x, y in zip(b, c))
            df, dchi = (max(abs(x[k] - y[k]) for x, y in zip(b, c)) for k in (1, 2))
            rep.check(len(b) == len(c) and dT <= tT and max(df, dchi) <= tf,
                      '%s: %s agree with %s (T %.1e relative, f %.1e, chi %.1e; %s)'
                      % (label, what, other, dT, df, dchi,
                         'value for value' if tT == 0.0 else 'within %.0e and %.0e' % (tT, tf)))


def main():
    rep = Report('O. Solidification, steady stream')
    for leg, fn in (('O0', o0), ('O1', o1)):
        if N.selected(leg):
            fn(rep)
    fronts = o2(rep) if N.selected('O2') or N.selected('O3') else {}
    if N.selected('O3'):
        o3(rep, fronts.get(240))
    for leg, fn in (('O3b', o3b), ('O4', o4), ('O5', o5), ('O6', o6)):
        if N.selected(leg):
            fn(rep)
    rep.close(WORK)


if __name__ == '__main__':
    main()
