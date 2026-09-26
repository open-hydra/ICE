"""L. Several materials - each family on its own material, tokens, inlet records per family.

INPUT/part-phase.txt names the materials, one line each, "<name> <groups> [key=value ...]",
and the families map onto them in that order (material-major). The property table gives
one zone per material. Every leg below runs two or three MK families side by side in one
case, so a family that reads another family's material, model or inlet is seen directly.

  L1   heating. Two families with the same state (rho_p, n) at rest, T_p = 300 K, in a gas
       at rest at 400 K; heat-transfer = Stokes (Nu = 2), drag = NoDrag, emissivity 0.
       Materials A (cp 1000, rho 2000) and B (cp 2000, rho 1000), phase "A 1" / "B 1".
       ICE recovers Rp = (3 rho_p / (4 pi n rho_mat))^(1/3) from the state, so each family
       relaxes as T_g + (T_0 - T_g) exp(-t/tau) with its own
           tau = rho_p cp / (4 pi k Rp n),
       tau_B / tau_A = 2 (1/2)^(1/3) = 1.587. At t = tau_A a family that heats with the
       other material is off by 0.165 of the initial gap.
  L1b  the same with "A 2" / "B 1" and three families: families 1 and 2 are A, 3 is B.
  L1c  control: "A 2", one zone (A), two families; every binary passes it.
  L2   tokens. "A 1 evaporation=CEM" / "B 1", [ICE-Physics] evaporation absent (default
       none): family A follows the CEM d-squared slope of case F (frozen T_p, cp 1e12),
       family B keeps rho_p and n to the 15 digits the solution file prints. L2b: the same case with evaporation = none
       written out gives the same part-field.tec byte for byte. L2d: evaporation = CEM in
       the INI and "B 1 evaporation=none": A evaporates, B does not.
  L2c  the token alone: one material, "A 1 evaporation=CEM", no table and no INI model.
  L2e  control: one material, no token, evaporation = CEM in the INI.
  L2f  boiling temperature per material. "A 1 evaporation=CEM" / "B 1 evaporation=CEM", two
       zones of the same material, boiling-temperature = 373.15 403.15. L2g: the same values
       under the alias Tboil give L2f's part-field.tec byte for byte. L2h: control, the two
       values swapped change it, so each value reaches its own material.
  L3   inlet records per family. One material, "A 2", an empty domain, code-402 inlets
       on face 1 (mass flux 0.01 kg/(s m^2), 10 m/s, 300 K) in ATLAS's order: the records
       of family 1 (radius 10 um), then those of family 2 (20 um). No coupling, so both
       families move identically and only n differs: n_1/n_2 = (r_2/r_1)^3 = 8 in the
       cells the inlet has filled. L3c: control, both radii 10 um, ratio 1.
  L3b  two mesh blocks (20 and 12 cells), ATLAS's order: mesh block, then family, then
       faces. Block 1 has radii 10/20 um, block 2 20/10 um: ratio 8 and 1/8.
       L3d: control, the same two blocks with one record per face (fanned out),
       10 um in block 1 and 20 um in block 2 for both families: ratio 1 in each.
  L4   inlet density. "A 1" / "B 1", rho 2000 and 1000 in the table, no INI density, one
       inlet record per face for both families (10 um): each family's n_p is its
       rho_p / (rho_mat 4/3 pi r^3) with its own rho_mat, so n_1/n_2 = 1/2.
       L4c: control, "A 2" with one zone (rho 2000).

Tolerances.
  L1: RK2 on y' = -y/tau is exact to (dt/tau)^3/6 per step, so after t = tau with
      dt = tau_A/1000 the error is (t/tau)(dt/tau)^2/6 |T - T_g| <= 6e-8 of the initial
      gap; the tolerance, 1e-6 of the gap, is 16 times that and 1.6e5 below the signal.
  L2: case F's 1e-6 on the d-squared slope, which F's own CEM leg meets at 1.6e-8 (the
      same arithmetic: a constant table cp returns the INI value's number). B: exact.
      L2g, L2h: byte comparisons of two runs of one binary.
  L3, L4: r_2 = 2 r_1 and rho_A = 2 rho_B are exact in binary, the flux and the limiter
      are homogeneous of degree 1 in n, and the inlet formula divides by rho_mat r^3, so
      n_1/n_2 is exact up to the initial content (1e-20 kg/m^3 against 1e-3 at the inlet).
      The per-family value n rho_mat 4/3 pi r^3 / rho_p = 1 carries a few roundings.
      1e-12 bounds both. A cell counts as filled at rho_p >= rho_in / 2.
"""
import importlib.util
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import (Case, Report, bc_rows, inlet_payload, linspace, read_solution,   # noqa: E402
                    run_ice, try_run, write_tec_blocks)

# Case F's evaporation oracle, loaded under its own name
_spec = importlib.util.spec_from_file_location('F_evaporation', str(HERE.parent / 'F-evaporation' / 'run.py'))
F = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(F)

WORK = HERE / 'work'
RHO, U, V, W, T, N = range(6)
NV = 6                                   # MK variables per family

# --- L1: heating ------------------------------------------------------------------
TG1, TP1 = 400.0, 300.0
RHO_G1, R1, GAM1, K1, MU1 = 1.2, 287.05, 1.4, 0.026, 1.8e-5
DP1, ALPHA1, RHO_REF1 = 1.0e-4, 1.0e-3, 2000.0
MAT_A = dict(cp=1000.0, rho=2000.0, zone='A')
MAT_B = dict(cp=2000.0, rho=1000.0, zone='B')
TOL_T = 1.0e-6

# --- L3, L4: inlets ---------------------------------------------------------------
GP, U_IN, T_IN = 0.01, 10.0, 300.0      # mass flux [kg/(s m^2)], speed, temperature
RHO_IN = GP / U_IN
R_SMALL, R_LARGE = 1.0e-5, 2.0e-5
RHO_E = 1.0e-20                          # the initial content of the domain
RHO_M3 = 2000.0
DX3 = 0.05
T_END3 = 0.05
TOL_N = 1.0e-12

# --- L2: evaporation, case F's material and gas -------------------------------------
TOL_K = 1.0e-6


def as_written(x):
    """A value as the %.15E of the IC file gives it back."""
    return float('%.15E' % x)


def as_output(x):
    """An IC value as the solution file prints it, to 15 significant digits."""
    return float('%.14E' % as_written(x))


def vector(*values):
    """An [ICE-Physics] value per material."""
    return ' '.join('%r' % v for v in values)


# ---------------------------------------------------------------------------
#  L1 - heating, each family on its own material
# ---------------------------------------------------------------------------

def heating(rep, label, phase, zones, fam_mat):
    """Families at rest heating in a gas at rest; family p is material fam_mat[p]."""
    case = Case(WORK / label.lower(), nx=8, Lx=1.0)
    written = case.properties_zones(250, 450, zones)
    rho_p = as_written(ALPHA1 * RHO_REF1)
    n = as_written(rho_p / (RHO_REF1 * 4.0 / 3.0 * math.pi * (0.5 * DP1) ** 3))
    taus = []
    for z in written:
        cp, rho_m = z['cp'][0], z['rho'][0]
        Rp = (0.75 * rho_p / (n * math.pi * rho_m)) ** (1.0 / 3.0)
        taus.append(rho_p * cp / (4.0 * math.pi * K1 * Rp * n))
    t_end, dt_max = taus[0], min(taus) / 1000.0
    case.families(*[(rho_p, 0.0, 0.0, 0.0, TP1, n)] * len(fam_mat))
    case.gas(rho=RHO_G1, u=0.0, v=0.0, w=0.0, T=TG1, R=R1, gam=GAM1, k=K1, mu=MU1)
    case.boundaries('extrapolation')
    case.phase(*phase)
    physics = {'emissivity': vector(*[0.0] * len(zones))} if len(zones) > 1 else None
    case.ini(t_end=t_end, drag='NoDrag', heat='Stokes', dt_max=dt_max, physics=physics,
             closures=('MK',) * len(fam_mat))
    print('   %-3s phase %s, tau = %s s' % (label, ' / '.join(phase), ', '.join('%.6e' % t for t in taus)))

    sol = try_run(rep, case.run)
    if sol is None:
        return
    t = sol['time']
    for p, m in enumerate(fam_mat):
        temp = sol['var'][p * NV + T]
        ref = TG1 + (TP1 - TG1) * math.exp(-t / taus[m])
        err = max(abs(x - ref) for x in temp) / (TG1 - TP1)
        others = [TG1 + (TP1 - TG1) * math.exp(-t / tau) for q, tau in enumerate(taus) if q != m]
        miss = ', '.join('%.3e' % (abs(o - ref) / (TG1 - TP1)) for o in others) or '-'
        print('       family %d (material %s): T_p ICE %.9f K, oracle %.9f K, error %.2e of the gap '
              '(other material: %s)' % (p + 1, zones[m]['zone'], temp[0], ref, err, miss))
        rep.check(err <= TOL_T, '%s: family %d heats as material %s to %.1e of the gap'
                  % (label, p + 1, zones[m]['zone'], err))


# ---------------------------------------------------------------------------
#  L2 - per-material evaporation tokens (case F's material and gas)
# ---------------------------------------------------------------------------

def evaporating(label, phase, zones=None, model=None, boiling=None):
    """Case F's frozen-temperature cloud, one family per material line of phase.

    zones: the table (one per material), or None for the INI constants of one material.
    model: [ICE-Physics] evaporation, or None to leave the key out.
    boiling: (key, values), the boiling temperature per material under that key, or None for case F's.
    """
    case = Case(WORK / label.lower(), nx=8, Lx=1.0)
    nmat = len(phase)
    if zones is not None:
        case.properties_zones(250, 1000, zones)
    case.families(*[(F.RHO_P0, 0.0, 0.0, 0.0, F.TP0, F.NDENS)] * nmat)
    case.gas(rho=F.RHO_G, u=0.0, v=0.0, w=0.0, T=F.TG, R=F.RGAS, gam=F.GAM, k=F.KG, mu=F.MU)
    case.boundaries('extrapolation')
    case.phase(*phase)
    phys = dict((k, vector(*[v] * nmat)) for k, v in F.VAPOUR.items())
    if boiling is not None:
        del phys['boiling-temperature']
        phys[boiling[0]] = vector(*boiling[1])
    if nmat > 1:
        phys['emissivity'] = vector(*[0.0] * nmat)
    if model is not None:
        phys['evaporation'] = model
    case.ini(t_end=0.04, cfl=0.8, rk='RK2', drag='Stokes', heat='Stokes', rho_al=F.RHO_L,
             cs=F.CS_FROZEN, dt_max=2.0e-5, physics=phys, closures=('MK',) * nmat)
    return case


def evap_legs(rep):
    K_ref = F.slope('CEM')

    def slope_check(label, sol, p):
        rho = sol['var'][p * NV + RHO]
        K_ice = (F.D0 ** 2 - F.diameter(rho[0]) ** 2) / sol['time']
        err = abs(K_ice - K_ref) / K_ref
        print('       family %d: K ICE %.8e, CEM oracle %.8e, rel. error %.2e'
              % (p + 1, K_ice, K_ref, err))
        rep.check(err <= TOL_K and max(rho) - min(rho) <= 1.0e-12 * F.RHO_P0,
                  '%s: family %d evaporates as CEM, d-squared slope to %.1e' % (label, p + 1, err))

    def still_check(label, sol, p):
        rho, n = sol['var'][p * NV + RHO], sol['var'][p * NV + N]
        rho0, n0 = as_output(F.RHO_P0), as_output(F.NDENS)
        print('       family %d: rho_p %.15e (IC %.15e), n %.15e' % (p + 1, rho[0], rho0, n[0]))
        rep.check(all(x == rho0 for x in rho) and all(x == n0 for x in n),
                  '%s: family %d does not evaporate (rho_p and n kept to the 15 digits printed)' % (label, p + 1))

    def log_check(label, case, texts):
        log = (case.dir / 'log').read_text()
        rep.check(all(s in log for s in texts), '%s: the setup reports %s' % (label, '; '.join(texts)))

    two = [dict(cp=F.CS_FROZEN, rho=F.RHO_L, zone='A'), dict(cp=2000.0, rho=2000.0, zone='B')]

    print('   L2  "A 1 evaporation=CEM" / "B 1", no [ICE-Physics] evaporation')
    case = evaporating('L2', ('A 1 evaporation=CEM', 'B 1'), two)
    sol = try_run(rep, case.run)
    if sol is not None:
        slope_check('L2', sol, 0)
        still_check('L2', sol, 1)
        log_check('L2', case, ('material 1 (A): CEM', 'material 2 (B): none'))
        first = (case.dir / 'OUTPUT/part-field.tec').read_bytes()
        print('   L2b the same with evaporation = none in the INI')
        case_b = evaporating('L2b', ('A 1 evaporation=CEM', 'B 1'), two, model='none')
        if try_run(rep, case_b.run) is not None:
            rep.check((case_b.dir / 'OUTPUT/part-field.tec').read_bytes() == first,
                      'L2b: evaporation = none in the INI gives the same part-field.tec byte for byte')

    print('   L2d "A 1" / "B 1 evaporation=none", evaporation = CEM in the INI')
    case = evaporating('L2d', ('A 1', 'B 1 evaporation=none'), two, model='CEM')
    sol = try_run(rep, case.run)
    if sol is not None:
        slope_check('L2d', sol, 0)
        still_check('L2d', sol, 1)
        log_check('L2d', case, ('material 1 (A): CEM', 'material 2 (B): none'))

    print('   L2c one material, "A 1 evaporation=CEM", no table, no [ICE-Physics] evaporation')
    sol = try_run(rep, evaporating('L2c', ('A 1 evaporation=CEM',)).run)
    if sol is not None:
        slope_check('L2c', sol, 0)

    print('   L2e control: one material, no token, evaporation = CEM in the INI')
    sol = try_run(rep, evaporating('L2e', ('A 1',), model='CEM').run)
    if sol is not None:
        slope_check('L2e', sol, 0)

    same = [dict(cp=F.CS_FROZEN, rho=F.RHO_L, zone='A'), dict(cp=F.CS_FROZEN, rho=F.RHO_L, zone='B')]
    tb = (F.TBOIL, F.TBOIL + 30.0)
    out = {}
    for label, key, values in (('L2f', 'boiling-temperature', tb), ('L2g', 'Tboil', tb),
                               ('L2h', 'boiling-temperature', tb[::-1])):
        print('   %s "A 1 evaporation=CEM" / "B 1 evaporation=CEM", %s = %s' % (label, key, vector(*values)))
        case = evaporating(label, ('A 1 evaporation=CEM', 'B 1 evaporation=CEM'), same, boiling=(key, values))
        if try_run(rep, case.run) is not None:
            out[label] = (case.dir / 'OUTPUT/part-field.tec').read_bytes()
    if 'L2f' in out and 'L2g' in out:
        rep.check(out['L2g'] == out['L2f'], 'L2g: the alias Tboil gives L2f\'s part-field.tec byte for byte')
    if 'L2f' in out and 'L2h' in out:
        rep.check(out['L2h'] != out['L2f'], 'L2h: the two values swapped change part-field.tec, '
                  'so each reaches its own material')


# ---------------------------------------------------------------------------
#  L3, L4 - inlet records per family, inlet density per material
# ---------------------------------------------------------------------------

def empty_state(nfam):
    """The initial content: 1e-20 kg/m^3 moving at the inlet speed, the same for every family."""
    n = RHO_E / (RHO_M3 * 4.0 / 3.0 * math.pi * R_SMALL ** 3)
    return [(RHO_E, U_IN, 0.0, 0.0, T_IN, n)] * nfam


def inlet_check(rep, label, sol, radii, rho_mats):
    """The filled cells of one mesh block: each family's n rho_mat 4/3 pi r^3 / rho_p = 1,
    and n_1/n_2 = (r_2/r_1)^3 rho_mat2/rho_mat1."""
    rho1 = sol['var'][RHO]
    filled = [c for c in range(len(rho1)) if rho1[c] >= 0.5 * RHO_IN]
    expect = (radii[1] / radii[0]) ** 3 * rho_mats[1] / rho_mats[0]
    ratio_err, own_err = 0.0, 0.0
    for c in filled:
        n1, n2 = sol['var'][N][c], sol['var'][NV + N][c]
        ratio_err = max(ratio_err, abs(n1 / n2 / expect - 1.0))
        for p in (0, 1):
            vol = 4.0 / 3.0 * math.pi * radii[p] ** 3
            own_err = max(own_err, abs(sol['var'][p * NV + N][c] * rho_mats[p] * vol
                                       / sol['var'][p * NV + RHO][c] - 1.0))
    seen = sol['var'][N][filled[0]] / sol['var'][NV + N][filled[0]] if filled else float('nan')
    print('       %2d cells filled of %d; n_1/n_2 = %.15g (expected %.15g); rel. errors: ratio %.2e, '
          'per family %.2e' % (len(filled), len(rho1), seen, expect, ratio_err, own_err))
    rep.check(len(filled) >= 3, '%s: the inlet fills %d cells' % (label, len(filled)))
    if filled:
        rep.check(ratio_err <= TOL_N, '%s: n_1/n_2 = %.6g to %.1e' % (label, expect, ratio_err))
        rep.check(own_err <= TOL_N, '%s: each family injects n = rho_p / (rho_mat 4/3 pi r^3) to %.1e'
                  % (label, own_err))


def inlet_case(rep, label, phase, inlets, radii, rho_mats, zones=None):
    """One mesh block, uncoupled, two families; inlets as Case.boundaries takes them."""
    case = Case(WORK / label.lower(), nx=int(round(1.0 / DX3)), Lx=1.0)
    if zones is not None:
        case.properties_zones(250, 450, zones)
    case.families(*empty_state(2))
    case.boundaries('extrapolation', inlets=inlets)
    case.phase(*phase)
    physics = {'emissivity': vector(0.0, 0.0)} if zones is not None and len(zones) > 1 else None
    case.ini(t_end=T_END3, drag='NoDrag', heat='NoHeat', rho_al=RHO_M3, physics=physics,
             closures=('MK', 'MK'))
    print('   %-3s phase %s, %d record block(s), radii %s um, rho_mat %s' % (
        label, ' / '.join(phase), len(inlets), '/'.join('%g' % (r * 1e6) for r in radii),
        '/'.join('%g' % r for r in rho_mats)))
    sol = try_run(rep, case.run)
    if sol is not None:
        inlet_check(rep, label, sol, radii, rho_mats)


def two_blocks(rep, label, radii, atlas):
    """Two disjoint mesh blocks of 20 and 12 cells, each with its face-1 inlet; radii[b][p].
    atlas: one block of records per family inside each mesh block; otherwise one record per
    face, fanned out over the families (radii[b][0] then serves both)."""
    case = Case(WORK / label.lower(), nx=20, Lx=1.0)
    blocks = []
    rows = []
    names = ['%s%d' % (s, p) for p in (1, 2) for s in ('rp', 'up', 'vp', 'wp', 'Tp', 'np')]
    for b, (nx, y0) in enumerate(((20, 0.0), (12, 1.0)), 1):
        xn = linspace(0.0, nx * DX3, nx + 1)
        cells = [[v] * nx for state in empty_state(2) for v in state]
        blocks.append((xn, [y0, y0 + DX3], [0.0, DX3], cells, 'B%d-CD' % b))
        pays = radii[b - 1] if atlas else radii[b - 1][:1]
        rows += [r for rp in pays for r in bc_rows(nx, 1, 1, block=b, inlet=inlet_payload(GP, U_IN, T_IN, rp))]
    write_tec_blocks(case.dir / 'INPUT/part-ic.tec', names, blocks)
    (case.dir / 'INPUT/part-bc.txt').write_text('\n'.join(rows) + '\n')
    case.phase('A 2')
    case.ini(t_end=0.03, drag='NoDrag', heat='NoHeat', rho_al=RHO_M3, closures=('MK', 'MK'))
    print('   %-3s two mesh blocks, %s, radii %s um' % (
        label, 'ATLAS order' if atlas else 'one record per face',
        ' | '.join('/'.join('%g' % (r * 1e6) for r in rb) for rb in radii)))
    if try_run(rep, lambda: run_ice(case.dir) or True) is None:
        return
    for b in (0, 1):
        sol = read_solution(case.dir / 'OUTPUT/part-field.tec', zone=b)
        print('       mesh block %d:' % (b + 1))
        inlet_check(rep, '%s block %d' % (label, b + 1), sol, radii[b], (RHO_M3, RHO_M3))


def main():
    rep = Report('L. Several materials')

    heating(rep, 'L1', ('A 1', 'B 1'), [MAT_A, MAT_B], [0, 1])
    heating(rep, 'L1b', ('A 2', 'B 1'), [MAT_A, MAT_B], [0, 0, 1])
    heating(rep, 'L1c', ('A 2',), [MAT_A], [0, 0])

    evap_legs(rep)

    pay = lambda r: inlet_payload(GP, U_IN, T_IN, r)                          # noqa: E731
    inlet_case(rep, 'L3', ('A 2',), [pay(R_SMALL), pay(R_LARGE)], (R_SMALL, R_LARGE), (RHO_M3, RHO_M3))
    inlet_case(rep, 'L3c', ('A 2',), [pay(R_SMALL), pay(R_SMALL)], (R_SMALL, R_SMALL), (RHO_M3, RHO_M3))
    two_blocks(rep, 'L3b', ((R_SMALL, R_LARGE), (R_LARGE, R_SMALL)), atlas=True)
    two_blocks(rep, 'L3d', ((R_SMALL, R_SMALL), (R_LARGE, R_LARGE)), atlas=False)

    inlet_case(rep, 'L4', ('A 1', 'B 1'), [pay(R_SMALL)], (R_SMALL, R_SMALL), (2000.0, 1000.0),
               zones=[dict(cp=1000.0, rho=2000.0, zone='A'), dict(cp=1000.0, rho=1000.0, zone='B')])
    inlet_case(rep, 'L4c', ('A 2',), [pay(R_SMALL)], (R_SMALL, R_SMALL), (2000.0, 2000.0),
               zones=[dict(cp=1000.0, rho=2000.0, zone='A')])

    rep.close(WORK)


if __name__ == '__main__':
    main()
