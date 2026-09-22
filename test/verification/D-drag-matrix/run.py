"""D. Every drag correlation, against an independent integration.

Case A verifies one law (Stokes) against a closed form. ICE offers thirteen, and the
rest have no coverage at all. Each one is run in the same uniform box and checked
twice:

  * low Reynolds number - the laws that are built as a Stokes law plus a correction
    must reproduce the exact Stokes exponential once the correction is negligible.
    This is independent of how the correlation is coded in ICE.
  * finite Reynolds number - the relaxation is compared with an RK4 integration of
    du/dt = (u_g - u)/tau(u) using the same correlation, evaluated in Python from
    Lib_Drag.f90's documented forms. This catches a law that is mis-wired, that
    diverges, or whose integration in the solver goes wrong; it cannot catch a
    formula that is wrong in the same way in both places.

The particles start at the gas temperature, so T_p (and with it the temperature ratio
the compressible laws use) stays fixed and only the momentum equation is exercised.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report                            # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)

LAWS = ('Newton', 'Stokes', 'Schlichting', 'Schiller-Naumann', 'Chang', 'Wen-Yu',
        'Putnam', 'Clift-Gauvin', 'Morsi-Alexander', 'Carlson-Hoglund', 'Henderson',
        'Crowe', 'Hermsen')

# Laws that are a Stokes law times a correction which vanishes as Re -> 0
STOKES_LIMIT = ('Stokes', 'Schlichting', 'Schiller-Naumann', 'Chang', 'Wen-Yu',
                'Putnam', 'Clift-Gauvin', 'Morsi-Alexander')

T_END = 0.03


def drag_case(ph, law, t_end, name, Lx=1.0, Ly=None):
    case = Case(WORK / name, nx=8, Lx=Lx, Ly=Ly)
    case.particles(rho=ph.rho_p, u=0.0, v=ph.vg, w=0.0, T=ph.Tg, n=ph.n)
    case.gas(rho=ph.rho_g, u=ph.ug, v=ph.vg, w=0.0, T=ph.Tg,
             R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=t_end, drag=law, heat='Stokes', rho_al=ph.rho_al, cs=ph.cs)
    return case.run()


def main():
    rep = Report('D. Drag correlations')

    # --- Low Reynolds number: the Stokes limit -------------------------------------
    slow = Physics(ug=1.0e-3)
    Re = 2.0 * slow.rho_g * slow.rp * slow.ug / slow.mu
    print('   Stokes limit (Re = %.1e), relative to the exact Stokes solution:' % Re)
    for law in STOKES_LIMIT:
        sol = drag_case(slow, law, T_END, 'lowre-%s' % law)
        ref = slow.ug * (1.0 - math.exp(-sol['time'] / slow.tau_stokes))
        err = abs(sol['var'][U][0] - ref) / slow.ug
        print('     %-18s %.3e' % (law, err))
        rep.check(err <= 0.01, '%s tends to the Stokes law (%.1e)' % (law, err))

    # --- Finite Reynolds number: against an RK4 integration of the same law ---------
    ph = Physics()
    Re = 2.0 * ph.rho_g * ph.rp * ph.ug / ph.mu
    print('   Finite slip (Re = %.1f at t = 0), u_p/u_g after %.3f s:' % (Re, T_END))
    print('     %-18s %-14s %-14s %s' % ('law', 'ICE', 'RK4 reference', 'difference'))
    for law in LAWS:
        sol = drag_case(ph, law, T_END, 'finite-%s' % law)
        u = sol['var'][U][0]
        ref, _ = ph.reference(sol['time'], drag=law, heat='Stokes')
        err = abs(u - ref) / ph.ug
        print('     %-18s %-14.9f %-14.9f %.2e' % (law, u / ph.ug, ref / ph.ug, err))
        rep.check(err <= 1e-3, '%s matches the reference integration (%.1e)' % (law, err))
        rep.check(max(abs(x - ph.Tg) for x in sol['var'][T]) <= 1e-9 * ph.Tg,
                  '%s leaves the temperature untouched' % law)

    # --- Above Re = 1000, where the branched correlations switch formula -----------
    # A transverse velocity shared with the gas keeps dt small enough for the much
    # shorter relaxation time; the slip, and so Re, stays purely axial.
    fast = Physics(ug=200.0, vg=200.0)
    Re = 2.0 * fast.rho_g * fast.rp * fast.ug / fast.mu
    t_fast = 4.0e-3
    print('   High slip (Re = %.0f, Ma = %.2f at t = 0), u_p/u_g after %.4f s:'
          % (Re, min(fast.ug / fast.sound, 1.0), t_fast))
    for law in LAWS:
        sol = drag_case(fast, law, t_fast, 'highre-%s' % law, Lx=0.01, Ly=1.25e-3)
        u = sol['var'][U][0]
        ref, _ = fast.reference(sol['time'], drag=law, heat='Stokes')
        err = abs(u - ref) / fast.ug
        print('     %-18s %-14.9f %-14.9f %.2e' % (law, u / fast.ug, ref / fast.ug, err))
        rep.check(err <= 1e-3, '%s matches the reference integration above Re 1000 (%.1e)'
                  % (law, err))

    rep.close(WORK)


if __name__ == '__main__':
    main()
