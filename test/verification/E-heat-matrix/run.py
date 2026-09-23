"""E. Every Nusselt correlation, against an independent integration.

The counterpart of case D for convective heat exchange. The cloud starts hot and
slipping, so the velocity and the temperature relax together: Re and Ma change while
the particles cool, and the Nusselt number follows them. The reference is an RK4
integration of the same coupled pair,

    du/dt = (u_g - u)/tau_Stokes,
    dT/dt = 2 Nu(Re, Pr, Ma) kg pi Rp n (T_g - T_p) / (rho_p cs),

with Nu evaluated in Python from Lib_Heat.f90's documented forms.

The laws that tend to Nu = 2 as the slip vanishes are additionally checked against the
exact relaxation of case B, which is independent of ICE.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report                            # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)

LAWS = ('Stokes', 'JAXA1', 'JAXA2', 'JAXA3', 'Chang', 'Ranz-Marshall', 'Kavanau-Drake')

# Laws whose Nusselt number tends to 2 at vanishing slip. JAXA1 tends to 0 and the
# two laws with a Mach correction keep a finite Ma/Re ratio there, so they are out.
NU2_LIMIT = ('Stokes', 'JAXA2', 'Chang', 'Ranz-Marshall')

TP0 = 400.0
T_END = 0.03


def heat_case(ph, law, t_end, name):
    case = Case(WORK / name, nx=8, Lx=1.0)
    case.particles(rho=ph.rho_p, u=0.0, v=0.0, w=0.0, T=TP0, n=ph.n)
    case.gas(rho=ph.rho_g, u=ph.ug, v=0.0, w=0.0, T=ph.Tg,
             R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    case.ini(t_end=t_end, drag='Stokes', heat=law, rho_al=ph.rho_al, cs=ph.cs)
    return case.run()


def main():
    rep = Report('E. Nusselt correlations')

    # --- Vanishing slip: the laws that must give Nu = 2 -----------------------------
    slow = Physics(ug=1.0e-3)
    Re = 2.0 * slow.rho_g * slow.rp * slow.ug / slow.mu
    print('   Nu = 2 limit (Re = %.1e), relative to the exact Nu = 2 relaxation:' % Re)
    for law in NU2_LIMIT:
        sol = heat_case(slow, law, T_END, 'lownu-%s' % law)
        ref = slow.Tg + (TP0 - slow.Tg) * math.exp(-sol['time'] / slow.tau_thermal)
        err = abs(sol['var'][T][0] - ref) / (TP0 - slow.Tg)
        print('     %-16s %.3e' % (law, err))
        rep.check(err <= 0.05, '%s tends to Nu = 2 (%.1e)' % (law, err))

    # --- Finite slip: against an RK4 integration of the coupled system --------------
    ph = Physics()
    Re = 2.0 * ph.rho_g * ph.rp * ph.ug / ph.mu
    print('   Finite slip (Re = %.1f at t = 0), T_p after %.3f s:' % (Re, T_END))
    print('     %-16s %-14s %-14s %s' % ('law', 'ICE [K]', 'RK4 ref [K]', 'difference'))
    for law in LAWS:
        sol = heat_case(ph, law, T_END, 'finite-%s' % law)
        temp, u = sol['var'][T][0], sol['var'][U][0]
        uref, tref = ph.reference(sol['time'], drag='Stokes', heat=law, Tp0=TP0)
        err = abs(temp - tref) / (TP0 - ph.Tg)
        print('     %-16s %-14.7f %-14.7f %.2e' % (law, temp, tref, err))
        rep.check(err <= 1e-3, '%s matches the reference integration (%.1e)' % (law, err))
        rep.check(abs(u - uref) / ph.ug <= 1e-3,
                  '%s leaves the Stokes velocity relaxation intact' % law)
        rep.check(max(temp_i for temp_i in sol['var'][T]) - min(sol['var'][T]) <= 1e-10 * TP0,
                  '%s stays uniform' % law)

    rep.close(WORK)


if __name__ == '__main__':
    main()
