"""J. The time step bounded by the particle relaxation time (tau-factor).

The relaxation sources are integrated explicitly, and the CFL condition knows nothing
about them: a cloud at rest in MK has no signal speed, so before the tau limit the only
bound on its step was dt-max. These particles are small enough that

    tau_p = rho_al dp^2 / (18 mu) = 2.2e-5 s = dt-max / 4.5,

so a step of dt-max puts the SSP-RK3 relaxation at z = -4.5, outside its stability
interval: the amplification R(z) = 1 + z + z^2/2 + z^3/6 is -8.56 and the slip grows
by that factor every step. With tau-factor = 1 (the default) the step is min(CFL step,
dt-max, tau_p), R(-1) = 1/3, and the cloud follows u_p = u_g (1 - exp(-t/tau_p)). The
first step is the exception: tau is not known before the first source evaluation, so
that step is dt-max alone and overshoots to u_p = 95.6; the steps after it are exactly
tau_p, which the end time shows as (t - dt-max)/tau_p being an integer.

J1      one uniform cloud, time accurate: the limit holds and every step after the
        first is tau_p. An accuracy variant sets tau-factor = 0.2 and dt-max = 5e-6 and
        is compared with the exponential at the time reached.
J2      two clouds, one per half of the slab in y (so that no particle crosses from one
        to the other), with dp = 2 and 20 um. Time accurate, the global step follows the
        smaller tau_p, and the large particles, stepped at dt/tau = 0.01, stay on their
        own exponential. With local steps every cell keeps its own limit instead: the
        small particles take tau-sized steps while the large ones keep dt-max.
J3      one half holds a dilute cloud below the density at which a cell counts as empty
        (rho_p = 1e-7 < 1e-6) whose tau_p is 0.61 times the populated one. Its tau must
        not set the step: the steps stay tau_p of the populated half, which a step of
        0.61 tau_p would break, and the dilute half keeps rho_p and n bit for bit.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Physics, Report, try_run            # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)
DT_MAX = 1.0e-4              # the default dt-max: the step before tau is known
NX, NY = 8, 2                # the two-population slabs: one population per row of cells
Y_SPLIT = 0.5 * NY / NX      # between the two rows
RHO_DILUTE = 1.0e-7          # below the empty-cell density of 1e-6
N_DILUTE = 2.5e7             # Rp = 0.78 um, tau_p = 0.61 tau_p of the 1 um cloud


def physics(dp=2.0e-6):
    return Physics(dp=dp, rho_al=2000.0, mu=2e-5, rho_g=1.0, ug=10.0, kg=0.025, cs=1000.0)


def exact(ph, t, tau):
    return ph.ug * (1.0 - math.exp(-t / tau))


def finite(sol):
    return all(math.isfinite(x) for var in sol['var'] for x in var)


def is_integer(k, tol=1e-6):
    return abs(k - round(k)) <= tol


def slab(name, ph, rho, n, ny=1):
    """Cloud at rest in a uniform gas stream at u_g = 10 m/s; ny = 2 splits it in y."""
    if ny == 1:
        case = Case(WORK / name, nx=NX, Lx=1.0)
    else:
        case = Case(WORK / name, nx=NX, ny=ny, Lx=1.0, Ly=ny / NX, Lz=1.0 / NX)
    case.particles(rho=rho, u=0.0, v=0.0, w=0.0, T=ph.Tg, n=n)
    case.gas(rho=ph.rho_g, u=ph.ug, v=0.0, w=0.0, T=ph.Tg, R=ph.R, gam=ph.gam, k=ph.kg, mu=ph.mu)
    case.boundaries('extrapolation')
    return case


def rows(sol, v):
    """Values of variable v on the lower and the upper row of an ny = 2 slab."""
    return sol['var'][v][:NX], sol['var'][v][NX:2 * NX]


def main():
    rep = Report('J. Relaxation-time limit')
    ph = physics()
    tau = ph.tau_stokes
    print('   tau_p = rho_al dp^2 / (18 mu) = %.6e s = dt-max / %.2f' % (tau, DT_MAX / tau))

    # --- J1: the limit, time accurate --------------------------------------------------
    case = slab('j1', ph, ph.rho_p, ph.n)
    case.ini(t_end=20.0 * tau, cfl=0.8, rk='RK3', drag='Stokes', heat='Stokes',
             rho_al=ph.rho_al, cs=ph.cs)
    sol = try_run(rep, case.run)
    if sol is not None:
        u, t = sol['var'][U], sol['time']
        dev = max(abs(x / ph.ug - 1.0) for x in u) if finite(sol) else float('inf')
        k = (t - DT_MAX) / tau
        print('   J1   t = %.6e s, (t - dt-max)/tau = %.8f, max|u_p/u_g - 1| = %.3e' % (t, k, dev))
        rep.check(dev < 1e-3, 'J1: the cloud has relaxed, u_p/u_g = 1 to %.1e' % dev)
        rep.check(k >= 1.0 and is_integer(k),
                  'J1: every step after the first is tau_p, (t - dt-max)/tau = %.6f' % k)

    # --- J1 accuracy variant: tau-factor 0.2, a first step of 0.225 tau ---------------
    case = slab('j1-accuracy', ph, ph.rho_p, ph.n)
    case.ini(t_end=2.0 * tau, cfl=0.8, rk='RK3', drag='Stokes', heat='Stokes',
             rho_al=ph.rho_al, cs=ph.cs, numerics={'tau-factor': 0.2, 'dt-max': 5e-6})
    sol = try_run(rep, case.run)
    if sol is not None:
        u, t = sol['var'][U], sol['time']
        err = max(abs(x - exact(ph, t, tau)) for x in u) / ph.ug if finite(sol) else float('inf')
        k = (t - 5e-6) / (0.2 * tau)
        print('   J1a  t = %.6e s = %.4f tau, (t - 5e-6)/(0.2 tau) = %.8f, max|u_p - exact|/u_g = %.3e'
              % (t, t / tau, k, err))
        rep.check(err < 5e-3, 'J1a: tau-factor 0.2 follows the exponential to %.1e' % err)
        rep.check(k >= 1.0 and is_integer(k), 'J1a: the steps are tau-factor x tau, %.6f of them' % k)

    # --- J2: two radii, time accurate: the global step follows the smaller tau_p ------
    small, large = physics(2.0e-6), physics(2.0e-5)
    tauS, tauL = small.tau_stokes, large.tau_stokes
    two = lambda x, y: small.n if y < Y_SPLIT else large.n
    case = slab('j2', small, small.rho_p, two, ny=NY)
    case.ini(t_end=20.0 * tauS, cfl=0.8, rk='RK3', drag='Stokes', heat='Stokes',
             rho_al=small.rho_al, cs=small.cs)
    sol = try_run(rep, case.run)
    if sol is not None:
        (us, ul), t = rows(sol, U), sol['time']
        ok = finite(sol)
        ds = max(abs(x / small.ug - 1.0) for x in us) if ok else float('inf')
        dl = max(abs(x - exact(large, t, tauL)) for x in ul) / large.ug if ok else float('inf')
        k = (t - DT_MAX) / tauS
        print('   J2   t = %.6e s, (t - dt-max)/tau_S = %.8f, small: max|u_p/u_g - 1| = %.3e, '
              'large: max|u_p - exact|/u_g = %.3e' % (t, k, ds, dl))
        rep.check(ds < 1e-3, 'J2: the small particles have relaxed, to %.1e' % ds)
        rep.check(dl < 1e-4, 'J2: the large particles follow their own exponential, to %.1e' % dl)
        rep.check(k >= 1.0 and is_integer(k),
                  'J2: the global step is the smaller tau_p, (t - dt-max)/tau_S = %.6f' % k)

    # --- J2 local twin: every cell keeps its own limit ---------------------------------
    case = slab('j2-local', small, small.rho_p, two, ny=NY)
    case.ini(iters=20, cfl=0.8, rk='RK3', drag='Stokes', heat='Stokes',
             rho_al=small.rho_al, cs=small.cs, numerics={'time-accurate': False})
    sol = try_run(rep, case.run)
    if sol is not None:
        us, ul = rows(sol, U)
        ok = finite(sol)
        ds = max(abs(x / small.ug - 1.0) for x in us) if ok else float('inf')
        # tau_L > dt-max, so the large particles take 20 steps of dt-max whatever the small ones do
        dl = max(abs(x - exact(large, 20 * DT_MAX, tauL)) for x in ul) / large.ug if ok else float('inf')
        print('   J2l  20 local steps: small: max|u_p/u_g - 1| = %.3e, large: max|u_p - exact(20 dt-max)|/u_g = %.3e'
              % (ds, dl))
        rep.check(ds < 1e-6, 'J2 local: the small particles have relaxed, to %.1e' % ds)
        rep.check(dl < 1e-4, 'J2 local: the large particles kept dt-max, to %.1e' % dl)

    # --- J3: a dilute half below the empty-cell density does not set the step ---------
    rp_dilute = (0.75 * RHO_DILUTE / (N_DILUTE * math.pi * small.rho_al)) ** (1.0 / 3.0)
    print('   J3   dilute cloud: Rp = %.4e m, tau_p = %.4f tau_S' % (rp_dilute, (rp_dilute / small.rp) ** 2))
    rho3 = lambda x, y: small.rho_p if y < Y_SPLIT else RHO_DILUTE
    n3 = lambda x, y: small.n if y < Y_SPLIT else N_DILUTE
    case = slab('j3', small, rho3, n3, ny=NY)
    case.ini(t_end=20.0 * tauS, iters=40, cfl=0.8, rk='RK3', drag='Stokes', heat='Stokes',
             rho_al=small.rho_al, cs=small.cs)
    sol = try_run(rep, case.run)
    if sol is not None:
        (us, ud), t = rows(sol, U), sol['time']
        rho_d, n_d = rows(sol, RHO)[1], rows(sol, N)[1]
        ok = finite(sol)
        ds = max(abs(x / small.ug - 1.0) for x in us) if ok else float('inf')
        k = (t - DT_MAX) / tauS
        print('   J3   t = %.6e s (t_end %.6e), (t - dt-max)/tau_S = %.8f, populated: max|u_p/u_g - 1| = %.3e'
              % (t, 20.0 * tauS, k, ds))
        print('        dilute half: rho_p %s, n %s' % (sorted(set(rho_d)), sorted(set(n_d))))
        rep.check(ok and t >= 20.0 * tauS * (1.0 - 1e-12), 'J3: the run reached t_end within 40 steps')
        rep.check(ds < 1e-3, 'J3: the populated half has relaxed, to %.1e' % ds)
        rep.check(k >= 1.0 and is_integer(k),
                  'J3: the dilute half did not set the step, (t - dt-max)/tau_S = %.6f' % k)
        rep.check(all(r == RHO_DILUTE for r in rho_d) and all(x == N_DILUTE for x in n_d),
                  'J3: the dilute half keeps rho_p and n bit for bit')

    rep.close(WORK)


if __name__ == '__main__':
    main()
