"""F. Evaporating droplet cloud - the d-squared law and the model matrix.

A uniform cloud of droplets sits at rest in a uniform, hotter gas at rest. There is
no slip, so Re = 0 and every correlation collapses onto its stagnant-film value, and
there are no spatial gradients, so the state of every cell follows the same pair of
ODEs. ICE transports the bulk density rho_p and the number density n; n has no source
term, so the droplets shrink but are never destroyed, and the diameter comes back out
of the state as

    d = 2 (0.75 rho_p / (n pi rho_l))^(1/3).

The mass equation is drho_p/dt = n mdot, and since every model here gives mdot
proportional to d, the diameter obeys

    d(d^2)/dt = 4 mdot / (pi d rho_l) = -K,   K constant,

which is the d-squared law. The first part of this case freezes the droplet
temperature - an enormous specific heat, the idealisation the law is derived under -
so K really is constant and d^2(t) = d0^2 - K t is exact for every model. K itself is
computed here from the published correlations, independently of ICE.

Freezing the temperature says nothing about the energy equation, so the second part
releases it and integrates the coupled (rho_p, T_p) system with RK4. That is what
exercises the latent sink and the gas-side heat that ASM and TC compute for
themselves in place of the Nusselt one.

The third part is exact again, and does not need the frozen temperature. Put the
d-squared law together with Nu = 2 and the Miller-Harstad-Bellan blowing factor and
the Stefan number collapses:

    b = -(3/2) Pr tau_d mdot / m_p = ln(1 + B_T)   identically,

because every property in tau_d, Pr and m_p cancels against the ones in mdot. So
f2 = b/(e^b - 1) = ln(1+B_T)/B_T, and the blown convective heat 2 pi d k (T_g-T_p) f2
equals the latent sink -mdot L_v exactly, for any gas, any material and any droplet
temperature. That combination is therefore rigorously isothermal, and the droplet
follows the plain d-squared law with a physical specific heat and nothing frozen.

The fourth part holds the droplet above its boiling point, where the saturation
pressure exceeds the gas pressure and the surface sits on the boiling clamp
X_s = 1 - xsCap. B_M is then about 6e11 but finite, and at frozen temperature each
model follows its d-squared slope; with the Langmuir-Knudsen interface the reference is
the integrated ODE, as in the first part.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import Case, Report                                     # noqa: E402

WORK = Path(__file__).resolve().parent / 'work'
RHO, U, V, W, T, N = range(6)

RU = 8314.46          # universal gas constant [J/(kmol K)], as Lib_Evaporation has it
PATM = 101325.0
XS_CAP = 1.0e-12      # Lib_Evaporation's xsCap: the boiling clamp holds X_s this far below 1

# --- The carrier gas: air at 800 K and one atmosphere, at rest -------------------
TG, RGAS, GAM, MU, KG = 800.0, 287.05, 1.4, 1.8e-5, 0.026
RHO_G = PATM / (RGAS * TG)
CPG = GAM * RGAS / (GAM - 1.0)
PR = MU * CPG / KG
MGAS = RU / RGAS

# --- The droplets: water, at rest, 1e-4 m across ---------------------------------
RHO_L, CS_L, LV = 1000.0, 4180.0, 2.26e6
MV, TBOIL, LE, YINF, ALPHAE = 18.015, 373.15, 1.0, 0.0, 1.0
D0, TP0, ALPHA = 1.0e-4, 300.0, 1.0e-3

RHO_P0 = ALPHA * RHO_L
NDENS = RHO_P0 / (RHO_L * (4.0 / 3.0) * math.pi * (0.5 * D0) ** 3)

# A specific heat this large freezes the droplet temperature: the net heating is
# order 1e6 W/m^3, so T_p moves by under 1e-4 K over the whole of the first part.
CS_FROZEN = 1.0e12

MODELS = ('d2-law', 'CEM', 'CEM-B', 'ASM', 'TC')

VAPOUR = {'latent-heat': LV, 'vapour-molar-mass': MV, 'boiling-temperature': TBOIL,
          'lewis-number': LE, 'vapour-mass-fraction': YINF, 'vapour-specific-heat': 0.0,
          'evaporation-coefficient': ALPHAE}


# ---------------------------------------------------------------------------
#  An independent implementation of Lib_Evaporation at Re = 0
# ---------------------------------------------------------------------------

def surface_state(Tp):
    """Clausius-Clapeyron saturation, the boiling clamp, and the Spalding number."""
    psat = PATM * math.exp(-(LV * MV / RU) * (1.0 / Tp - 1.0 / TBOIL))
    p = RHO_G * RGAS * TG
    Xs = min(psat / p, 1.0 - XS_CAP)
    Ys = molar2mass(Xs)
    BM = (Ys - YINF) / (1.0 - Ys) if Ys > YINF else 0.0
    return p, Xs, Ys, BM


def molar2mass(Xs):
    return Xs * MV / (Xs * MV + (1.0 - Xs) * MGAS)


def mass2molar(Ys):
    return Ys * MGAS / (Ys * MGAS + (1.0 - Ys) * MV)


def _F(B):
    """Abramzon-Sirignano F(B) = (1+B)^0.7 ln(1+B)/B."""
    return (1.0 + B) ** 0.7 * math.log(1.0 + B) / B if B > 1.0e-10 else 1.0


def gas_side_rate(model, d, Tp, Ys, BM):
    """mdot [kg/s] (<= 0) and the gas-side heat [W], at Re = 0.

    Returns (mdot, Qdot, override). Qdot is meaningful only when override is set:
    ASM and TC resolve the heat inside their own film, the rest leave it to the
    Nusselt correlation. The latent sink belongs to neither and is added by the
    closure, so it does not appear here.
    """
    Dv = KG / (RHO_G * CPG * LE)
    lnBM = math.log(1.0 + BM)

    if model == 'd2-law':
        BT = CPG * (TG - Tp) / LV
        if BT <= 0.0:
            return 0.0, 0.0, False
        return -2.0 * math.pi * d * KG / CPG * math.log(1.0 + BT), 0.0, False

    if model == 'CEM':
        return -math.pi * d * RHO_G * Dv * 2.0 * lnBM, 0.0, False

    if model == 'CEM-B':
        Tf = Tp + (TG - Tp) / 3.0                  # 1/3 rule film
        rho_f = RHO_G * TG / Tf
        k_f = KG * (Tf / TG) ** 0.7
        Dv_f = k_f / (rho_f * CPG * LE)
        return -math.pi * d * rho_f * Dv_f * 2.0 * lnBM, 0.0, False

    if model == 'ASM':
        # Sh0 = 2 at Re = 0, so Sh* = 2 + (Sh0-2)/F(BM) = 2 and the mass rate is CEM's
        mdot = -math.pi * d * RHO_G * Dv * 2.0 * lnBM
        Nu_star = 2.0
        BT = BM
        for _ in range(5):
            phi = CPG / CPG * 2.0 / Nu_star / LE
            BT = (1.0 + BM) ** phi - 1.0
            if BT <= 0.0:
                BT, Nu_star = BM, 2.0
                break
            Nu_star = 2.0 + (2.0 - 2.0) / _F(BT)
        Qdot = -mdot * CPG * (TG - Tp) / BT if BT > 0.0 else 0.0
        return mdot, Qdot, True

    if model == 'TC':
        Xs = min(mass2molar(Ys), 1.0 - XS_CAP)
        Xinf = mass2molar(YINF)
        Minf = Xinf * MV + (1.0 - Xinf) * MGAS
        rhs0 = MV / Minf * math.log((1.0 - Xinf) / (1.0 - Xs))
        if rhs0 <= 0.0:
            return 0.0, 0.0, False
        Tts = Tp / TG
        Lev = KG / (CPG * Dv * RHO_G)
        mhat = _tc_solve(rhs0, Tts, Lev)
        mdot = -math.pi * d * RHO_G * Dv * 2.0 * mhat
        chi = min(-mdot * CPG / (math.pi * d * KG * 2.0), 700.0)
        if chi > 1.0e-9:
            Qdot = -mdot * CPG * (TG - Tp) / (math.exp(chi) - 1.0)
        else:
            Qdot = math.pi * d * KG * 2.0 * (TG - Tp)
        return mdot, Qdot, True

    raise ValueError('unknown evaporation model ' + model)


def _tc_solve(rhs0, Tts, Lev):
    """Newton with a bisection safeguard on m + (Ts-1) Lev (f(m/Lev) - 1) = rhs0."""
    mlo = rhs0 / max(Tts, 1.0)
    mhi = rhs0 / min(max(Tts, 1.0e-3), 1.0)
    mhat = rhs0
    for _ in range(30):
        x = mhat / Lev
        if x > 1.0e-6:
            E = math.exp(-x)
            f = x / (1.0 - E)
            fp = (1.0 - E - x * E) / (1.0 - E) ** 2
        else:
            f = 1.0 + 0.5 * x + x * x / 12.0
            fp = 0.5 + x / 6.0
        G = mhat + (Tts - 1.0) * Lev * (f - 1.0) - rhs0
        dG = 1.0 + (Tts - 1.0) * fp
        if abs(G) <= 1.0e-13 * max(1.0, rhs0):
            break
        if G > 0.0:
            mhi = mhat
        else:
            mlo = mhat
        mhat -= G / dG
        if mhat <= mlo or mhat >= mhi:
            mhat = 0.5 * (mlo + mhi)
    return mhat


def evaporate(model, d, Tp, interface='VLE'):
    """The whole of Lib_Evaporation's `evaporation` at Re = 0."""
    if model == 'none':
        return 0.0, 0.0, False
    p, Xs, Ys, BM = surface_state(Tp)
    if Ys <= YINF or BM <= 0.0:
        return 0.0, 0.0, False
    mdot, Qdot, ovr = gas_side_rate(model, d, Tp, Ys, BM)
    if interface == 'LK':
        mdot, Qdot, ovr = _lk(model, d, Tp, p, Xs, mdot, Qdot, ovr)
    return mdot, Qdot, ovr


def _lk(model, d, Tp, p, XsEq, md, Qd, ovr):
    """Langmuir-Knudsen non-equilibrium interface: Xs depressed, Picard on mdot."""
    Sc = PR * LE
    LKd2 = 2.0 * MU * math.sqrt(2.0 * math.pi * Tp * RU / MV) / (ALPHAE * Sc * p) / d
    diff_prev = float('inf')
    mdNew, QdNew, ovrNew = md, Qd, ovr
    for _ in range(30):
        beta = -md * PR / (2.0 * math.pi * MU * d)
        XsN = max(XsEq - LKd2 * beta, 0.0)
        YsN = molar2mass(XsN)
        mdNew, QdNew, ovrNew = 0.0, 0.0, False
        if YsN > YINF:
            BMN = (YsN - YINF) / (1.0 - YsN)
            if BMN > 0.0:
                mdNew, QdNew, ovrNew = gas_side_rate(model, d, Tp, YsN, BMN)
        diff = abs(mdNew - md)
        if diff <= 1.0e-12 * abs(mdNew):
            return mdNew, QdNew, ovrNew
        md = mdNew if diff < diff_prev else 0.5 * (md + mdNew)
        diff_prev = diff
    return mdNew, QdNew, ovrNew


def blowing_factor(d, mp, mdot):
    """Miller-Harstad-Bellan f2 = b/(exp(b)-1), the Stefan-blowing heat reduction."""
    if mdot >= 0.0 or mp <= 0.0:
        return 1.0
    taud = RHO_L * d * d / (18.0 * MU)
    b = -1.5 * PR * taud * mdot / mp
    if b <= 1.0e-10:
        return 1.0
    if b > 500.0:
        return 0.0
    return b / (math.exp(b) - 1.0)


# ---------------------------------------------------------------------------
#  Reference solutions
# ---------------------------------------------------------------------------

def slope(model, interface='VLE', Tp=TP0):
    """K in d^2(t) = d0^2 - K t, at frozen droplet temperature.

    mdot is proportional to d for every model, so K does not depend on d; it is
    evaluated at d0 and divided back out.
    """
    mdot, _, _ = evaporate(model, D0, Tp, interface)
    return -4.0 * mdot / (math.pi * D0 * RHO_L)


def diameter(rho_p):
    return 2.0 * (0.75 * rho_p / (NDENS * math.pi * RHO_L)) ** (1.0 / 3.0)


def coupled(model, t_end, cs, blowing='none', interface='VLE', steps=40000, tp0=TP0):
    """RK4 on (rho_p, T_p) with the temperature released.

    rho_p' = n mdot and, after the enthalpy the leaving mass carries cancels,
    rho_p cs T_p' = n (Qconv + mdot Lv): the convective heat in, the latent heat out.
    """
    h = t_end / steps

    def rhs(rho_p, Tp):
        d = diameter(rho_p)
        mdot, Qevap, ovr = evaporate(model, d, Tp, interface)
        if ovr:
            Qconv = Qevap
        else:
            Qconv = 2.0 * 2.0 * KG * math.pi * (0.5 * d) * (TG - Tp)   # Nu = 2
            if blowing == 'LK':
                mp = RHO_L * math.pi / 6.0 * d ** 3
                Qconv *= blowing_factor(d, mp, mdot)
        return NDENS * mdot, NDENS * (Qconv + mdot * LV) / (rho_p * cs)

    rho_p, Tp = RHO_P0, tp0
    for _ in range(steps):
        k1 = rhs(rho_p, Tp)
        k2 = rhs(rho_p + 0.5 * h * k1[0], Tp + 0.5 * h * k1[1])
        k3 = rhs(rho_p + 0.5 * h * k2[0], Tp + 0.5 * h * k2[1])
        k4 = rhs(rho_p + h * k3[0], Tp + h * k3[1])
        rho_p += h / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        Tp += h / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    return rho_p, Tp


# ---------------------------------------------------------------------------
#  Running ICE
# ---------------------------------------------------------------------------

def run(name, model, t_end, dt_max, cs, blowing='none', interface='VLE', tp0=TP0):
    case = Case(WORK / name, nx=8, Lx=1.0)
    case.particles(rho=RHO_P0, u=0.0, v=0.0, w=0.0, T=tp0, n=NDENS)
    case.gas(rho=RHO_G, u=0.0, v=0.0, w=0.0, T=TG,
             R=RGAS, gam=GAM, k=KG, mu=MU)
    case.boundaries('extrapolation')
    phys = dict(VAPOUR, evaporation=model)
    phys['evaporation-interface'] = interface
    phys['evaporation-blowing'] = blowing
    case.ini(t_end=t_end, cfl=0.8, rk='RK2', drag='Stokes', heat='Stokes',
             rho_al=RHO_L, cs=cs, dt_max=dt_max, physics=phys)
    return case.run()


def try_run(rep, *args, **kw):
    """run(), with a failed run recorded as a FAIL so the legs after it still run."""
    try:
        return run(*args, **kw)
    except RuntimeError as exc:
        rep.check(False, str(exc))
        return None


def main():
    rep = Report('F. Droplet evaporation')

    # --- Part 1: the d-squared law, one run per model, temperature frozen ---------
    t_end, dt_max = 0.04, 2.0e-5
    print('   frozen T_p = %.1f K, T_g = %.1f K, d0 = %.0f um, t_end = %g s'
          % (TP0, TG, D0 * 1e6, t_end))
    print('   model      K (ICE)        K (reference)   rel. error   d^2 left')
    for model in MODELS:
        sol = run('d2-%s' % model, model, t_end, dt_max, CS_FROZEN)
        t = sol['time']
        d2 = diameter(sol['var'][RHO][0]) ** 2
        K_ice = (D0 ** 2 - d2) / t
        K_ref = slope(model)
        err = abs(K_ice - K_ref) / K_ref
        print('   %-9s  %.8e  %.8e  %9.2e  %8.4f'
              % (model, K_ice, K_ref, err, d2 / D0 ** 2))
        rep.check(err <= 1.0e-6, '%s: d-squared slope to %.1e' % (model, err))

        # the field must stay uniform and the droplet count must not move
        rho = sol['var'][RHO]
        rep.check(max(rho) - min(rho) <= 1.0e-12 * RHO_P0, '%s: uniform' % model)
        rep.check(max(abs(x - NDENS) for x in sol['var'][N]) <= 1.0e-9 * NDENS,
                  '%s: number density unchanged' % model)
        rep.check(max(abs(x - TP0) for x in sol['var'][T]) <= 1.0e-3,
                  '%s: temperature held by the large specific heat' % model)

    # --- evaporation = none must leave the mass alone ------------------------------
    sol = run('none', 'none', t_end, dt_max, CS_FROZEN)
    rep.check(abs(sol['var'][RHO][0] - RHO_P0) <= 1.0e-14 * RHO_P0,
              'evaporation = none: bulk density untouched')

    # --- At Re = 0 the ASM film correction is inactive, so its mass rate is CEM's ---
    rep.check(abs(slope('ASM') - slope('CEM')) <= 1.0e-12 * slope('CEM'),
              'ASM and CEM agree at Re = 0 (Sh* = Sh0 = 2)')

    # --- The Langmuir-Knudsen interface depresses the surface fraction -------------
    # The Knudsen-layer thickness goes as 1/d while the rate that drives it does not,
    # so the depression deepens as the droplet shrinks and d^2 stops being linear in
    # time. The reference is therefore the integrated ODE, not a slope.
    sol = run('lk', 'CEM', t_end, dt_max, CS_FROZEN, interface='LK')
    rho_ref, _ = coupled('CEM', sol['time'], CS_FROZEN, interface='LK')
    err = abs(sol['var'][RHO][0] - rho_ref) / RHO_P0
    K_lk, K_vle = slope('CEM', interface='LK'), slope('CEM')
    print('   LK interface: initial K = %.8e vs VLE %.8e (%.3f %% slower)'
          % (K_lk, K_vle, 100.0 * (1.0 - K_lk / K_vle)))
    rep.check(err <= 1.0e-8, 'CEM + LK: bulk density to %.1e' % err)
    rep.check(K_lk < K_vle, 'LK evaporates more slowly than VLE')
    rep.check(abs(sol['var'][RHO][0] - run('lk-vle', 'CEM', t_end, dt_max,
                                           CS_FROZEN)['var'][RHO][0]) > 0.0,
              'the interface treatment reaches the solver')

    # --- Part 2: the temperature released -----------------------------------------
    # tau_thermal = rho_l cs d^2 / (12 k) is 1.3e-4 s here, so the run is short and
    # the step small enough to resolve it. This is the part that tests the latent
    # sink and the two heat paths that part 1 cannot reach.
    t_end, dt_max = 0.05, 5.0e-6
    tau_th = RHO_L * CS_L * D0 ** 2 / (12.0 * KG)
    print('   released T_p: tau_thermal = %.3e s, t_end = %g s' % (tau_th, t_end))
    print('   case       T_p (ICE)    T_p (RK4)     rho_p err')
    free = {}
    for model in ('d2-law', 'ASM', 'TC'):
        sol = run('hot-%s' % model, model, t_end, dt_max, CS_L)
        rho_ice, T_ice = sol['var'][RHO][0], sol['var'][T][0]
        rho_ref, T_ref = coupled(model, sol['time'], CS_L)
        eT = abs(T_ice - T_ref) / (TG - TP0)
        eR = abs(rho_ice - rho_ref) / RHO_P0
        free[model] = T_ice
        print('   %-9s  %10.6f  %10.6f   %9.2e' % (model, T_ice, T_ref, eR))
        rep.check(eT <= 1.0e-8, '%s: T_p to %.1e of the initial gap' % (model, eT))
        rep.check(eR <= 1.0e-9, '%s: bulk density to %.1e' % (model, eR))

    # --- Part 3: d-squared law + Nu = 2 + blowing is exactly isothermal -----------
    # Nothing is frozen here: the specific heat is the real one, and the droplet
    # holds its temperature because the identity above makes the two heat terms
    # cancel term by term. Without the blowing factor the same run heats by ~20 K.
    sol = run('isothermal', 'd2-law', t_end, dt_max, CS_L, blowing='LK')
    drift = max(abs(x - TP0) for x in sol['var'][T])
    K_ice = (D0 ** 2 - diameter(sol['var'][RHO][0]) ** 2) / sol['time']
    K_ref = slope('d2-law')
    err = abs(K_ice - K_ref) / K_ref
    print('   d2-law + blowing: T_p drift %.2e K (without blowing: %+.3f K),'
          ' d-squared slope to %.1e' % (drift, free['d2-law'] - TP0, err))
    rep.check(drift <= 1.0e-8, 'd2-law + blowing is isothermal to %.1e K' % drift)
    rep.check(err <= 1.0e-6, 'd2-law + blowing: d-squared slope to %.1e' % err)
    rep.check(free['d2-law'] - TP0 > 1.0, 'without blowing the same droplet heats up')

    # --- Part 4: above the boiling point ------------------------------------------
    # psat(380 K) = 1.27 atm > p, so every model evaluates at X_s = 1 - xsCap. The run
    # lasts until CEM has lost 10 % of its mass. One ulp of Y_s at the cap moves the
    # rate by 2.5e-6, hence 1e-5 on K; d2-law never reads B_M and TC caps on its own.
    tp_boil = 380.0
    p, Xs, _, BM = surface_state(tp_boil)
    t_end = (1.0 - 0.9 ** (2.0 / 3.0)) * D0 ** 2 / slope('CEM', Tp=tp_boil)
    dt_max = t_end / 1000.0
    print('   boiling: T_p = %.0f K, psat/p = %.3f, B_M = %.3e, t_end = %.4e s'
          % (tp_boil, PATM * math.exp(-(LV * MV / RU) * (1.0 / tp_boil - 1.0 / TBOIL)) / p,
             BM, t_end))
    print('   model      K (ICE)        K (reference)   rel. error')
    for model in MODELS:
        sol = try_run(rep, 'boil-%s' % model, model, t_end, dt_max, CS_FROZEN, tp0=tp_boil)
        if sol is None:
            continue
        K_ice = (D0 ** 2 - diameter(sol['var'][RHO][0]) ** 2) / sol['time']
        K_ref = slope(model, Tp=tp_boil)
        err = abs(K_ice - K_ref) / K_ref
        print('   %-9s  %.8e  %.8e  %9.2e' % (model, K_ice, K_ref, err))
        rep.check(err <= 1.0e-5, '%s boiling: d-squared slope to %.1e' % (model, err))

    # The Langmuir-Knudsen depression takes X_s off the cap, to 0.98 here, so this leg
    # is not ulp-sensitive and keeps the 1e-8 of the leg at 300 K.
    sol = try_run(rep, 'boil-lk', 'CEM', t_end, dt_max, CS_FROZEN, interface='LK', tp0=tp_boil)
    if sol is not None:
        rho_ref, _ = coupled('CEM', sol['time'], CS_FROZEN, interface='LK', steps=4000, tp0=tp_boil)
        err = abs(sol['var'][RHO][0] - rho_ref) / RHO_P0
        print('   CEM + LK boiling: rho_p %.10e vs %.10e (%.2e)' % (sol['var'][RHO][0], rho_ref, err))
        rep.check(err <= 1.0e-8, 'CEM + LK boiling: bulk density to %.1e' % err)

    rep.close(WORK)


if __name__ == '__main__':
    main()
