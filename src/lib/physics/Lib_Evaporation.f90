module ICE_Lib_Evaporation
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    private
    public :: assign_evaporation, assign_interface, assign_blowing, assign_liquid, assign_boiling
    public :: blowingFactor
    public :: evaporation
    public :: psat_CC, invTboil_at
    public :: log1pS, expm1S, chiETC
    public :: boilingTemperatureCC, alphaBoil, zgrBoiling, heatFedRate

    !> Vapour properties, packed contiguously by Assign_Setup and indexed by the
    !  parameters below.
    integer(I4), parameter, public :: nep = 8
    integer(I4), parameter, public :: iMv=1, iLv=2, icpv=3, iLe=4
    integer(I4), parameter, public :: iYinf=5, iLvMvOverRu=6, iinvTboil=7
    !> alphaE is consumed by the Langmuir-Knudsen interface.
    integer(I4), parameter, public :: ialphaE=8

    real(R8), parameter :: pi   = acos(-1._R8)
    real(R8), parameter, public :: Ru = 8314.46_R8  ! universal gas constant [J/(kmol K)]
    real(R8), parameter, public :: Patm = 101325._R8  ! atmospheric pressure [Pa]
    !> The boiling clamp holds Xs this far below 1, so BM = (Ys-Yinf)/(1-Ys) stays finite.
    real(R8), parameter, public :: xsCap = 1.e-12_R8

contains

    subroutine assign_evaporation(evaporation_word,evapSelect)
        implicit none
        character(len=*), intent(in)  :: evaporation_word
        integer(kind=I4), intent(out) :: evapSelect

        select case (evaporation_word)
        case('none')
            evapSelect = 0
        case('d2-law')
            evapSelect = 1
        case('CEM')
            evapSelect = 2
        case('CEM-B')
            evapSelect = 3
        case('ASM')
            evapSelect = 4
        case('TC')
            evapSelect = 5
        case default
            write(*,*)
            write(*,*)
            write(*,*) "Wrong evaporation model input ---> "//trim(evaporation_word)
            write(*,*) "Choose one of the following :"
            write(*,*) "- none  (no mass transfer)"
            write(*,*) "- d2-law"
            write(*,*) "- CEM   (Classical Evaporation Model)"
            write(*,*) "- CEM-B (Corrected CEM) "
            write(*,*) "- ASM   (Abramzon-Sirignano Model)"
            write(*,*) "- TC    (Tonini-Cossali analytical model)"
            write(*,*)
            error stop '[ERROR ICE] unknown evaporation model'
        end select

    end subroutine assign_evaporation

    !> Composable phase-change axes: word -> selector mappers.

    !> Interface: 0 VLE equilibrium, 1 Langmuir-Knudsen non-equilibrium.
    subroutine assign_interface(word, intfSelect)
        implicit none
        character(len=*), intent(in)  :: word
        integer(kind=I4), intent(out) :: intfSelect
        select case (word)
        case('VLE'); intfSelect = 0
        case('LK');  intfSelect = 1
        case default
            write(*,*) "Wrong interface model input ---> "//trim(word)
            write(*,*) "Choose one of: VLE, LK"
            error stop '[ERROR ICE] unknown interface model'
        end select
    end subroutine assign_interface

    !> Liquid-side conduction: 0 ITC (infinite conductivity), 1 P2T (parabolic temperature profile).
    subroutine assign_liquid(word, liqSelect)
        implicit none
        character(len=*), intent(in)  :: word
        integer(kind=I4), intent(out) :: liqSelect
        select case (word)
        case('ITC'); liqSelect = 0
        case('P2T'); liqSelect = 1
        case default
            write(*,*) "Wrong liquid-conduction model input ---> "//trim(word)
            write(*,*) "Choose one of: ITC, P2T"
            error stop '[ERROR ICE] unknown liquid-conduction model'
        end select
    end subroutine assign_liquid

    !> Boiling: 0 Xs=1 clamp, 1 Zuo-Gomes-Rutland branch above the boiling temperature at the local pressure.
    subroutine assign_boiling(word, boilSelect)
        implicit none
        character(len=*), intent(in)  :: word
        integer(kind=I4), intent(out) :: boilSelect
        select case (word)
        case('clamp'); boilSelect = 0
        case('ZGR');   boilSelect = 1
        case default
            write(*,*) "Wrong boiling model input ---> "//trim(word)
            write(*,*) "Choose one of: clamp, ZGR"
            error stop '[ERROR ICE] unknown boiling model'
        end select
    end subroutine assign_boiling

    !> Evaporative heat-transfer reduction f2: 0 none (f2=1), 1 MHB98 eq.19 f2 = b/(exp(b)-1).
    subroutine assign_blowing(word, blowSelect)
        implicit none
        character(len=*), intent(in)  :: word
        integer(kind=I4), intent(out) :: blowSelect
        select case (word)
        case('none'); blowSelect = 0
        case('LK');   blowSelect = 1
        case default
            write(*,*) "Wrong blowing model input ---> "//trim(word)
            write(*,*) "Choose one of: none, LK"
            error stop '[ERROR ICE] unknown blowing model'
        end select
    end subroutine assign_blowing

    !> MHB98 eq.19: f2 = b/(exp(b)-1) with b = -(3/2)*Pr_G*tau_d*mdot/m_d (their eq.17), the
    !  quasi-steady reduction of convective heat transfer by the outgoing vapour (Stefan blowing).
    !  b->0 gives f2->1 exactly; returns 1 whenever the state cannot produce a valid b.
    pure function blowingFactor(gamma, Rg, mug, kg, rhop, dp, mp, mdot) result(f2)
        implicit none
        real(R8), intent(in) :: gamma, Rg, mug, kg, rhop, dp, mp, mdot
        real(R8) :: f2, cpg, Pr, taud, b
        f2 = 1._R8
        if (mdot >= 0._R8 .or. mp <= 0._R8 .or. mug <= 0._R8 .or. kg <= 0._R8) return
        if (dp <= 0._R8 .or. rhop <= 0._R8 .or. gamma <= 1._R8) return
        cpg  = gamma*Rg/(gamma-1._R8)
        Pr   = mug*cpg/kg
        taud = rhop*dp*dp/(18._R8*mug)
        b    = -1.5_R8*Pr*taud*mdot/mp
        if (b /= b .or. b <= 1.e-10_R8) return          ! NaN or vanishing blowing -> f2=1
        if (b > 500._R8) then; f2 = 0._R8; return; endif ! exp overflow guard
        f2 = b/(exp(b)-1._R8)
    end function blowingFactor

    pure subroutine evaporation(rhog, Tg, gamma, Rg, mug, kg, &
                                 Tp, dp, Re, evapSelect, intfSelect, ep, &
                                 mdot, Qdot_evap, override_Qdot, psatExt, BMout)
        implicit none
        real(R8), intent(in)  :: rhog, Tg, gamma, Rg, mug, kg
        real(R8), intent(in)  :: Tp, dp, Re
        integer(kind=I4), intent(in) :: evapSelect, intfSelect
        real(R8), intent(in)  :: ep(nep)
        real(R8), intent(out) :: mdot, Qdot_evap
        logical,  intent(out) :: override_Qdot
        real(R8), intent(in), optional :: psatExt   ! saturation pressure from the property table [Pa]
        real(R8), intent(out), optional :: BMout    ! the Spalding number of the surface state (before LK), 0 if none
        real(R8) :: cpg, p, Mg, Pr, Sc, Re05, psat, Xs, Ys, BM

        mdot = 0._R8
        Qdot_evap = 0._R8
        override_Qdot = .false.
        if (present(BMout)) BMout = 0._R8
        if (evapSelect == 0) return

        !> Common gas properties (computed once)
        cpg  = gamma * Rg / (gamma - 1._R8)
        if (cpg <= 0._R8) error stop '[ERROR ICE] evaporation: cp(gas) <= 0'
        p    = rhog * Rg * Tg
        Mg   = Ru / Rg

        !> Saturation pressure: the property table's Psat when given, else Clausius-Clapeyron
        if (present(psatExt)) then
            psat = psatExt
        else
            psat = psat_CC(Tp, ep(iLvMvOverRu), ep(iinvTboil))
        endif

        !> Surface vapor fraction
        if (psat >= p) then
            Xs = 1._R8 - xsCap  ! boiling regime: clamp
        else
            Xs = min(psat / p, 1._R8 - xsCap)
        endif
        Ys = molar2mass(Xs, ep(iMv), Mg)

        !> Spalding mass transfer number
        if (Ys <= ep(iYinf)) return  ! no evaporation
        BM = (Ys - ep(iYinf)) / (1._R8 - Ys)
        if (BM <= 0._R8) return
        if (present(BMout)) BMout = BM

        !> Re^0.5 (shared by Sh and Nu correlations)
        Re05 = sqrt(max(Re, 0._R8))
        Pr   = mug * cpg / kg
        Sc   = Pr * ep(iLe)

        call gasSideRate(Ys, BM, mdot, Qdot_evap, override_Qdot)

        !> Non-equilibrium interface: depress Xs per Langmuir-Knudsen, Picard on mdot
        if (intfSelect == 1) call lkCorrection(Xs, mdot, Qdot_evap, override_Qdot)

    contains

        !> Gas-side rate at a given surface state; common props via host association.
        pure subroutine gasSideRate(Ysrf, BMsrf, md, Qd, ovr)
            real(R8), intent(in)  :: Ysrf, BMsrf
            real(R8), intent(out) :: md, Qd
            logical,  intent(out) :: ovr
            md = 0._R8; Qd = 0._R8; ovr = .false.
            select case (evapSelect)
            case (1)
                call d2law(Tg, kg, Tp, dp, cpg, ep(iLv), md)
            case (2)
                call CEM_model(rhog, kg, dp, cpg, Sc, Re05, BMsrf, ep(iLe), md)
            case (3)
                call CEMB_model(rhog, Tg, mug, kg, Tp, dp, cpg, Re05, BMsrf, ep(iLe), md)
            case (4)
                call ASM_model(rhog, Tg, kg, Tp, dp, cpg, Pr, Sc, Re05, BMsrf, &
                               ep(iLe), ep(icpv), md, Qd, ovr)
            case (5)
                call TC_model(rhog, Tg, kg, Tp, dp, Mg, cpg, Pr, Sc, Re05, Ysrf, &
                              ep(iLe), ep(icpv), ep(iYinf), ep(iMv), md, Qd, ovr)
            end select
        end subroutine gasSideRate

        !> Langmuir-Knudsen non-equilibrium interface (Miller-Harstad-Bellan 1998, model M2):
        !  Xs_neq = Xs_eq - (2 L_K/d)·beta, beta = -mdot·Pr/(2π·mu_g·d) — implicit, bounded Picard.
        pure subroutine lkCorrection(XsEq, md, Qd, ovr)
            real(R8), intent(in)    :: XsEq
            real(R8), intent(inout) :: md, Qd
            logical,  intent(inout) :: ovr
            integer,  parameter :: maxPicard = 30
            real(R8), parameter :: tolPicard = 1.e-12_R8
            real(R8) :: LKd2, beta, XsN, YsN, BMN, mdNew, QdNew, diff, diffPrev
            logical  :: ovrNew
            integer  :: it

            LKd2 = 2._R8 * mug * sqrt(2._R8*pi*Tp*Ru/ep(iMv)) / (ep(ialphaE) * Sc * p) / dp
            diffPrev = huge(1._R8)
            do it = 1, maxPicard
                beta = -md * Pr / (2._R8 * pi * mug * dp)
                XsN  = max(XsEq - LKd2*beta, 0._R8)
                YsN  = molar2mass(XsN, ep(iMv), Mg)
                mdNew = 0._R8; QdNew = 0._R8; ovrNew = .false.
                if (YsN > ep(iYinf)) then
                    BMN = (YsN - ep(iYinf)) / (1._R8 - YsN)
                    if (BMN > 0._R8) call gasSideRate(YsN, BMN, mdNew, QdNew, ovrNew)
                endif
                diff = abs(mdNew - md)
                if (diff <= tolPicard*abs(mdNew)) then
                    md = mdNew; Qd = QdNew; ovr = ovrNew
                    return
                endif
                !> Plain while contracting; damp on expansion (risk R4 oscillation guard)
                if (diff < diffPrev) then; md = mdNew; else; md = 0.5_R8*(md + mdNew); endif
                diffPrev = diff
            enddo
            md = mdNew; Qd = QdNew; ovr = ovrNew  ! unconverged: last full evaluation (bounded, no NaN)
        end subroutine lkCorrection

    end subroutine evaporation


    !> Surface mole fraction -> mass fraction.
    pure function molar2mass(Xs, Mv, Mg) result(Ys)
        implicit none
        real(R8), intent(in) :: Xs, Mv, Mg
        real(R8) :: Ys
        Ys = Xs * Mv / (Xs * Mv + (1._R8 - Xs) * Mg)
    end function molar2mass


    !> Mass fraction -> mole fraction (exact inverse of molar2mass).
    pure function mass2molar(Ys, Mv, Mg) result(Xs)
        implicit none
        real(R8), intent(in) :: Ys, Mv, Mg
        real(R8) :: Xs
        Xs = Ys * Mg / (Ys * Mg + (1._R8 - Ys) * Mv)
    end function mass2molar


    !===================================================================!
    !  Model 1: d²-law (Godsave 1953)                                  !
    !  Stagnant film (Sh=2, Nu=2), quasi-steady                        !
    !  Ref: Godsave (1953) eq.(13), Spalding (1954)                    !
    !===================================================================!
    pure subroutine d2law(Tg, kg, Tp, dp, cpg, Lv, mdot)
        implicit none
        real(R8), intent(in)  :: Tg, kg, Tp, dp, cpg, Lv
        real(R8), intent(out) :: mdot
        real(R8) :: BT

        BT = cpg * (Tg - Tp) / Lv
        if (BT <= 0._R8) then
            mdot = 0._R8
            return
        endif
        mdot = -2._R8 * pi * dp * kg / cpg * log(1._R8 + BT)  ! Nu=2 stagnant film (Godsave/Spalding)

    end subroutine d2law


    !===================================================================!
    !  Model 2: CEM — Classical Evaporation Model (Spalding 1954)      !
    !  Convective correction via Ranz-Marshall Sh                      !
    !  Ref: Spalding (1954), Fluent diffusion-controlled model         !
    !===================================================================!
    pure subroutine CEM_model(rhog, kg, dp, cpg, Sc, Re05, BM, Le, mdot)
        implicit none
        real(R8), intent(in)  :: rhog, kg, dp, cpg, Sc, Re05, BM, Le
        real(R8), intent(out) :: mdot
        real(R8) :: Dv, Sh

        Dv = kg / (rhog * cpg * Le)
        Sh = 2._R8 + 0.6_R8 * Re05 * Sc**(1._R8/3._R8)
        mdot = -pi * dp * rhog * Dv * Sh * log(1._R8 + BM)

    end subroutine CEM_model


    !===================================================================!
    !  Model 3: CEM-B — CEM with 1/3 rule film correction             !
    !  Ref: Hubbard-Denny-Mills (1975), Abramzon-Sirignano eq.(4)      !
    !===================================================================!
    pure subroutine CEMB_model(rhog, Tg, mug, kg, Tp, dp, cpg, Re05, BM, Le, mdot)
        implicit none
        real(R8), intent(in)  :: rhog, Tg, mug, kg, Tp, dp, cpg, Re05, BM, Le
        real(R8), intent(out) :: mdot
        real(R8) :: Tf, rho_f, mu_f, kl_f, cpg_f, Dv_f, Sc_f, Re_f, Sh_f

        !> Film temperature (1/3 rule)
        Tf = Tp + (Tg - Tp) / 3._R8

        !> Film properties (power-law approximation)
        rho_f = rhog * Tg / Tf             ! ideal gas, p=const
        mu_f  = mug  * (Tf / Tg)**0.7_R8   ! Sutherland approx
        kl_f  = kg   * (Tf / Tg)**0.7_R8
        cpg_f = cpg                         ! weak T-dependence
        Dv_f  = kl_f / (rho_f * cpg_f * Le)
        Sc_f  = mu_f * cpg_f * Le / kl_f
        Re_f  = Re05 * Re05 * rho_f / rhog * mug / mu_f  ! Re_f = Re * rho_f/rho * mu/mu_f

        Sh_f = 2._R8 + 0.6_R8 * sqrt(max(Re_f, 0._R8)) * Sc_f**(1._R8/3._R8)
        mdot = -pi * dp * rho_f * Dv_f * Sh_f * log(1._R8 + BM)

    end subroutine CEMB_model


    !===================================================================!
    !  Model 4: ASM — Abramzon-Sirignano Model (1989)                  !
    !  Extended film with Stefan flow correction                       !
    !  Ref: Abramzon & Sirignano (1989) eq.(8)-(24)                    !
    !===================================================================!
    pure subroutine ASM_model(rhog, Tg, kg, Tp, dp, cpg, Pr, Sc, Re05, BM, &
                              Le, cpv, mdot, Qdot_evap, override_Qdot)
        implicit none
        real(R8), intent(in)  :: rhog, Tg, kg, Tp, dp, cpg, Pr, Sc, Re05, BM
        real(R8), intent(in)  :: Le, cpv
        real(R8), intent(out) :: mdot, Qdot_evap
        logical,  intent(out) :: override_Qdot
        real(R8) :: Dv, lnBM, FM, FT, Sh0, Sh_star, Nu0, Nu_star
        real(R8) :: BT, phi, cpv_eff
        integer  :: iter

        override_Qdot = .true.
        Dv   = kg / (rhog * cpg * Le)
        lnBM = log(1._R8 + BM)

        !> Step 1: Modified Sherwood (mass transfer)
        FM  = F_correction(BM)
        Sh0 = 2._R8 + 0.552_R8 * Re05 * Sc**(1._R8/3._R8)  ! Frossling
        Sh_star = 2._R8 + (Sh0 - 2._R8) / FM

        !> Mass transfer rate [eq. 8]
        mdot = -pi * dp * rhog * Dv * Sh_star * lnBM

        !> Step 2: Modified Nusselt (heat transfer)
        !  cpv_eff: vapor specific heat (from INI or approximation)
        if (cpv > 0._R8) then
            cpv_eff = cpv
        else
            cpv_eff = cpg  ! approximation
        endif

        Nu0 = 2._R8 + 0.552_R8 * Re05 * Pr**(1._R8/3._R8)  ! Frossling

        !> Iterative phi-BT-FT-Nu* loop [eq. 22-23]
        !  Initialize: phi=1 (Le=1, cpv=cpg limit)
        Nu_star = Nu0
        do iter = 1, 5
            phi = cpv_eff / cpg * Sh_star / Nu_star / Le  ! eq. 22
            BT  = (1._R8 + BM)**phi - 1._R8               ! eq. 23
            if (BT <= 0._R8) then
                BT = BM  ! fallback
                Nu_star = Nu0
                exit
            endif
            FT = F_correction(BT)
            Nu_star = 2._R8 + (Nu0 - 2._R8) / FT         ! eq. 11
        enddo

        !> Step 3: gas-side heat Q_G [eq. 24, without L], positive into the droplet; the F(7)
        !  assembly adds the latent sink mdot*Lv itself, so Lv must NOT appear here.
        if (BT > 0._R8) then
            Qdot_evap = -mdot * cpv_eff * (Tg - Tp) / BT
        else
            Qdot_evap = 0._R8
        endif

    end subroutine ASM_model


    !===================================================================!
    !  Model 5: TC — Tonini-Cossali analytical model (2012)            !
    !  Variable-density Stefan-Fuchs: molar Stefan flow through the    !
    !  exact quiescent T(r); single monotone transcendental for m̂      !
    !  Ref: Tonini-Cossali IJTS 57 (2012) 45; Antonov et al. IJMF 179  !
    !  (2024) 104922 eq.(9); derivation: plan-bucket/f2-tc-derivation  !
    !===================================================================!
    pure subroutine TC_model(rhog, Tg, kg, Tp, dp, Mg, cpg, Pr, Sc, Re05, Ys, &
                             Le, cpv, Yinf, Mv, mdot, Qdot_evap, override_Qdot)
        implicit none
        real(R8), intent(in)  :: rhog, Tg, kg, Tp, dp, Mg, cpg, Pr, Sc, Re05, Ys
        real(R8), intent(in)  :: Le, cpv, Yinf, Mv
        real(R8), intent(out) :: mdot, Qdot_evap
        logical,  intent(out) :: override_Qdot
        integer,  parameter :: maxNewton = 30
        real(R8), parameter :: tolG = 1.e-13_R8
        real(R8) :: Dv, cpv_eff, Xs, Xinf, Minf, rhs0, Tts, Lev
        real(R8) :: mhat, mlo, mhi, x, E, f, fp, G, dG
        real(R8) :: Sh, Nu, chi
        integer  :: iter

        mdot = 0._R8
        Qdot_evap = 0._R8
        override_Qdot = .false.

        Dv = kg / (rhog * cpg * Le)
        if (cpv > 0._R8) then
            cpv_eff = cpv
        else
            cpv_eff = cpg  ! approximation (as ASM)
        endif

        !> Molar (partial-pressure) frame; cap Xs below 1 so the boiling clamp stays finite
        Xs   = min(mass2molar(Ys, Mv, Mg), 1._R8 - xsCap)
        Xinf = mass2molar(Yinf, Mv, Mg)
        Minf = Xinf*Mv + (1._R8 - Xinf)*Mg
        rhs0 = Mv/Minf * log((1._R8 - Xinf)/(1._R8 - Xs))  ! = -p_cr·ln[(p_cr-p_vs)/(p_cr-Yinf)]
        if (rhs0 <= 0._R8) return

        !> Solve m̂ + (T̃s-1)·Lev·(f(m̂/Lev)-1) = rhs0, f(x)=x/(1-e^-x); G'>0 for all T̃s>0
        Tts = Tp / Tg
        Lev = kg / (cpv_eff * Dv * rhog)          ! vapor-cp Lewis number (= Le·cpg/cpv)
        mlo = rhs0 / max(Tts, 1._R8)
        mhi = rhs0 / min(max(Tts, 1.e-3_R8), 1._R8)
        mhat = rhs0                               ! isothermal exact solution as first guess
        do iter = 1, maxNewton
            x = mhat / Lev
            if (x > 1.e-6_R8) then
                E  = exp(-x)
                f  = x / (1._R8 - E)
                fp = (1._R8 - E - x*E) / (1._R8 - E)**2
            else
                f  = 1._R8 + 0.5_R8*x + x*x/12._R8
                fp = 0.5_R8 + x/6._R8
            endif
            G  = mhat + (Tts - 1._R8)*Lev*(f - 1._R8) - rhs0
            dG = 1._R8 + (Tts - 1._R8)*fp
            if (abs(G) <= tolG*max(1._R8, rhs0)) exit
            if (G > 0._R8) then; mhi = mhat; else; mlo = mhat; endif
            mhat = mhat - G/dG
            if (mhat <= mlo .or. mhat >= mhi) mhat = 0.5_R8*(mlo + mhi)  ! bisection safeguard
        enddo

        !> Ranz-Marshall bolt-on (TC2012 itself is quiescent Sh=Nu=2; engineering extension)
        Sh = 2._R8 + 0.6_R8 * Re05 * Sc**(1._R8/3._R8)
        mdot = -pi * dp * rhog * Dv * Sh * mhat   ! = (Sh/2)·4π·Rd·Dv·ρ∞·m̂, ≤ 0

        !> Gas-side heat with the model's own BT = e^χ-1 (from the T(r) first integral);
        !  latent sink added by the F(7) assembly, so Lv must NOT appear here (as ASM).
        override_Qdot = .true.
        Nu  = 2._R8 + 0.6_R8 * Re05 * Pr**(1._R8/3._R8)
        chi = min(-mdot * cpv_eff / (pi * dp * kg * Nu), 700._R8)
        if (chi > 1.e-9_R8) then
            Qdot_evap = -mdot * cpv_eff * (Tg - Tp) / (exp(chi) - 1._R8)
        else
            Qdot_evap = pi * dp * kg * Nu * (Tg - Tp)  ! χ→0 conduction limit
        endif

    end subroutine TC_model


    !===================================================================!
    !  F(B) correction factor [Abramzon-Sirignano eq. 17]              !
    !  F(B) = (1+B)^0.7 * ln(1+B) / B                                 !
    !===================================================================!
    pure function F_correction(B) result(F)
        implicit none
        real(R8), intent(in) :: B
        real(R8) :: F

        if (B > 1.0e-10_R8) then
            F = (1._R8 + B)**0.7_R8 * log(1._R8 + B) / B
        else
            F = 1._R8  ! limit for B→0
        endif

    end function F_correction


    !> ln(1+x): the series below |x| = 1e-2 (truncation <= x^8/9 relative).
    pure function log1pS(x) result(y)
        implicit none
        real(R8), intent(in) :: x
        real(R8) :: y

        if (abs(x) >= 1.e-2_R8) then
            y = log(1._R8 + x)
        else
            y = x*(1._R8 - x*(0.5_R8 - x*(1._R8/3._R8 - x*(0.25_R8 - x*(0.2_R8 - x*(1._R8/6._R8 &
                  - x*(1._R8/7._R8 - x*0.125_R8)))))))
        endif

    end function log1pS


    !> e^y - 1: the series below |y| = 1e-2 (truncation <= y^8/9! relative).
    pure function expm1S(y) result(e)
        implicit none
        real(R8), intent(in) :: y
        real(R8) :: e

        if (abs(y) >= 1.e-2_R8) then
            e = exp(y) - 1._R8
        else
            e = y*(1._R8 + y*(0.5_R8 + y*(1._R8/6._R8 + y*(1._R8/24._R8 + y*(1._R8/120._R8 &
                  + y*(1._R8/720._R8 + y*(1._R8/5040._R8 + y/40320._R8)))))))
        endif

    end function expm1S


    !> Effective-conductivity factor chi(Pe_L) = 1.86 + 0.86 tanh(2.245 log10(Pe_L/30)) (Abramzon-Sirignano 1989 eq. 39);
    !  1 with no internal circulation (Pe_L <= 1e-300).
    pure function chiETC(PeL) result(chi)
        implicit none
        real(R8), intent(in) :: PeL
        real(R8) :: chi

        if (PeL <= 1.e-300_R8) then
            chi = 1._R8
        else
            chi = 1.86_R8 + 0.86_R8*tanh(2.245_R8*log10(PeL/30._R8))
        endif

    end function chiETC


    !> Boiling temperature at the pressure p on the Clausius-Clapeyron line: psat_CC(TB) = p.
    pure function boilingTemperatureCC(p, LvMvOverRu, invTboil) result(TB)
        implicit none
        real(R8), intent(in) :: p, LvMvOverRu, invTboil
        real(R8) :: TB

        TB = 1._R8/(invTboil - log(p/Patm)/LvMvOverRu)

    end function boilingTemperatureCC


    !> Boiling heat-transfer coefficient alpha_s [W/m2/K] at the superheat dT, floored at 0.5 K (Adachi et al. 1997,
    !  through Zuo-Gomes-Rutland 2000 eqs. 8-9 and OpenFOAM's LiquidEvaporationBoil).
    pure function alphaBoil(dT) result(alpha)
        implicit none
        real(R8), intent(in) :: dT
        real(R8) :: alpha, dTf

        dTf = max(dT, 0.5_R8)
        if (dTf < 5._R8) then
            alpha = 760._R8*dTf**0.26_R8
        elseif (dTf < 25._R8) then
            alpha = 27._R8*dTf**2.33_R8
        else
            alpha = 13800._R8*dTf**0.39_R8
        endif

    end function alphaBoil


    !> Zuo-Gomes-Rutland boiling of a drop at Tbar above its boiling temperature TB: the flash rate Gf fed by the
    !  superheat (alpha_s, or the parabolic conduction 10 pi kCond d dT with kCond), the heat-fed rate G of
    !  G (1 + Gf/G) = B ln(1 + A (1 + Gf/G)), A = cpg (Tg - TB)/Lv, B = pi kf d Nu/cpg with the 1/3-rule film
    !  conductivity; mdot = -(G + Gf) and Qdot [W] the heat through the blown film.
    pure subroutine zgrBoiling(Tbar, TB, d, kg, cpg, Nu, Tg, Lv, mdot, Qdot, kCond)
        implicit none
        real(R8), intent(in)  :: Tbar, TB, d, kg, cpg, Nu, Tg, Lv
        real(R8), intent(out) :: mdot, Qdot
        real(R8), intent(in), optional :: kCond
        real(R8) :: dTf, Gf, A, Tf, kf, B, G, S, Qs

        dTf = max(Tbar - TB, 0.5_R8)
        if (present(kCond)) then
            Gf = 10._R8*pi*kCond*d*dTf/Lv
        else
            Gf = alphaBoil(Tbar - TB)*dTf*pi*d*d/Lv
        endif
        A  = cpg*(Tg - TB)/Lv
        Tf = TB + (Tg - TB)/3._R8
        kf = kg*(Tf/Tg)**0.7_R8
        B  = pi*kf*d*Nu/cpg
        if (A <= 0._R8) then
            G = 0._R8
        elseif (Gf == 0._R8) then
            G = B*log1pS(A)
        else
            G = heatFedRate(Gf, A, B)
        endif
        S = G + Gf
        if (S/B > 0._R8) then
            Qs = S*cpg*(Tg - TB)/expm1S(S/B)
        else
            Qs = B*cpg*(Tg - TB)
        endif
        mdot = -S
        Qdot = Qs

    end subroutine zgrBoiling


    !> The heat-fed rate G = Gf/x of zgrBoiling: x > 0 the root of H(x) = Gf (1 + x)/x - B ln(1 + A (1 + x)), strictly
    !  decreasing, by Brent in u = ln x on [-60, 60] to 4 eps; G = B ln(1 + A) or 0 when the root lies beyond an end.
    pure function heatFedRate(Gf, A, B) result(G)
        implicit none
        real(R8), intent(in) :: Gf, A, B
        real(R8) :: G
        integer,  parameter :: maxEval = 200
        real(R8), parameter :: uLo = -60._R8, uHi = 60._R8, tol4 = 4._R8*epsilon(1._R8)
        real(R8) :: ua, ub, uc, fa, fb, fc, dd, e, s, p, q, r, tol1, xm
        integer  :: n

        ua = uLo; fa = hRes(ua)
        ub = uHi; fb = hRes(ub)
        if (fa <= 0._R8) then; G = B*log1pS(A); return; endif   ! Gf negligible: the unblown heat-fed rate
        if (fb >= 0._R8) then; G = 0._R8; return; endif        ! Gf dominant: no heat-fed vapour
        uc = ua; fc = fa; dd = ub - ua; e = dd
        do n = 1, maxEval
            if (abs(fc) < abs(fb)) then
                ua = ub; ub = uc; uc = ua; fa = fb; fb = fc; fc = fa
            endif
            tol1 = 0.5_R8*tol4*(abs(ub) + 1._R8)
            xm   = 0.5_R8*(uc - ub)
            if (abs(xm) <= tol1 .or. fb == 0._R8) exit
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                s = fb/fa
                if (ua == uc) then
                    p = 2._R8*xm*s; q = 1._R8 - s
                else
                    q = fa/fc; r = fb/fc
                    p = s*(2._R8*xm*q*(q - r) - (ub - ua)*(r - 1._R8))
                    q = (q - 1._R8)*(r - 1._R8)*(s - 1._R8)
                endif
                if (p > 0._R8) then; q = -q; else; p = -p; endif
                s = e; e = dd
                if (2._R8*p >= 3._R8*xm*q - abs(tol1*q) .or. p >= abs(0.5_R8*s*q)) then
                    dd = xm; e = dd
                else
                    dd = p/q
                endif
            else
                dd = xm; e = dd
            endif
            ua = ub; fa = fb
            if (abs(dd) > tol1) then; ub = ub + dd; else; ub = ub + sign(tol1, xm); endif
            fb = hRes(ub)
            if ((fb > 0._R8 .and. fc > 0._R8) .or. (fb < 0._R8 .and. fc < 0._R8)) then
                uc = ua; fc = fa; dd = ub - ua; e = dd
            endif
        enddo
        G = Gf*exp(-ub)

    contains

        !> H at u = ln x, divided by B.
        pure function hRes(u) result(h)
            real(R8), intent(in) :: u
            real(R8) :: h
            h = (Gf/B)*(1._R8 + exp(-u)) - log1pS(A*(1._R8 + exp(u)))
        end function hRes

    end function heatFedRate


    !> 1/T of the one-atmosphere boiling point on the Clausius-Clapeyron line through (Tboil, pBoil).
    pure function invTboil_at(Tboil, pBoil, LvMvOverRu) result(invT)
        implicit none
        real(R8), intent(in) :: Tboil, pBoil, LvMvOverRu
        real(R8) :: invT

        invT = 1._R8/Tboil + log(pBoil/Patm)/LvMvOverRu

    end function invTboil_at


    !===================================================================!
    !  Clausius-Clapeyron saturation pressure                          !
    !  psat(T) = Patm * exp(-LvMv/Ru * (1/T - 1/Tboil))              !
    !===================================================================!
    pure function psat_CC(T, LvMvOverRu, invTboil) result(psat)
        implicit none
        real(R8), intent(in) :: T, LvMvOverRu, invTboil
        real(R8) :: psat

        psat = Patm * exp(-LvMvOverRu * (1._R8/T - invTboil))

    end function psat_CC


end module ICE_Lib_Evaporation