!> Unit test of liquid-conduction = P2T, the parabolic temperature profile inside the drop, through the cell
!> source (ICE_Mod_Sources' compute_source_ and cell_surface, ICE_Lib_Drop_Surface's p2t_surface,
!> ICE_Lib_Evaporation's chiETC), against closed forms and identities typed here. One MK cell, the drops at
!> rest and the gas moving at the slip; IGLOO's gas (rho 1.2, mu 1.8e-5, k 0.026, gamma 1.4, R 287) and vapour
!> (Lv 2.26e6, Mv 18, T_boil 373.15, cp_v 1900, Le 1); drag Stokes, heat Ranz-Marshall unless said. S1 is the
!> mass slot, S5 the energy slot, F7 = (S5 - c T S1)/(rho_p c) the mean temperature's rate.
!>   L1  Re 0, no mass transfer, zeta = Nu k_g/(2 k_l) = 1: S5 = n pi d k_g Nu 5 (T_g - T)/(5 + zeta), S1 = 0
!>   L2  L1's surface temperature (5 T + zeta T_g)/(5 + zeta) = 350 K
!>   L3  k-liq 1e12: the sources of the infinite-conductivity (ITC) material
!>   L4  chi(Pe_L) = 1.86 + 0.86 tanh(2.245 log10(Pe_L/30)) [AS89] eq. 39, 1 without circulation
!>   L5  CEM at Re 20: the surface balance, S1 = n mdot(Ts), F7 = (Q_g + mdot (Lv + c (Ts - T)))/(m c), with Q_g and
!>       mdot typed here at Ts (the vapour leaves at the surface temperature)
!>   L6  a grid of states and gas-side models, with the boiling clamp and with boiling = ZGR: every slot and Ts
!>       finite, every surface balance solved
!>   L7  informational: chi at T and at Ts, Ts - T at L5's state
!>   L8  boiling = ZGR with P2T at A = 0 (T_g = T_B = 600 K), water: Ts = T_B, the flash fed by the parabolic
!>       conduction F8 = -10 pi k_l d dT_f/Lv, and (F8/F7)/m = c/(Lv - c dT) (the vapour leaves at T_B)
!>   L9  L1 with emissivity 0.8: the radiation inside the surface balance and in the source
!>   L10 drag Henderson at Re 1, Ma 0.26: the drag's temperature ratio taken at Ts
!>   L11 cp, Lv and Psat tabulated: S1 = n mdot(Ts), S5 = n (Q_g(Ts) + mdot(Ts) (e(Ts) + Lv(Ts)))
!>   L13 an empty cell (rho_p <= rho_empty), P2T and ZGR: Ts = T and the ITC sources, bit for bit
!>   L14 a mean temperature at the 1e-25 K floor, no bracket: p2t_surface's stat 0, Ts and the energy slot not finite
!>   L12 L5 with a Conductivity column constant at 0.6 and k-liq 0.3 beside it: L5's sources bit for bit
!>   L15 L8 with a Conductivity column constant at 0.6 and k-liq 0: L8's sources bit for bit (the flash's k_l)
!>   L16 mu-liq = 0 at Re 20 (the inputs refuse it): p2t_surface's stat 0 and Ts NaN, the source's Ts NaN
!>   L17 a Viscosity column mu_l(T) = 1e-3 exp((320 - T)/20) at L5's state: Ts the root of the surface balance with
!>       mu_l at the surface temperature [AS89] (U_s eqs. 30, 32 and chi eq. 39 typed here, the column interpolated
!>       here, B_M, k_l and c_l at T), by bisection; the root with mu_l frozen at T must differ
!> Every leg checks that cell_surface's Ts equals compute_source_'s bit for bit (the probe reads the former).
!> Prints every assertion and exits non-zero if any failed.
program test_liquid_p2t
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_next_after
  implicit none
  real(R8), parameter :: PI = 3.14159265358979323846_R8, RU = 8314.46_R8, PATM = 101325._R8
  real(R8), parameter :: SIGMA_SB = 5.670374419e-8_R8
  real(R8), parameter :: RHOG = 1.2_R8, MUG = 1.8e-5_R8, KG = 0.026_R8, GAM = 1.4_R8, RG = 287._R8
  real(R8), parameter :: LV = 2.26e6_R8, MV = 18._R8, TBOIL = 373.15_R8, CPV = 1900._R8
  real(R8), parameter :: CPG = GAM*RG/(GAM - 1._R8), PRG = MUG*CPG/KG
  logical :: ok

  ok = .true.
  call setup_model('Stokes', 'Ranz-Marshall')
  call legs_rest(ok)
  call leg_chi(ok)
  call leg_cem(ok)
  call leg_henderson(ok)
  call setup_model('Stokes', 'Ranz-Marshall')
  call leg_tables(ok)
  call leg_empty(ok)
  call leg_floor(ok)
  call leg_no_viscosity(ok)
  call leg_flash(ok)
  call leg_grid(ok)
  call leg_columns(ok)
  call leg_viscosity(ok)
  if (.not. ok) then
    write(*,'(A)') 'test_liquid_p2t: FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_liquid_p2t: PASS'

contains

  subroutine check(cond, what, ok)
    logical,          intent(in)    :: cond
    character(len=*), intent(in)    :: what
    logical,          intent(inout) :: ok
    if (cond) then
      write(*,'(A)') '   OK   '//what
    else
      write(*,'(A)') '   FAIL '//what
      ok = .false.
    endif
  end subroutine check


  !> Family 1 is MK with the given drag and heat laws, no blowing correction.
  subroutine setup_model(drag, heat)
    use ICE_Config_Types_m, only: obj_time_scheme
    use ICE_Global_m,       only: solid_of
    use ICE_Lib_Model,      only: assign_all
    use ICE_Lib_Drag,       only: assign_drag
    use ICE_Lib_Heat,       only: assign_heat
    character(len=*), intent(in) :: drag, heat
    if (.not. allocated(obj_time_scheme%model)) then
      allocate(obj_time_scheme%model(1))
      obj_time_scheme%model(1) = 'MK'
      allocate(solid_of(1))
      solid_of = .false.
      call assign_all(1)
    endif
    call assign_drag(drag, obj_time_scheme%dragSelect)
    call assign_heat(heat, obj_time_scheme%heatSelect)
    obj_time_scheme%blowSelect = 0
  end subroutine setup_model


  !> A material of density rho and cp c with the evaporation model word, Y_inf yinf, P2T with k-liq kl and mu-liq mul.
  subroutine material(mat, rho, c, kl, mul, word, yinf)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: assign_evaporation
    type(condensed_phase_t), intent(out) :: mat
    real(R8),                intent(in)  :: rho, c, kl, mul, yinf
    character(len=*),        intent(in)  :: word
    call assign_evaporation(word, mat%evapSelect)
    mat%intfSelect = 0; mat%liqSelect = 1; mat%boilSelect = 0
    mat%rho_al = rho; mat%cs_al = c; mat%emiss = 0._R8; mat%kLiq = kl; mat%muLiq = mul
    mat%lv_al = LV; mat%Mv = MV; mat%cpv = CPV; mat%Le = 1._R8; mat%Yinf = yinf; mat%Tboil = TBOIL; mat%alphaE = 1._R8
    call pack(mat)
  end subroutine material


  !> The packed vapour data, as Setup_Materials builds them.
  subroutine pack(mat)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: iMv, iLv, icpv, iLe, iYinf, iLvMvOverRu, iinvTboil, ialphaE
    type(condensed_phase_t), intent(inout) :: mat
    mat%ep = 0._R8
    mat%ep(iMv)         = mat%Mv
    mat%ep(iLv)         = mat%lv_al
    mat%ep(icpv)        = mat%cpv
    mat%ep(iLe)         = mat%Le
    mat%ep(iYinf)       = mat%Yinf
    mat%ep(iLvMvOverRu) = mat%lv_al*mat%Mv/RU
    mat%ep(iinvTboil)   = 1._R8/mat%Tboil
    mat%ep(ialphaE)     = mat%alphaE
  end subroutine pack


  !> The cell's sources and surface temperature; cell_surface's Ts must be compute_source_'s, bit for bit.
  subroutine sources(mat, cond, gas, s, ts, name, ok, Rp, Re, Nu)
    use ICE_Config_Types_m, only: condensed_phase_t
    use ICE_Mod_Sources,    only: compute_source_, cell_surface
    type(condensed_phase_t), intent(in)    :: mat
    real(R8),                intent(in)    :: cond(6), gas(5)
    real(R8),                intent(out)   :: s(6), ts
    character(len=*),        intent(in)    :: name
    logical,                 intent(inout) :: ok
    real(R8), optional,      intent(out)   :: Rp, Re, Nu
    real(R8) :: tau, rho_mat, rp_, re_, ma_, pr_, nu_, tsc, tb_
    logical  :: boil_
    tau = 0._R8
    call compute_source_(cond, tau, gas, RG, GAM, KG, MUG, s, mat, 6, ts)
    call cell_surface(cond, gas, RG, GAM, KG, MUG, mat, 6, rho_mat, rp_, re_, ma_, pr_, nu_, tsc, tb_, boil_)
    call check(tsc == ts .or. (ieee_is_nan(tsc) .and. ieee_is_nan(ts)), &
               name//": cell_surface's Ts is the source's, bit for bit", ok)
    if (present(Rp)) Rp = rp_
    if (present(Re)) Re = re_
    if (present(Nu)) Nu = nu_
  end subroutine sources


  !> The diameter of the cell's drops, as ICE forms it.
  pure function diameter(cond, rho) result(d)
    real(R8), intent(in) :: cond(6), rho
    real(R8) :: d
    d = 2._R8*(0.75_R8*cond(1)/(cond(6)*PI*rho))**(1._R8/3._R8)
  end function diameter


  !> L1, L2, L3, L9: Re 0, Y_inf = 1 (no mass transfer), T 300 K, T_g 600 K, d 11.89 um, rho 2950, c 1250.
  subroutine legs_rest(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    logical, intent(inout) :: ok
    real(R8), parameter :: T0 = 300._R8, TG = 600._R8, RHO = 2950._R8, C = 1250._R8, KL = 0.026_R8
    type(condensed_phase_t) :: mat, itc
    real(R8) :: cond(6), gas(5), s(6), si(6), ts, tsi, d, zeta, ref, f7, rel, nu, a, f, fp, dt
    integer  :: it

    cond = [2.596367726746e-3_R8, 0._R8, 0._R8, 0._R8, T0, 1.e9_R8]
    gas  = [RHOG, 0._R8, 0._R8, 0._R8, TG]
    call material(mat, RHO, C, KL, 1.e-3_R8, 'CEM', 1._R8)
    call sources(mat, cond, gas, s, ts, 'L1', ok, Nu=nu)
    d    = diameter(cond, RHO)
    zeta = nu*KG/(2._R8*KL)
    ref  = cond(6)*PI*d*KG*nu*5._R8*(TG - T0)/(5._R8 + zeta)
    rel  = abs(s(5)/ref - 1._R8)
    f7   = (s(5) - C*T0*s(1))/(cond(1)*C)
    write(*,'(A,ES22.14,A,ES22.14,A,ES10.3,A,F10.6)') '   L1: S5 ', s(5), ' closed form ', ref, ' rel ', rel, ' zeta ', zeta
    write(*,'(A,ES22.14)') '   L1: F7 ', f7
    call check(s(1) == 0._R8, 'L1 no mass transfer: S1 = 0 exactly', ok)
    call check(rel <= 1.e-13_R8, 'L1 S5 = n pi d k_g Nu 5 (T_g - T)/(5 + zeta) to 1e-13', ok)
    ref = (5._R8*T0 + zeta*TG)/(5._R8 + zeta)
    rel = abs(ts/ref - 1._R8)
    write(*,'(A,F20.12,A,F20.12,A,ES10.3)') '   L2: Ts ', ts, ' closed form ', ref, ' rel ', rel
    call check(rel <= 1.e-13_R8, 'L2 Ts = (5 T + zeta T_g)/(5 + zeta) to 1e-13', ok)

    ! L3: the infinite-conductivity limit
    mat%kLiq = 1.e12_R8
    call sources(mat, cond, gas, s, ts, 'L3', ok)
    call material(itc, RHO, C, KL, 1.e-3_R8, 'CEM', 1._R8)
    itc%liqSelect = 0
    call sources(itc, cond, gas, si, tsi, 'L3 (ITC)', ok)
    rel = maxval(abs(s - si)/max(abs(si), tiny(1._R8)), mask=si /= 0._R8)
    write(*,'(A,ES12.4,A,ES10.3)') '   L3: Ts - T ', ts - T0, ', worst relative slot difference ', rel
    call check(all((s == 0._R8) .eqv. (si == 0._R8)) .and. rel <= 1.e-12_R8, 'L3 k-liq 1e12 gives the ITC sources to 1e-12', ok)

    ! L9: emissivity 0.8 inside the surface balance (the root by Newton, typed here) and in the source
    call material(mat, RHO, C, KL, 1.e-3_R8, 'CEM', 1._R8)
    mat%emiss = 0.8_R8
    call sources(mat, cond, gas, s, ts, 'L9', ok, Nu=nu)
    a = 350._R8
    do it = 1, 50
      f  = 10._R8*PI*KL*d*(a - T0) - PI*d*KG*nu*(TG - a) - 0.8_R8*SIGMA_SB*PI*d*d*(TG**4 - a**4)
      fp = 10._R8*PI*KL*d + PI*d*KG*nu + 4._R8*0.8_R8*SIGMA_SB*PI*d*d*a**3
      dt = f/fp
      a  = a - dt
      if (abs(dt) <= 1.e-15_R8*a) exit
    enddo
    ref = cond(6)*(PI*d*KG*nu*(TG - a) + 0.8_R8*SIGMA_SB*PI*d*d*(TG**4 - a**4))
    write(*,'(A,F20.12,A,F20.12,A,ES22.14,A,ES22.14)') '   L9: Ts ', ts, ' root ', a, ' S5 ', s(5), ' oracle ', ref
    write(*,'(A,ES10.3,A,ES10.3)') '   L9: relative: Ts ', abs(ts/a - 1._R8), ' S5 ', abs(s(5)/ref - 1._R8)
    call check(abs(ts/a - 1._R8) <= 1.e-13_R8, 'L9 Ts the root of the balance with radiation to 1e-13', ok)
    call check(abs(s(5)/ref - 1._R8) <= 1.e-13_R8, 'L9 S5 = n (Q_g(Ts) + eps sigma pi d^2 (T_g^4 - Ts^4)) to 1e-13', ok)
  end subroutine legs_rest


  !> L4: chiETC.
  subroutine leg_chi(ok)
    use ICE_Lib_Evaporation, only: chiETC
    logical, intent(inout) :: ok
    real(R8) :: c(4)
    c = [chiETC(0._R8), chiETC(30._R8), chiETC(300._R8), chiETC(1.e12_R8)]
    write(*,'(A,4F20.15)') '   L4: chi ', c
    call check(all(abs(c - [1._R8, 1.86_R8, 2.70091464264437_R8, 2.72_R8]) <= 1.e-14_R8), &
               'L4 chi(Pe_L) at 0, 30, 300, 1e12: 1, 1.86, 2.70091464264437, 2.72 to 1e-14', ok)
  end subroutine leg_chi


  !> The CEM rate [kg/s] at the surface temperature T (Clausius-Clapeyron, or the given saturation pressure).
  function cem_rate(T, Tg, d, Re, mv_, ps) result(md)
    real(R8), intent(in) :: T, Tg, d, Re, mv_
    real(R8), intent(in), optional :: ps
    real(R8) :: md, p, psat, xs, ys, bm, sh
    p = RHOG*RG*Tg
    if (present(ps)) then
      psat = ps
    else
      psat = PATM*exp(-LV*mv_/RU*(1._R8/T - 1._R8/TBOIL))
    endif
    xs = min(psat/p, 1._R8 - 1.e-12_R8)
    ys = xs*mv_/(xs*mv_ + (1._R8 - xs)*RU/RG)
    bm = ys/(1._R8 - ys)
    sh = 2._R8 + 0.6_R8*sqrt(Re)*PRG**(1._R8/3._R8)
    md = -PI*d*RHOG*(KG/(RHOG*CPG))*sh*log(1._R8 + bm)
  end function cem_rate


  !> L5, L7: CEM, T 320 K, T_g 800 K, u_g 6 m/s (Re 20), d 50 um, rho 1000, c 4180, k-liq 0.6, mu-liq 1e-3, Y_inf 0.
  subroutine leg_cem(ok)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: chiETC
    logical, intent(inout) :: ok
    real(R8), parameter :: T0 = 320._R8, TG = 800._R8, UG = 6._R8, RHO = 1000._R8, C = 4180._R8, KL = 0.6_R8, MUL = 1.e-3_R8
    real(R8), parameter :: D0 = 50.e-6_R8
    type(condensed_phase_t) :: mat
    real(R8) :: cond(6), gas(5), s(6), ts, d, re, nu, m, md, qg, k, bm, us, chi0, chis, res, f7, ref, xs, ys
    cond = [1._R8, 0._R8, 0._R8, 0._R8, T0, 1._R8/(RHO*PI*D0**3/6._R8)]
    gas  = [RHOG, UG, 0._R8, 0._R8, TG]
    call material(mat, RHO, C, KL, MUL, 'CEM', 0._R8)
    call sources(mat, cond, gas, s, ts, 'L5', ok, Re=re, Nu=nu)
    d = diameter(cond, RHO)
    m = cond(1)/cond(6)
    ! B_M at the mean temperature, K formed here (A-S eqs. 30, 32, 39)
    xs   = PATM*exp(-LV*MV/RU*(1._R8/T0 - 1._R8/TBOIL))/(RHOG*RG*TG)
    ys   = xs*MV/(xs*MV + (1._R8 - xs)*RU/RG)
    bm   = ys/(1._R8 - ys)
    us   = (12.69_R8/32._R8)*UG*(MUG/MUL)*re**(1._R8/3._R8)/(1._R8 + bm)
    chi0 = 1.86_R8 + 0.86_R8*tanh(2.245_R8*log10(us*d*RHO*C/KL/30._R8))
    k    = 10._R8*PI*chi0*KL*d
    md   = cem_rate(ts, TG, d, re, MV)
    qg   = PI*d*KG*nu*(TG - ts)
    res  = abs(ts - T0 - (qg + md*LV)/k)
    f7   = (s(5) - C*T0*s(1))/(cond(1)*C)
    ref  = (qg + md*(LV + C*(ts - T0)))/(m*C)
    write(*,'(A,F10.6,A,F14.10,A,F20.12,A,ES10.3)') '   L5: Re ', re, ' Nu ', nu, ' Ts ', ts, ' residual [K] ', res
    write(*,'(A,ES22.14,A,ES22.14,A,ES22.14,A,ES22.14)') '   L5: S1 ', s(1), ' n mdot(Ts) ', cond(6)*md, &
      ' F7 ', f7, ' oracle ', ref
    write(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') '   L5: relative: Ts ', abs(ts/323.682316697455_R8 - 1._R8), &
      ' S1 ', abs(s(1)/(cond(6)*md) - 1._R8), ' F7 ', abs(f7/ref - 1._R8)
    call check(res <= 4._R8*epsilon(1._R8)*T0 .and. T0 < ts .and. ts < TG, &
               'L5 control: |Ts - T - Q_net(Ts)/K| <= 4 eps T, T < Ts < T_g', ok)
    call check(abs(ts/323.682316697455_R8 - 1._R8) <= 1.e-12_R8, 'L5 Ts = 323.682316697455 K (mpmath) to 1e-12', ok)
    call check(abs(s(1)/(cond(6)*md) - 1._R8) <= 1.e-13_R8, 'L5 S1 = n mdot(Ts) to 1e-13', ok)
    call check(abs(f7/ref - 1._R8) <= 1.e-13_R8, 'L5 F7 = (Q_g + mdot (Lv + c (Ts - T)))/(m c) to 1e-13', ok)
    ! L7: informational, chi with B_M at Ts
    xs   = PATM*exp(-LV*MV/RU*(1._R8/ts - 1._R8/TBOIL))/(RHOG*RG*TG)
    ys   = xs*MV/(xs*MV + (1._R8 - xs)*RU/RG)
    chis = chiETC((12.69_R8/32._R8)*UG*(MUG/MUL)*re**(1._R8/3._R8)/(1._R8 + ys/(1._R8 - ys))*d*RHO*C/KL)
    write(*,'(A,F12.8,A,F12.8,A,F12.8)') '   L7: chi(T) ', chi0, ', chi(Ts) ', chis, ', Ts - T ', ts - T0
  end subroutine leg_cem


  !> Henderson's subsonic drag, typed from test/verification/common.py.
  pure function henderson(Re, Ma, G, Tr) result(cd)
    real(R8), intent(in) :: Re, Ma, G, Tr
    real(R8), parameter :: TOLL = 1.e-20_R8
    real(R8) :: cd, s
    s  = Ma*(0.5_R8*G)**0.5_R8
    cd = 24._R8*(Re + s*(4.33_R8 + (3.65_R8 - 1.53_R8*Tr)/(1._R8 + 0.353_R8*Tr)*exp(-0.247_R8*Re/(s + TOLL))) &
         + TOLL)**(-1._R8) + exp(-0.5_R8*Ma/(Re**0.5_R8 + TOLL))*((4.5_R8 + 0.38_R8*(0.03_R8*Re + 0.48_R8*Re**0.5_R8)) &
         /(1._R8 + 0.03_R8*Re + 0.48_R8*Re**0.5_R8) + 0.1_R8*Ma**2 + 0.2_R8*Ma**8) + 0.6_R8*s*(1._R8 - exp(-Ma/(Re + TOLL)))
  end function henderson


  !> L10: drag Henderson, CEM, T 320 K, T_g 800 K, u_g 150 m/s, d 100 nm.
  subroutine leg_henderson(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    logical, intent(inout) :: ok
    real(R8), parameter :: T0 = 320._R8, TG = 800._R8, UG = 150._R8, RHO = 1000._R8, D0 = 100.e-9_R8
    type(condensed_phase_t) :: mat
    real(R8) :: cond(6), gas(5), s(6), ts, rp, re, ma, cd, tau, ref
    call setup_model('Henderson', 'Ranz-Marshall')
    cond = [1._R8, 0._R8, 0._R8, 0._R8, T0, 1._R8/(RHO*PI*D0**3/6._R8)]
    gas  = [RHOG, UG, 0._R8, 0._R8, TG]
    call material(mat, RHO, 4180._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    call sources(mat, cond, gas, s, ts, 'L10', ok, Rp=rp, Re=re)
    ma  = UG/sqrt(GAM*RG*TG)
    cd  = henderson(re, ma, GAM, ts/TG)
    tau = 8._R8*RHO*rp/(3._R8*RHOG*cd*UG + 1.e-20_R8)
    ref = cond(1)*UG/tau
    write(*,'(A,F14.10,A,F14.10,A,F20.12,A,F18.12,A,ES22.14,A,ES22.14)') '   L10: Re ', re, ' Ma ', ma, ' Ts ', ts, &
      ' Cd ', cd, ' source(2) ', s(2), ' oracle ', ref
    write(*,'(A,ES10.3)') '   L10: relative source(2) ', abs(s(2)/ref - 1._R8)
    call check(abs(ts/324.397176409744_R8 - 1._R8) <= 1.e-12_R8, 'L10 Ts = 324.397176409744 K (mpmath) to 1e-12', ok)
    call check(abs(s(2)/ref - 1._R8) <= 1.e-13_R8, 'L10 source(2) = rho_p u_g/tau with Cd(Ts/T_g) to 1e-13', ok)
  end subroutine leg_henderson


  !> L11: T 320.5 K, T_g 800 K, u_g 6 m/s, d 50 um, rho 1000; on 280..400 K cp = 4180 + 2 (T - 300) (e = 4180 T +
  !> (T - 300)^2), Lv = 2.26e6 - 2500 (T - 300), Psat exact for that Lv through 1 atm at 373.15 K, Mv 18.015; CEM.
  subroutine leg_tables(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    logical, intent(inout) :: ok
    real(R8), parameter :: T0 = 320.5_R8, TG = 800._R8, UG = 6._R8, RHO = 1000._R8, D0 = 50.e-6_R8, MV11 = 18.015_R8
    type(condensed_phase_t) :: mat
    real(R8) :: et(280:400), ct(280:400), lt(280:400), pt(280:400), cond(6), gas(5), s(6), ts, d, re, nu, md, ref
    real(R8) :: a, b
    integer  :: i
    a = 2.26e6_R8 + 2500._R8*300._R8; b = -2500._R8
    do i = 280, 400
      ct(i) = 4180._R8 + 2._R8*real(i - 300, R8)
      et(i) = 4180._R8*real(i, R8) + real(i - 300, R8)**2
      lt(i) = 2.26e6_R8 - 2500._R8*real(i - 300, R8)
      pt(i) = PATM*exp(MV11/RU*(a*(1._R8/TBOIL - 1._R8/real(i, R8)) + b*log(real(i, R8)/TBOIL)))
    enddo
    call material(mat, RHO, 4180._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    mat%Mv = MV11
    mat%T_min = 280; mat%T_max = 400; mat%use_table = .true.; mat%cs_varies = .true.
    allocate(mat%rho_tab(280:400), source=[(RHO, i = 280, 400)])
    allocate(mat%cs_tab(280:400), source=ct)
    allocate(mat%e_tab(280:400), source=et)
    allocate(mat%h_tab(280:400), source=et)
    allocate(mat%lv_tab(280:400), source=lt); mat%use_lv = .true.
    allocate(mat%psat_tab(280:400), source=pt); mat%use_psat = .true.
    mat%lv_al = value_at(lt, 298.15_R8)
    call pack(mat)
    cond = [1._R8, 0._R8, 0._R8, 0._R8, T0, 1._R8/(RHO*PI*D0**3/6._R8)]
    gas  = [RHOG, UG, 0._R8, 0._R8, TG]
    call sources(mat, cond, gas, s, ts, 'L11', ok, Re=re, Nu=nu)
    d   = diameter(cond, RHO)
    md  = cem_rate(ts, TG, d, re, MV11, value_at(pt, ts))
    ref = cond(6)*(PI*d*KG*nu*(TG - ts) + md*(value_at(et, ts) + value_at(lt, ts)))
    write(*,'(A,F20.12,A,ES22.14,A,ES22.14,A,ES22.14)') '   L11: Ts ', ts, ' S1 ', s(1), ' S5 ', s(5), ' oracle ', ref
    write(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') '   L11: relative: Ts ', abs(ts/324.101360088569_R8 - 1._R8), &
      ' S1 ', abs(s(1)/(cond(6)*md) - 1._R8), ' S5 ', abs(s(5)/ref - 1._R8)
    call check(abs(ts/324.101360088569_R8 - 1._R8) <= 1.e-12_R8, 'L11 Ts = 324.101360088569 K (mpmath) to 1e-12', ok)
    call check(abs(s(1)/(cond(6)*md) - 1._R8) <= 1.e-13_R8, 'L11 S1 = n mdot(Ts) to 1e-13', ok)
    call check(abs(s(5)/ref - 1._R8) <= 1.e-13_R8, 'L11 S5 = n (Q_g(Ts) + mdot(Ts) (e(Ts) + Lv(Ts))) to 1e-13', ok)
  end subroutine leg_tables


  !> L13: L5's state at rho_p 1e-7 kg/m^3 (n scaled, d 50 um): an empty cell keeps the mean temperature.
  subroutine leg_empty(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    logical, intent(inout) :: ok
    real(R8), parameter :: RHO = 1000._R8, D0 = 50.e-6_R8
    type(condensed_phase_t) :: mat, itc
    real(R8) :: cond(6), gas(5), s(6), si(6), ts, tsi
    cond = [1.e-7_R8, 0._R8, 0._R8, 0._R8, 320._R8, 1.e-7_R8/(RHO*PI*D0**3/6._R8)]
    gas  = [RHOG, 6._R8, 0._R8, 0._R8, 800._R8]
    call material(mat, RHO, 4180._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    mat%boilSelect = 1
    call material(itc, RHO, 4180._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    itc%liqSelect = 0
    call sources(mat, cond, gas, s, ts, 'L13', ok)
    call sources(itc, cond, gas, si, tsi, 'L13 (ITC)', ok)
    write(*,'(A,F20.12)') '   L13: Ts ', ts
    call check(ts == 320._R8, 'L13 empty cell: Ts = T bit for bit', ok)
    call check(all(s == si), 'L13 empty cell: the ITC sources bit for bit', ok)
  end subroutine leg_empty


  !> L14: a Psat column 101.325 T on 1..5000 K, NoHeat, CEM, Y_inf 0, d 50 um, T at the 1e-25 K floor: no bracket.
  subroutine leg_floor(ok)
    use ICE_Config_Types_m,   only: condensed_phase_t
    use ICE_Mod_Sources,      only: cell_surface
    use ICE_Lib_Drop_Surface, only: p2t_surface
    logical, intent(inout) :: ok
    real(R8), parameter :: RHO = 1000._R8, D0 = 50.e-6_R8
    type(condensed_phase_t) :: mat
    real(R8) :: cond(6), gas(5), s(6), ts, rho_mat, rp, re, ma, pr, nu, tsd, tb
    logical  :: boil
    integer  :: i, stat
    call setup_model('Stokes', 'NoHeat')
    call material(mat, RHO, 4180._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    mat%T_min = 1; mat%T_max = 5000; mat%use_table = .true.
    allocate(mat%rho_tab(1:5000), source=[(RHO, i = 1, 5000)])
    allocate(mat%cs_tab(1:5000), source=[(4180._R8, i = 1, 5000)])
    allocate(mat%psat_tab(1:5000), source=[(101.325_R8*real(i, R8), i = 1, 5000)]); mat%use_psat = .true.
    cond = [1._R8, 0._R8, 0._R8, 0._R8, 1.e-25_R8, 1._R8/(RHO*PI*D0**3/6._R8)]
    gas  = [RHOG, 0._R8, 0._R8, 0._R8, 800._R8]
    call cell_surface(cond, gas, RG, GAM, KG, MUG, mat, 6, rho_mat, rp, re, ma, pr, nu, tsd, tb, boil)
    call p2t_surface(cond(5), rho_mat, rp, re, nu, 0._R8, gas, RG, GAM, KG, MUG, mat, tsd, stat)
    call sources(mat, cond, gas, s, ts, 'L14', ok)
    write(*,'(A,I0,A,ES12.4,A,ES12.4)') '   L14: stat ', stat, ' Ts ', ts, ' S5 ', s(5)
    call check(stat == 0 .and. ieee_is_nan(tsd), "L14 no bracket: p2t_surface's stat 0 and Ts NaN", ok)
    call check(.not. ieee_is_finite(ts) .and. .not. ieee_is_finite(s(5)), &
               "L14 no bracket: the source's Ts and energy slot not finite", ok)
    call setup_model('Stokes', 'Ranz-Marshall')
  end subroutine leg_floor


  !> L16: L5's state with mu-liq 0: no liquid viscosity reaches the surface balance.
  subroutine leg_no_viscosity(ok)
    use ICE_Config_Types_m,   only: condensed_phase_t
    use ICE_Mod_Sources,      only: cell_surface
    use ICE_Lib_Drop_Surface, only: p2t_surface
    logical, intent(inout) :: ok
    real(R8), parameter :: RHO = 1000._R8, D0 = 50.e-6_R8
    type(condensed_phase_t) :: mat
    real(R8) :: cond(6), gas(5), s(6), ts, rho_mat, rp, re, ma, pr, nu, tsd, tb
    logical  :: boil
    integer  :: stat
    call material(mat, RHO, 4180._R8, 0.6_R8, 0._R8, 'CEM', 0._R8)
    cond = [1._R8, 0._R8, 0._R8, 0._R8, 320._R8, 1._R8/(RHO*PI*D0**3/6._R8)]
    gas  = [RHOG, 6._R8, 0._R8, 0._R8, 800._R8]
    call cell_surface(cond, gas, RG, GAM, KG, MUG, mat, 6, rho_mat, rp, re, ma, pr, nu, tsd, tb, boil)
    call p2t_surface(cond(5), rho_mat, rp, re, nu, 6._R8, gas, RG, GAM, KG, MUG, mat, tsd, stat)
    call sources(mat, cond, gas, s, ts, 'L16', ok)
    write(*,'(A,I0,A,ES12.4)') '   L16: stat ', stat, ' Ts ', ts
    call check(stat == 0 .and. ieee_is_nan(tsd), "L16 mu-liq 0: p2t_surface's stat 0 and Ts NaN", ok)
    call check(ieee_is_nan(ts), "L16 mu-liq 0: the source's Ts NaN", ok)
  end subroutine leg_no_viscosity


  !> The one-atmosphere anchor stepped one ulp at a time until boilingTemperatureCC gives TB exactly at the gas
  !> pressure of Tg.
  function anchor_for(TB, Tg) result(invTb)
    use ICE_Lib_Evaporation, only: boilingTemperatureCC
    real(R8), intent(in) :: TB, Tg
    real(R8) :: invTb, p, t
    integer  :: k
    p = RHOG*RG*Tg
    invTb = 1._R8/TB + log(p/PATM)/(LV*MV/RU)
    do k = 1, 100000
      t = boilingTemperatureCC(p, LV*MV/RU, invTb)
      if (t == TB) return
      if (t < TB) then
        invTb = ieee_next_after(invTb, -1._R8)
      else
        invTb = ieee_next_after(invTb, 1._R8)
      endif
    enddo
    write(*,'(A)') '   FAIL no one-atmosphere anchor gives T_B exactly'
    error stop 1
  end function anchor_for


  !> L8: water (rho 958, c 4186, d 50 um, k-liq 0.6) at rest, T_g = T_B = 600 K, T = T_B + 30, + 10, + 1, + 0.2 K.
  subroutine leg_flash(ok)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: iinvTboil
    logical, intent(inout) :: ok
    real(R8), parameter :: RHO = 958._R8, C = 4186._R8, KL = 0.6_R8, D0 = 50.e-6_R8
    type(condensed_phase_t) :: mat
    real(R8) :: dts(4), cond(6), gas(5), s(6), ts, m, d, f7, f8, rq, rf, w8, wf
    logical  :: atb
    integer  :: k
    call material(mat, RHO, C, KL, 1.e-3_R8, 'CEM', 0._R8)
    mat%boilSelect = 1
    mat%ep(iinvTboil) = anchor_for(600._R8, 600._R8)
    dts = [30._R8, 10._R8, 1._R8, 0.2_R8]
    m   = RHO*PI*D0**3/6._R8
    gas = [RHOG, 0._R8, 0._R8, 0._R8, 600._R8]
    w8 = 0._R8; wf = 0._R8; atb = .true.
    do k = 1, 4
      cond = [1._R8, 0._R8, 0._R8, 0._R8, 600._R8 + dts(k), 1._R8/m]
      call sources(mat, cond, gas, s, ts, 'L8', ok)
      d  = diameter(cond, RHO)
      atb = atb .and. ts == 600._R8
      f8 = s(1)/cond(6)
      f7 = (s(5) - C*cond(5)*s(1))/(cond(1)*C)
      rq = -10._R8*PI*KL*d*max(dts(k), 0.5_R8)/LV
      rf = C/(LV - C*dts(k))
      w8 = max(w8, abs(f8/rq - 1._R8)); wf = max(wf, abs((f8/f7)/m/rf - 1._R8))
      write(*,'(A,F5.1,A,ES22.14,A,ES20.12,A,ES18.10)') '   L8: dT ', dts(k), ' K: S1 ', s(1), ' F8 ', f8, &
        ' (F8/F7)/m ', (f8/f7)/m
    enddo
    write(*,'(A,ES10.3,A,ES10.3)') '   L8: worst relative: F8 ', w8, ', (F8/F7)/m ', wf
    call check(atb, 'L8 Ts = T_B bit for bit', ok)
    call check(w8 <= 1.e-13_R8, 'L8 F8 = -10 pi k_l d dT_f/Lv (the 0.5 K floor) to 1e-13', ok)
    call check(wf <= 1.e-13_R8, 'L8 (F8/F7)/m = c/(Lv - c dT) to 1e-13', ok)
  end subroutine leg_flash


  !> L6: T_g 300, 600, 1500 K x T 280, 350, T_boil - 0.1, T_boil + 5 K x d 1, 100 um x k-liq 0.01, 0.6, u_g 6 m/s, for
  !> d2-law, CEM, CEM-B, ASM, TC, CEM + LK, CEM with the blowing correction; the clamp, then ZGR; then one ZGR state
  !> with heat-transfer NoHeat in the regime.
  subroutine leg_grid(ok)
    use ICE_Config_Types_m, only: condensed_phase_t, obj_time_scheme
    logical, intent(inout) :: ok
    character(len=8), parameter :: models(7) = [character(len=8) :: 'd2-law', 'CEM', 'CEM-B', 'ASM', 'TC', 'CEM', 'CEM']
    integer,          parameter :: intf(7) = [0, 0, 0, 0, 0, 1, 0], blow(7) = [0, 0, 0, 0, 0, 0, 1]
    real(R8), parameter :: RHO = 1000._R8
    type(condensed_phase_t) :: mat
    real(R8) :: tgs(3), tbs(4), ds(2), kls(2), cond(6), gas(5)
    integer  :: ib, im, ig, it, id, ik, cnt(4)
    tgs = [300._R8, 600._R8, 1500._R8]
    tbs = [280._R8, 350._R8, TBOIL - 0.1_R8, TBOIL + 5._R8]
    ds  = [1.e-6_R8, 100.e-6_R8]
    kls = [0.01_R8, 0.6_R8]
    cnt = 0
    do ib = 0, 1
      do im = 1, 7
        obj_time_scheme%blowSelect = blow(im)
        do ig = 1, 3; do it = 1, 4; do id = 1, 2; do ik = 1, 2
          call material(mat, RHO, 4180._R8, kls(ik), 1.e-3_R8, trim(models(im)), 0._R8)
          mat%intfSelect = intf(im); mat%boilSelect = ib
          cond = [1._R8, 0._R8, 0._R8, 0._R8, tbs(it), 1._R8/(RHO*PI*ds(id)**3/6._R8)]
          gas  = [RHOG, 6._R8, 0._R8, 0._R8, tgs(ig)]
          call grid_point(mat, cond, gas, cnt, ok)
        enddo; enddo; enddo; enddo
      enddo
    enddo
    obj_time_scheme%blowSelect = 0
    ! one ZGR state in the regime without convective heat (B = 0)
    call setup_model('Stokes', 'NoHeat')
    call material(mat, RHO, 4180._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    mat%boilSelect = 1
    cond = [1._R8, 0._R8, 0._R8, 0._R8, 600._R8, 1._R8/(RHO*PI*(50.e-6_R8)**3/6._R8)]
    gas  = [RHOG, 6._R8, 0._R8, 0._R8, 600._R8]
    call grid_point(mat, cond, gas, cnt, ok)
    call setup_model('Stokes', 'Ranz-Marshall')
    write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '   L6: ', cnt(1), ' states, ', cnt(2), ' not finite, ', cnt(3), &
      ' surface balances unsolved, ', cnt(4), ' in the boiling regime'
    call check(cnt(2) == 0, 'L6 every slot and Ts finite', ok)
    call check(cnt(3) == 0, 'L6 every surface balance solved', ok)
  end subroutine leg_grid


  !> One L6 state: cnt = states, not finite, unsolved, in the regime.
  subroutine grid_point(mat, cond, gas, cnt, ok)
    use ICE_Config_Types_m,   only: condensed_phase_t
    use ICE_Mod_Sources,      only: cell_surface
    use ICE_Lib_Drop_Surface, only: p2t_surface
    type(condensed_phase_t), intent(in)    :: mat
    real(R8),                intent(in)    :: cond(6), gas(5)
    integer,                 intent(inout) :: cnt(4)
    logical,                 intent(inout) :: ok
    real(R8) :: s(6), ts, rho_mat, rp, re, ma, pr, nu, tsd, tb
    logical  :: boil, ok1
    integer  :: stat
    ok1 = .true.
    call cell_surface(cond, gas, RG, GAM, KG, MUG, mat, 6, rho_mat, rp, re, ma, pr, nu, tsd, tb, boil)
    call sources_quiet(mat, cond, gas, s, ts, ok1)
    if (.not. ok1) ok = .false.
    cnt(1) = cnt(1) + 1
    if (.not. (all(ieee_is_finite(s)) .and. ieee_is_finite(ts))) cnt(2) = cnt(2) + 1
    if (boil) then
      cnt(4) = cnt(4) + 1
    else
      call p2t_surface(cond(5), rho_mat, rp, re, nu, norm2(gas(2:4) - cond(2:4)), gas, RG, GAM, KG, MUG, mat, tsd, stat)
      if (stat == 0) cnt(3) = cnt(3) + 1
    endif
  end subroutine grid_point


  !> sources without a line per call: false when cell_surface's Ts is not the source's.
  subroutine sources_quiet(mat, cond, gas, s, ts, ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    use ICE_Mod_Sources,    only: compute_source_, cell_surface
    type(condensed_phase_t), intent(in)    :: mat
    real(R8),                intent(in)    :: cond(6), gas(5)
    real(R8),                intent(out)   :: s(6), ts
    logical,                 intent(inout) :: ok
    real(R8) :: tau, rho_mat, rp_, re_, ma_, pr_, nu_, tsc, tb_
    logical  :: boil_
    tau = 0._R8
    call compute_source_(cond, tau, gas, RG, GAM, KG, MUG, s, mat, 6, ts)
    call cell_surface(cond, gas, RG, GAM, KG, MUG, mat, 6, rho_mat, rp_, re_, ma_, pr_, nu_, tsc, tb_, boil_)
    if (.not. (tsc == ts .or. (ieee_is_nan(tsc) .and. ieee_is_nan(ts)))) then
      write(*,'(A)') "   FAIL L6: cell_surface's Ts is not the source's"
      ok = .false.
    endif
  end subroutine sources_quiet


  !> L12, L15: a constant Conductivity column gives the token's sources; the flash takes the column's k_l.
  subroutine leg_columns(ok)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: iinvTboil
    logical, intent(inout) :: ok
    type(condensed_phase_t) :: tok, col
    real(R8) :: cond(6), gas(5), s(6), sc(6), ts, tsc
    integer  :: i
    ! L12: L5's state
    call material(tok, 1000._R8, 4180._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    call material(col, 1000._R8, 4180._R8, 0.3_R8, 1.e-3_R8, 'CEM', 0._R8)
    call constant_table(col, 280, 400, 1000._R8, 4180._R8)
    allocate(col%kl_tab(280:400), source=[(0.6_R8, i = 280, 400)]); col%use_kl = .true.
    cond = [1._R8, 0._R8, 0._R8, 0._R8, 320._R8, 1._R8/(1000._R8*PI*(50.e-6_R8)**3/6._R8)]
    gas  = [RHOG, 6._R8, 0._R8, 0._R8, 800._R8]
    call sources(tok, cond, gas, s, ts, 'L12 (k-liq 0.6)', ok)
    call sources(col, cond, gas, sc, tsc, 'L12 (the column)', ok)
    write(*,'(A,ES22.14,A,ES22.14)') '   L12: S5 token ', s(5), ' column ', sc(5)
    call check(all(s == sc) .and. ts == tsc, 'L12 a constant Conductivity column gives the k-liq sources bit for bit', ok)
    ! L15: L8's state at T_B + 30 K
    call material(tok, 958._R8, 4186._R8, 0.6_R8, 1.e-3_R8, 'CEM', 0._R8)
    tok%boilSelect = 1; tok%ep(iinvTboil) = anchor_for(600._R8, 600._R8)
    call material(col, 958._R8, 4186._R8, 0._R8, 1.e-3_R8, 'CEM', 0._R8)
    col%boilSelect = 1; col%ep(iinvTboil) = tok%ep(iinvTboil)
    call constant_table(col, 280, 1000, 958._R8, 4186._R8)
    allocate(col%kl_tab(280:1000), source=[(0.6_R8, i = 280, 1000)]); col%use_kl = .true.
    cond = [1._R8, 0._R8, 0._R8, 0._R8, 630._R8, 1._R8/(958._R8*PI*(50.e-6_R8)**3/6._R8)]
    gas  = [RHOG, 0._R8, 0._R8, 0._R8, 600._R8]
    call sources(tok, cond, gas, s, ts, 'L15 (k-liq 0.6)', ok)
    call sources(col, cond, gas, sc, tsc, 'L15 (the column)', ok)
    write(*,'(A,ES22.14,A,ES22.14)') '   L15: S1 token ', s(1), ' column ', sc(1)
    call check(all(s == sc) .and. ts == tsc, 'L15 the flash takes the Conductivity column: the k-liq sources bit for bit', &
               ok)
  end subroutine leg_columns


  !> A table of constant density and cp on lo..hi, as Load_Table holds it.
  subroutine constant_table(mat, lo, hi, rho, c)
    use ICE_Config_Types_m, only: condensed_phase_t
    type(condensed_phase_t), intent(inout) :: mat
    integer,                 intent(in)    :: lo, hi
    real(R8),                intent(in)    :: rho, c
    integer :: i
    mat%T_min = lo; mat%T_max = hi; mat%use_table = .true.
    allocate(mat%rho_tab(lo:hi), source=[(rho, i = lo, hi)])
    allocate(mat%cs_tab(lo:hi), source=[(c, i = lo, hi)])
  end subroutine constant_table


  !> L17: L5's state with a Viscosity column mu_l(T) = 1e-3 exp((320 - T)/20) on 280..400 K and mu-liq 0.
  subroutine leg_viscosity(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    logical, intent(inout) :: ok
    real(R8), parameter :: T0 = 320._R8, TG = 800._R8, UG = 6._R8, RHO = 1000._R8, C = 4180._R8, KL = 0.6_R8
    type(condensed_phase_t) :: mat
    real(R8) :: mut(280:400), cond(6), gas(5), s(6), ts, d, re, nu, xs, ys, bm, tor, tfr, md
    integer  :: i
    do i = 280, 400
      mut(i) = 1.e-3_R8*exp((320._R8 - real(i, R8))/20._R8)
    enddo
    call material(mat, RHO, C, KL, 0._R8, 'CEM', 0._R8)
    call constant_table(mat, 280, 400, RHO, C)
    allocate(mat%mul_tab(280:400), source=mut); mat%use_mul = .true.
    cond = [1._R8, 0._R8, 0._R8, 0._R8, T0, 1._R8/(RHO*PI*(50.e-6_R8)**3/6._R8)]
    gas  = [RHOG, UG, 0._R8, 0._R8, TG]
    call sources(mat, cond, gas, s, ts, 'L17', ok, Re=re, Nu=nu)
    d  = diameter(cond, RHO)
    xs = PATM*exp(-LV*MV/RU*(1._R8/T0 - 1._R8/TBOIL))/(RHOG*RG*TG)
    ys = xs*MV/(xs*MV + (1._R8 - xs)*RU/RG)
    bm = ys/(1._R8 - ys)
    tor = root_mu(.true., mut, re, nu, d, bm)
    tfr = root_mu(.false., mut, re, nu, d, bm)
    md  = cem_rate(ts, TG, d, re, MV)
    write(*,'(A,F20.12,A,F20.12,A,F20.12)') '   L17: Ts ', ts, ' oracle ', tor, ' mu_l frozen at T ', tfr
    write(*,'(A,ES10.3,A,ES10.3)') '   L17: relative Ts ', abs(ts/tor - 1._R8), ', S1 ', abs(s(1)/(cond(6)*md) - 1._R8)
    call check(abs(ts/tor - 1._R8) <= 1.e-13_R8, 'L17 Ts the root of the balance with mu_l(Ts) to 1e-13', ok)
    call check(abs(tor - tfr) > 1.e-9_R8, 'L17 guard: the root with mu_l at T differs by more than 1e-9 K', ok)
    call check(abs(s(1)/(cond(6)*md) - 1._R8) <= 1.e-13_R8, 'L17 S1 = n mdot(Ts) to 1e-13', ok)
  end subroutine leg_viscosity


  !> L17's oracle: the root of T - T0 - (Q_g(T) + mdot(T) Lv)/K(mu_l), mu_l at T (atTs) or at T0, by bisection
  !> on [T0, T_g]; L5's drop and gas, the column mut interpolated here.
  function root_mu(atTs, mut, re, nu, d, bm) result(T)
    logical,  intent(in) :: atTs
    real(R8), intent(in) :: mut(280:400), re, nu, d, bm
    real(R8), parameter :: T0 = 320._R8, TG = 800._R8
    real(R8) :: T, lo, hi, glo, gm
    integer  :: it
    lo = T0; hi = TG
    glo = g_mu(lo, atTs, mut, re, nu, d, bm)
    do it = 1, 300
      T = 0.5_R8*(lo + hi)
      if (T <= lo .or. T >= hi) exit
      gm = g_mu(T, atTs, mut, re, nu, d, bm)
      if (gm == 0._R8) exit
      if ((gm > 0._R8) .eqv. (glo > 0._R8)) then; lo = T; glo = gm; else; hi = T; endif
    enddo
  end function root_mu


  !> L17's residual at Tx [K]: U_s, Pe_L and chi (A-S eqs. 30, 32, 39) with mu_l at Tx (atTs) or at T0, the column
  !> linear between its nodes and its end value outside.
  function g_mu(Tx, atTs, mut, re, nu, d, bm) result(g)
    real(R8), intent(in) :: Tx, mut(280:400), re, nu, d, bm
    logical,  intent(in) :: atTs
    real(R8), parameter :: T0 = 320._R8, TG = 800._R8, UG = 6._R8, RHO = 1000._R8, C = 4180._R8, KL = 0.6_R8
    real(R8) :: g, tm, mu, us, chi
    tm = T0; if (atTs) tm = Tx
    mu = value_at(mut, tm)
    us = (12.69_R8/32._R8)*UG*(MUG/mu)*re**(1._R8/3._R8)/(1._R8 + bm)
    chi = 1.86_R8 + 0.86_R8*tanh(2.245_R8*log10(us*d*RHO*C/KL/30._R8))
    g = Tx - T0 - (PI*d*KG*nu*(TG - Tx) + cem_rate(Tx, TG, d, re, MV)*LV)/(10._R8*PI*chi*KL*d)
  end function g_mu


  !> A column at T: linear between the nodes, the end value outside (the lookup's rule, written here).
  pure function value_at(tab, T) result(v)
    real(R8), intent(in) :: tab(280:400), T
    real(R8) :: v, Tc
    integer  :: i
    Tc = min(max(T, 280._R8), 400._R8)
    if (Tc >= 400._R8) then
      v = tab(400)
    else
      i = int(Tc)
      v = tab(i) + (tab(i+1) - tab(i))*(Tc - real(i, R8))
    endif
  end function value_at

end program test_liquid_p2t
