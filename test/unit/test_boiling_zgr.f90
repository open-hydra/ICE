!> Unit test of boiling = ZGR, the Zuo-Gomes-Rutland boiling branch above the boiling temperature at the gas
!> pressure, through ICE_Lib_Evaporation's zgrBoiling, alphaBoil and boilingTemperatureCC, ICE_Lib_Drop_Surface's
!> boiling_temperature, ICE_Lib_Properties' mat_T_from_psat and the cell source (compute_source_, cell_surface),
!> against values and identities typed here. IGLOO's water (Lv 2.26e6, Mv 18, c 4186, rho 958, d 50 um) in the box
!> gas (rho 1.2, R 287, gamma 1.4, k 0.026, c_p 1004.5), Re 0, Nu 2, CEM below the boiling temperature T_B.
!>   Z1 the heat-fed rate G of G (1 + Gf/G) = B ln(1 + A (1 + Gf/G)), B with the 1/3-rule film conductivity
!>   Z2 the limits A < 0 (G = 0, the film heat negative) and A = 0 (no film heat)
!>   Z3 alpha_s: its three laws and the 0.5 K floor
!>   Z4 the adiabatic flash at A = 0 through the source: Ts = T_B, (F8/F7)/m = c/(Lv - c dT) (the vapour leaves at T_B)
!>   Z5 T_B(p): Clausius-Clapeyron's inverse, and the inverse of a Psat column inside its 1-K interval
!>   Z6 a Psat column 1.2 x Clausius-Clapeyron: T_B is the column's root, not Clausius-Clapeyron's
!>   Z7 the regime switch: below T_B the clamp's sources bit for bit, at T_B the boiling branch's rate
!> Prints every assertion and exits non-zero if any failed.
program test_boiling_zgr
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_next_after
  implicit none
  real(R8), parameter :: PI = 3.14159265358979323846_R8, RU = 8314.46_R8, PATM = 101325._R8
  real(R8), parameter :: RHOG = 1.2_R8, MUG = 1.8e-5_R8, KG = 0.026_R8, GAM = 1.4_R8, RG = 287._R8
  real(R8), parameter :: LV = 2.26e6_R8, MV = 18._R8, CPV = 1900._R8, RHOL = 958._R8, CL = 4186._R8, D0 = 50.e-6_R8
  real(R8), parameter :: CPGA = 1004.5_R8
  logical :: ok

  ok = .true.
  call setup_model()
  call leg_heat_fed(ok)
  call leg_limits(ok)
  call leg_alpha(ok)
  call leg_flash(ok)
  call leg_tb(ok)
  call leg_tb_column(ok)
  call leg_switch(ok)
  if (.not. ok) then
    write(*,'(A)') 'test_boiling_zgr: FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_boiling_zgr: PASS'

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


  !> Family 1 is MK, drag Stokes, heat-transfer Ranz-Marshall (Nu 2 at Re 0), no blowing correction.
  subroutine setup_model()
    use ICE_Config_Types_m, only: obj_time_scheme
    use ICE_Global_m,       only: solid_of
    use ICE_Lib_Model,      only: assign_all
    use ICE_Lib_Drag,       only: assign_drag
    use ICE_Lib_Heat,       only: assign_heat
    allocate(obj_time_scheme%model(1))
    obj_time_scheme%model(1) = 'MK'
    allocate(solid_of(1))
    solid_of = .false.
    call assign_all(1)
    call assign_drag('Stokes', obj_time_scheme%dragSelect)
    call assign_heat('Ranz-Marshall', obj_time_scheme%heatSelect)
    obj_time_scheme%blowSelect = 0
  end subroutine setup_model


  !> alpha_s [W/m2/K] at the superheat dT, typed from OpenFOAM's LiquidEvaporationBoil (the 0.5 K floor, three laws).
  pure function alpha_law(dT) result(a)
    real(R8), intent(in) :: dT
    real(R8) :: a, t
    t = max(dT, 0.5_R8)
    if (t < 5._R8) then
      a = 760._R8*t**0.26_R8
    elseif (t < 25._R8) then
      a = 27._R8*t**2.33_R8
    else
      a = 13800._R8*t**0.39_R8
    endif
  end function alpha_law


  !> Clausius-Clapeyron, typed: Patm exp(-(Lv Mv/Ru)(1/T - invT0)).
  pure function psat_cc(T, invT0) result(ps)
    real(R8), intent(in) :: T, invT0
    real(R8) :: ps
    ps = PATM*exp(-(LV*MV/RU)*(1._R8/T - invT0))
  end function psat_cc


  !> Z1: d 50 um, Nu 2, c_p 1004.5, Tg 600 K, T_B 450 K (A = 0.066670), dT 20, 5, 0.5 K.
  subroutine leg_heat_fed(ok)
    use ICE_Lib_Evaporation, only: zgrBoiling
    logical, intent(inout) :: ok
    real(R8) :: dts(3), gref(3), a, tf, kf, b, md, qd, gf, g, worst, worst2
    integer  :: k
    dts  = [20._R8, 5._R8, 0.5_R8]
    gref = [4.0108994750e-10_R8, 4.6131345055e-10_R8, 4.6190847194e-10_R8]
    a  = CPGA*(600._R8 - 450._R8)/LV
    tf = 450._R8 + (600._R8 - 450._R8)/3._R8
    kf = KG*(tf/600._R8)**0.7_R8
    b  = PI*kf*D0*2._R8/CPGA
    worst = 0._R8; worst2 = 0._R8
    do k = 1, 3
      call zgrBoiling(450._R8 + dts(k), 450._R8, D0, KG, CPGA, 2._R8, 600._R8, LV, md, qd)
      gf = alpha_law(dts(k))*max(dts(k), 0.5_R8)*PI*D0*D0/LV
      g  = -md - gf
      worst  = max(worst, abs(g/gref(k) - 1._R8))
      worst2 = max(worst2, abs(g*(1._R8 + gf/g) - b*log(1._R8 + a*(1._R8 + gf/g)))/b)
      write(*,'(A,F5.1,A,ES18.10,A,ES18.10)') '   Z1: dT ', dts(k), ' K: G ', g, ' Gf ', gf
    enddo
    write(*,'(A,ES10.3,A,ES10.3)') '   Z1: worst G/G_ref - 1 ', worst, ', worst residual/B ', worst2
    call check(worst <= 1.e-9_R8, 'Z1 G at dT 20, 5, 0.5 K against the ten-digit references to 1e-9', ok)
    call check(worst2 <= 1.e-13_R8, 'Z1 |G (1 + Gf/G) - B ln(1 + A (1 + Gf/G))|/B to 1e-13', ok)
  end subroutine leg_heat_fed


  !> Z2: Tg 440 K (A < 0) and Tg = T_B 450 K (A = 0), dT 20 and 0.5 K.
  subroutine leg_limits(ok)
    use ICE_Lib_Evaporation, only: zgrBoiling
    logical, intent(inout) :: ok
    real(R8) :: tf, kf, b, dt, md, qd, gf, qs, worst, worst2
    logical  :: sgn
    integer  :: k
    worst = 0._R8; worst2 = 0._R8; sgn = .true.
    tf = 450._R8 + (440._R8 - 450._R8)/3._R8
    kf = KG*(tf/440._R8)**0.7_R8
    b  = PI*kf*D0*2._R8/CPGA
    do k = 1, 2
      dt = merge(20._R8, 0.5_R8, k == 1)
      call zgrBoiling(450._R8 + dt, 450._R8, D0, KG, CPGA, 2._R8, 440._R8, LV, md, qd)
      gf = alpha_law(dt)*max(dt, 0.5_R8)*PI*D0*D0/LV
      qs = gf*CPGA*(440._R8 - 450._R8)/(2._R8*sinh(0.5_R8*gf/b)*exp(0.5_R8*gf/b))
      worst  = max(worst, abs(-md/gf - 1._R8))
      worst2 = max(worst2, abs(qd/qs - 1._R8))
      sgn = sgn .and. qd < 0._R8
      write(*,'(A,F5.1,A,ES16.8,A,ES16.8,A)') '   Z2: A < 0, dT ', dt, ' K: mdot ', md, ' Q_s ', qd, ' W'
      call zgrBoiling(450._R8 + dt, 450._R8, D0, KG, CPGA, 2._R8, 450._R8, LV, md, qd)
      sgn = sgn .and. qd == 0._R8
      worst = max(worst, abs(-md/gf - 1._R8))
    enddo
    call check(worst <= 1.e-14_R8, 'Z2 mdot = -Gf at A <= 0 (G = 0) to 1e-14', ok)
    call check(worst2 <= 1.e-14_R8, 'Z2 Q_s = Gf c_p (Tg - T_B)/(exp(Gf/B) - 1) at A < 0 to 1e-14', ok)
    call check(sgn, 'Z2 Q_s < 0 at A < 0, Q_s = 0 exactly at A = 0', ok)
  end subroutine leg_limits


  !> Z3: alpha_s at dT 0.2 (the floor), 1, 10, 30 K.
  subroutine leg_alpha(ok)
    use ICE_Lib_Evaporation, only: alphaBoil
    logical, intent(inout) :: ok
    real(R8) :: dtc(4), got(4), aref(4)
    integer  :: k
    dtc  = [0.2_R8, 1._R8, 10._R8, 30._R8]
    aref = [760._R8*0.5_R8**0.26_R8, 760._R8, 27._R8*10._R8**2.33_R8, 13800._R8*30._R8**0.39_R8]
    got  = [(alphaBoil(dtc(k)), k = 1, 4)]
    write(*,'(A,4F22.12)') '   Z3: alpha_s ', got
    call check(all(abs(got/aref - 1._R8) <= 1.e-14_R8), 'Z3 alpha_s at dT 0.2, 1, 10, 30 K to 1e-14', ok)
  end subroutine leg_alpha


  !> A water material of the box, CEM below T_B, boiling = ZGR, its one-atmosphere anchor invTb.
  subroutine material(mat, invTb, boil)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: assign_evaporation, iMv, iLv, icpv, iLe, iYinf, iLvMvOverRu, iinvTboil, ialphaE
    type(condensed_phase_t), intent(out) :: mat
    real(R8),                intent(in)  :: invTb
    integer,                 intent(in)  :: boil
    call assign_evaporation('CEM', mat%evapSelect)
    mat%intfSelect = 0; mat%liqSelect = 0; mat%boilSelect = boil
    mat%rho_al = RHOL; mat%cs_al = CL; mat%emiss = 0._R8
    mat%lv_al = LV; mat%Mv = MV; mat%cpv = CPV; mat%Le = 1._R8; mat%Yinf = 0._R8; mat%Tboil = 1._R8/invTb
    mat%alphaE = 1._R8
    mat%ep = 0._R8
    mat%ep(iMv) = MV; mat%ep(iLv) = LV; mat%ep(icpv) = CPV; mat%ep(iLe) = 1._R8; mat%ep(iYinf) = 0._R8
    mat%ep(iLvMvOverRu) = LV*MV/RU; mat%ep(iinvTboil) = invTb; mat%ep(ialphaE) = 1._R8
  end subroutine material


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


  !> Z4: the adiabatic flash at A = 0 (Tg = T_B = 600 K), rho_p 1, n = 1/m, T = T_B + 30, + 10, + 1 K.
  subroutine leg_flash(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    use ICE_Mod_Sources,    only: compute_source_
    logical, intent(inout) :: ok
    type(condensed_phase_t) :: mat
    real(R8) :: dtd(3), cond(6), gas(5), s(6), tau, ts, f7, f8, m, ref, worst
    logical  :: atb
    integer  :: k
    call material(mat, anchor_for(600._R8, 600._R8), 1)
    dtd = [30._R8, 10._R8, 1._R8]
    m = RHOL*PI*D0**3/6._R8
    gas = [RHOG, 0._R8, 0._R8, 0._R8, 600._R8]
    worst = 0._R8; atb = .true.
    do k = 1, 3
      cond = [1._R8, 0._R8, 0._R8, 0._R8, 600._R8 + dtd(k), 1._R8/m]
      tau = 0._R8
      call compute_source_(cond, tau, gas, RG, GAM, KG, MUG, s, mat, 6, ts)
      atb = atb .and. ts == 600._R8
      f8  = s(1)/cond(6)
      f7  = (s(5) - CL*cond(5)*s(1))/(cond(1)*CL)
      ref = CL/(LV - CL*dtd(k))
      worst = max(worst, abs((f8/f7)/m/ref - 1._R8))
      write(*,'(A,F5.1,A,ES22.14,A,ES22.14,A,ES18.10,A,ES18.10)') '   Z4: dT ', dtd(k), ' K: S1 ', s(1), ' S5 ', s(5), &
        ' (F8/F7)/m ', (f8/f7)/m, ' c/(Lv - c dT) ', ref
    enddo
    write(*,'(A,ES10.3)') '   Z4: worst relative ', worst
    call check(atb, 'Z4 Ts = T_B bit for bit', ok)
    call check(worst <= 1.e-13_R8, 'Z4 (F8/F7)/m = c/(Lv - c dT) at A = 0 to 1e-13', ok)
  end subroutine leg_flash


  !> Z5: no table, p 0.2, 1, 5 atm; then a Psat column of Clausius-Clapeyron on 1-K nodes 280..500 K.
  subroutine leg_tb(ok)
    use ICE_Config_Types_m,   only: condensed_phase_t
    use ICE_Lib_Drop_Surface, only: boiling_temperature
    use ICE_Lib_Properties,   only: mat_T_from_psat, mat_psat
    logical, intent(inout) :: ok
    type(condensed_phase_t) :: mat
    real(R8) :: pr(3), tb, x, worst, worst2
    logical  :: inside
    integer  :: i, k
    call material(mat, 1._R8/373.15_R8, 1)
    pr = [0.2_R8, 1._R8, 5._R8]*PATM
    worst = 0._R8
    do k = 1, 3
      tb = boiling_temperature(mat, pr(k))
      worst = max(worst, abs(psat_cc(tb, 1._R8/373.15_R8)/pr(k) - 1._R8))
      write(*,'(A,F6.2,A,F18.10,A)') '   Z5: p ', pr(k)/PATM, ' atm: T_B ', tb, ' K (Clausius-Clapeyron)'
    enddo
    call check(worst <= 1.e-14_R8, 'Z5 psat_CC(T_B(p)) = p at 0.2, 1, 5 atm to 1e-14', ok)
    mat%T_min = 280; mat%T_max = 500; mat%use_table = .true.; mat%use_psat = .true.
    allocate(mat%psat_tab(280:500), source=[(psat_cc(real(i, R8), 1._R8/373.15_R8), i = 280, 500)])
    worst2 = 0._R8; inside = .true.
    do k = 1, 3
      x = mat_T_from_psat(mat, pr(k))
      worst2 = max(worst2, abs(mat_psat(mat, x)/pr(k) - 1._R8))
      i = int(x)
      inside = inside .and. mat%psat_tab(i) <= pr(k) .and. pr(k) <= mat%psat_tab(min(i + 1, 500))
      inside = inside .and. boiling_temperature(mat, pr(k)) == x
      write(*,'(A,F6.2,A,F18.10,A)') '   Z5: p ', pr(k)/PATM, ' atm: T_B ', x, ' K (the column)'
    enddo
    write(*,'(A,ES10.3,A,ES10.3)') '   Z5: worst relative: Clausius-Clapeyron ', worst, ', the column ', worst2
    call check(worst2 <= 1.e-12_R8, 'Z5 psat(T_from_psat(p)) = p on the column to 1e-12', ok)
    call check(inside, 'Z5 the root inside its 1-K interval, and boiling_temperature takes the column', ok)
  end subroutine leg_tb


  !> Z6: a Psat column 1.2 x Clausius-Clapeyron on 280..500 K at p = 1 atm.
  subroutine leg_tb_column(ok)
    use ICE_Config_Types_m,   only: condensed_phase_t
    use ICE_Lib_Drop_Surface, only: boiling_temperature
    logical, intent(inout) :: ok
    type(condensed_phase_t) :: mat
    real(R8) :: tab(280:500), tb, lo, hi, x, v
    integer  :: i, it
    call material(mat, 1._R8/373.15_R8, 1)
    do i = 280, 500
      tab(i) = 1.2_R8*psat_cc(real(i, R8), 1._R8/373.15_R8)
    enddo
    mat%T_min = 280; mat%T_max = 500; mat%use_table = .true.; mat%use_psat = .true.
    allocate(mat%psat_tab(280:500), source=tab)
    tb = boiling_temperature(mat, PATM)
    ! the interpolant's root by bisection, typed here
    lo = 280._R8; hi = 500._R8
    do it = 1, 200
      x = 0.5_R8*(lo + hi)
      i = min(int(x), 499)
      v = tab(i) + (tab(i+1) - tab(i))*(x - real(i, R8))
      if (v < PATM) then; lo = x; else; hi = x; endif
    enddo
    write(*,'(A,F18.10,A,F18.10,A)') '   Z6: T_B ', tb, ' K, the column''s root ', x, ' K'
    call check(abs(tb/x - 1._R8) <= 1.e-12_R8, 'Z6 T_B is the Psat column''s root to 1e-12 (not 373.15 K)', ok)
  end subroutine leg_tb_column


  !> Z7: Tg 600 K, T_B 450 K; T = T_B (1 - 1e-12) and T = T_B.
  subroutine leg_switch(ok)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Mod_Sources,     only: compute_source_, cell_surface
    use ICE_Lib_Evaporation, only: zgrBoiling
    logical, intent(inout) :: ok
    type(condensed_phase_t) :: zgr, clamp
    real(R8) :: cond(6), gas(5), s(6), sc(6), tau, ts, m, rho_mat, rp, re, ma, pr, nu, tsc, tb, md, qd
    logical  :: boil
    m = RHOL*PI*D0**3/6._R8
    call material(zgr, anchor_for(450._R8, 600._R8), 1)
    call material(clamp, anchor_for(450._R8, 600._R8), 0)
    gas = [RHOG, 0._R8, 0._R8, 0._R8, 600._R8]
    cond = [1._R8, 0._R8, 0._R8, 0._R8, 450._R8*(1._R8 - 1.e-12_R8), 1._R8/m]
    tau = 0._R8
    call compute_source_(cond, tau, gas, RG, GAM, KG, MUG, s, zgr, 6, ts)
    call compute_source_(cond, tau, gas, RG, GAM, KG, MUG, sc, clamp, 6)
    write(*,'(A,ES22.14,A,ES22.14)') '   Z7: T = T_B (1 - 1e-12): S1 ZGR ', s(1), ' clamp ', sc(1)
    call check(all(s == sc), 'Z7 below T_B: the clamp''s sources bit for bit', ok)
    cond(5) = 450._R8
    call compute_source_(cond, tau, gas, RG, GAM, KG, MUG, s, zgr, 6, ts)
    call cell_surface(cond, gas, RG, GAM, KG, MUG, zgr, 6, rho_mat, rp, re, ma, pr, nu, tsc, tb, boil)
    call zgrBoiling(450._R8, tb, 2._R8*rp, KG, GAM*RG/(GAM - 1._R8), nu, 600._R8, LV, md, qd)
    write(*,'(A,ES22.14,A,ES22.14,A,F18.12)') '   Z7: T = T_B: S1 ', s(1), ' n mdot ', md*cond(6), ' T_B ', tb
    call check(boil .and. tb == 450._R8 .and. ts == 450._R8, 'Z7 at T_B: the regime, Ts = T_B', ok)
    call check(s(1) == -(-md*cond(6)), 'Z7 at T_B: S1 = n mdot of zgrBoiling bit for bit', ok)
  end subroutine leg_switch

end program test_boiling_zgr
