!>@brief The surface temperature of a family's drops, which the gas-side models take: the boiling temperature at
!> the gas pressure when the mean temperature reaches it (boiling = ZGR), else the root of the surface balance of
!> the parabolic temperature profile inside the drop (liquid-conduction = P2T), else the mean temperature (ITC),
!> also for a material that does not evaporate and in an empty cell. IGLOO's models, per drop.
module ICE_Lib_Drop_Surface
  use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private
  public :: drop_surface, p2t_surface, boiling_temperature

contains

  !> Ts of the cell's drops; NaN when the P2T surface balance has no root. boiling: the mean temperature at or above
  !  the boiling temperature TB at the gas pressure (ZGR; TB huge when not evaluated).
  subroutine drop_surface(cond_prim, gas_prim, gas_R, gas_gam, gas_k, gas_mu, mat, nb, rho_mat, Rp, Re, Nu, Ts, &
                          TB, boiling)
    use ICE_Global_m,       only: rho_empty
    use ICE_Config_Types_m, only: condensed_phase_t
    implicit none
    real(R8), dimension(:),  intent(in)  :: cond_prim, gas_prim
    real(R8),                intent(in)  :: gas_R, gas_gam, gas_k, gas_mu, rho_mat, Rp, Re, Nu
    type(condensed_phase_t), intent(in)  :: mat
    integer(I4),             intent(in)  :: nb
    real(R8),                intent(out) :: Ts, TB
    logical,                 intent(out) :: boiling
    integer :: stat

    Ts = cond_prim(nb-1)
    TB = huge(1._R8)
    boiling = .false.
    if (mat%evapSelect == 0 .or. cond_prim(1) <= rho_empty) return
    if (mat%boilSelect == 1) then
      TB = boiling_temperature(mat, gas_prim(1)*gas_R*gas_prim(size(gas_prim)))
      boiling = cond_prim(nb-1) >= TB
      if (boiling) then
        Ts = TB
        return
      endif
    endif
    if (mat%liqSelect == 1) &
      call p2t_surface(cond_prim(nb-1), rho_mat, Rp, Re, Nu, norm2(gas_prim(2:4)-cond_prim(2:4)), &
                       gas_prim, gas_R, gas_gam, gas_k, gas_mu, mat, Ts, stat)
  end subroutine drop_surface


  !> P2T surface temperature (parabolic profile inside the drop): the root of G(Ts) = Ts - tbar - Qnet(Ts)/K,
  !  K = 10 pi keff d, keff = chi(Pe_L) kl (Abramzon-Sirignano 1989 eqs. 25, 30, 32, 39): kl, cl and B_M at tbar;
  !  Qnet [W] the gas-side, radiative and latent heats at Ts; Brent on an exact bracket. stat: 1 the root
  !  bracketed by tbar and the first estimate, 2 after widening, 0 no bracket, no convergence or mu_l not
  !  positive (ts NaN).
  subroutine p2t_surface(tbar, rho_mat, Rp, Re, Nu, slip, gas_prim, gas_R, gas_gam, gas_k, gas_mu, mat, ts, stat)
    use ICE_Parameters_m,    only: pi, sigma_SB
    use ICE_Config_Types_m,  only: condensed_phase_t, obj_time_scheme
    use ICE_Lib_Evaporation, only: evaporation, blowingFactor, chiETC, nep, iLv, icpv
    use ICE_Lib_Properties,  only: mat_psat, mat_lv, mat_cpv, mat_emiss, mat_e, mat_cp_const
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    real(R8),                intent(in)  :: tbar, rho_mat, Rp, Re, Nu, slip
    real(R8), dimension(:),  intent(in)  :: gas_prim
    real(R8),                intent(in)  :: gas_R, gas_gam, gas_k, gas_mu
    type(condensed_phase_t), intent(in)  :: mat
    real(R8),                intent(out) :: ts
    integer,                 intent(out) :: stat
    integer,  parameter :: maxEval = 100
    real(R8), parameter :: tol4 = 4._R8*epsilon(1._R8)
    real(R8) :: d, kl, mul, cl, BM, K, q0, tLo, t1, g0, g1, Tg, Tmin, Tmax
    real(R8) :: a, b, c, fa, fb, fc, dd, e, s, p, q, r, tol1, xm
    integer  :: nEval, sw, ng

    stat = 0
    ts = ieee_value(ts, ieee_quiet_nan)
    ng = size(gas_prim)
    Tg = gas_prim(ng)
    d  = 2._R8*Rp
    kl  = mat%kLiq
    mul = mat%muLiq
    if (mat%cs_varies) then
      cl = mat_e(mat, tbar + 0.5_R8) - mat_e(mat, tbar - 0.5_R8)
    else
      cl = mat_cp_const(mat)
    endif
    Tmin = 0._R8; Tmax = 0._R8
    if (mat%use_table) then
      Tmin = real(mat%T_min, R8); Tmax = real(mat%T_max, R8)
    endif
    !> no liquid viscosity reached the solve (the inputs refuse mu_l <= 0): no root
    if (.not. mul > 0._R8) return
    q0 = qnet(tbar, BM)
    if (.not. ieee_is_finite(q0)) return
    if (q0 == 0._R8) then
      ts = tbar; stat = 1; return
    endif
    K = kAt(mul)
    !> bracket: tbar and the first estimate T1 = tbar + Qnet(tbar)/K, widened to Tg (heating) or the low end (cooling)
    tLo = max(1._R8, Tmin - 1._R8)
    g0  = -q0/K
    t1  = min(max(tbar + q0/K, tLo), max(Tg, Tmax))
    g1  = gfun(t1)
    nEval = 2
    if (ieee_is_finite(g1) .and. g0*g1 <= 0._R8) then
      a = tbar; fa = g0; b = t1; fb = g1; sw = 1
    else
      sw = 2
      if (ieee_is_finite(g1)) then; a = t1; fa = g1; else; a = tbar; fa = g0; endif
      if (q0 > 0._R8) then
        if (.not. Tg > tbar) return
        b = Tg
      else
        b = tLo
      endif
      if (b == a) return
      fb = gfun(b); nEval = nEval + 1
      if (.not. ieee_is_finite(fb)) return
      if (fa*fb > 0._R8) return
    endif
    !> Brent (zeroin): inverse quadratic / secant steps, bisection safeguard
    c = a; fc = fa; dd = b - a; e = dd
    do
      if (abs(fc) < abs(fb)) then
        a = b; b = c; c = a; fa = fb; fb = fc; fc = fa
      endif
      tol1 = 0.5_R8*tol4*abs(b)
      xm   = 0.5_R8*(c - b)
      if (abs(xm) <= tol1 .or. abs(fb) <= tol4*tbar) then
        ts = b; stat = sw; return
      endif
      if (nEval >= maxEval) return
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb/fa
        if (a == c) then
          p = 2._R8*xm*s; q = 1._R8 - s
        else
          q = fa/fc; r = fb/fc
          p = s*(2._R8*xm*q*(q - r) - (b - a)*(r - 1._R8))
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
      a = b; fa = fb
      if (abs(dd) > tol1) then; b = b + dd; else; b = b + sign(tol1, xm); endif
      fb = gfun(b); nEval = nEval + 1
      if (.not. ieee_is_finite(fb)) return
      if ((fb > 0._R8 .and. fc > 0._R8) .or. (fb < 0._R8 .and. fc < 0._R8)) then
        c = a; fc = fa; dd = b - a; e = dd
      endif
    enddo

  contains

    !> The surface balance's residual [K].
    function gfun(T) result(g)
      real(R8), intent(in) :: T
      real(R8) :: g
      g = T - tbar - qnet(T)/K
    end function gfun

    !> 10 pi keff d at the liquid viscosity mu: the internal circulation's surface velocity (A-S eqs. 30, 32),
    !  Pe_L = Re_L Pr_L, chi(Pe_L) (eq. 39).
    function kAt(mu) result(Kmu)
      real(R8), intent(in) :: mu
      real(R8) :: Kmu, Us, PeL
      Us = 0._R8
      if (Re > 0._R8) Us = (12.69_R8/32._R8)*slip*(gas_mu/mu)*Re**(1._R8/3._R8)/(1._R8 + BM)
      PeL = Us*d*rho_mat*cl/kl
      Kmu = 10._R8*pi*chiETC(PeL)*kl*d
    end function kAt

    !> The heat the liquid absorbs at the surface temperature T [W]; BMo, the Spalding number there.
    function qnet(T, BMo) result(qn)
      real(R8), intent(in) :: T
      real(R8), intent(out), optional :: BMo
      real(R8) :: qn, Qg, md, Qe, epT(nep)
      logical  :: ovr
      Qg = 2._R8*Nu*gas_k*pi*Rp*(Tg - T)
      epT = mat%ep
      if (mat%use_lv)  epT(iLv)  = mat_lv(mat, T)
      if (mat%use_cpv) epT(icpv) = mat_cpv(mat, T)
      if (mat%use_psat) then
        call evaporation(gas_prim(1), Tg, gas_gam, gas_R, gas_mu, gas_k, T, d, Re, mat%evapSelect, mat%intfSelect, &
                         epT, md, Qe, ovr, psatExt=mat_psat(mat, T), BMout=BMo)
      else
        call evaporation(gas_prim(1), Tg, gas_gam, gas_R, gas_mu, gas_k, T, d, Re, mat%evapSelect, mat%intfSelect, &
                         epT, md, Qe, ovr, BMout=BMo)
      endif
      if (ovr) then
        Qg = Qe
      elseif (obj_time_scheme%blowSelect == 1) then
        Qg = Qg*blowingFactor(gas_gam, gas_R, gas_mu, gas_k, rho_mat, d, 4._R8/3._R8*pi*Rp**3*rho_mat, md)
      endif
      Qg = Qg + mat_emiss(mat, T)*sigma_SB*4._R8*pi*Rp*Rp*(Tg**4._I4 - T**4._I4)
      qn = Qg + md*epT(iLv)
    end function qnet

  end subroutine p2t_surface


  !> The boiling temperature at the pressure p: the inverse of the material's Psat column, else of Clausius-Clapeyron.
  pure function boiling_temperature(mat, p) result(TB)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: boilingTemperatureCC, iLvMvOverRu, iinvTboil
    use ICE_Lib_Properties,  only: mat_T_from_psat
    implicit none
    type(condensed_phase_t), intent(in) :: mat
    real(R8),                intent(in) :: p
    real(R8) :: TB
    if (mat%use_psat) then
      TB = mat_T_from_psat(mat, p)
    else
      TB = boilingTemperatureCC(p, mat%ep(iLvMvOverRu), mat%ep(iinvTboil))
    endif
  end function boiling_temperature

end module ICE_Lib_Drop_Surface
