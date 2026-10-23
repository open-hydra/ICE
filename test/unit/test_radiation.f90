!> Unit test of the radiative exchange in the cell source (ICE_Mod_Sources' compute_source_): one cell of MK
!> particles of d = 30 um at rest at T_p = 2400 K in a gas at rest at T_g = 600 K, drag NoDrag, heat-transfer
!> NoHeat and no evaporation, so radiation is the only energy term.
!>   R1 at emissivity 0.8 the energy slot is 0.8 sigma pi d^2 n (T_g^4 - T_p^4), sigma = 5.670374419e-8 and
!>      pi d^2 the whole sphere surface, to 1e-14 relative (the solver's radius is a cube root of the state)
!>   R2 at emissivity 0 the energy slot is 0 and every other slot is R1's
!>   R3 a material whose emissivity is never set radiates nothing
!> Prints every assertion and exits non-zero if any failed.
program test_radiation
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  implicit none
  real(R8), parameter :: RHO_L = 2000._R8, CP_L = 1000._R8
  logical :: ok

  ok = .true.
  call setup_model()
  call legs(ok)
  if (.not. ok) then
    write(*,'(A)') 'test_radiation: FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_radiation: PASS'

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


  !> Family 1 is MK, drag NoDrag and heat-transfer NoHeat, no blowing correction.
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
    call assign_drag('NoDrag', obj_time_scheme%dragSelect)
    call assign_heat('NoHeat', obj_time_scheme%heatSelect)
    obj_time_scheme%blowSelect = 0
  end subroutine setup_model


  !> The same cell with a grey material, a material at emissivity 0 and one left at its default.
  subroutine legs(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    use ICE_Mod_Sources,    only: compute_source_
    logical, intent(inout) :: ok
    real(R8), parameter :: SIGMA = 5.670374419e-8_R8, PI = acos(-1._R8)
    real(R8), parameter :: D = 30.e-6_R8, RHO_P = 1._R8, TP = 2400._R8, TG = 600._R8
    real(R8), parameter :: RG = 287.05_R8, GAM = 1.4_R8, KG = 0.05_R8, MU = 3.e-5_R8
    type(condensed_phase_t) :: grey, off, fresh
    real(R8) :: cond(6), gas(5), s8(6), s0(6), sf(6), t8, t0, tf, n, q, rel

    n    = RHO_P/(RHO_L*PI/6._R8*D**3)
    cond = [RHO_P, 0._R8, 0._R8, 0._R8, TP, n]
    gas  = [101325._R8/(RG*TG), 0._R8, 0._R8, 0._R8, TG]
    call material(grey);  grey%emiss = 0.8_R8
    call material(off);   off%emiss  = 0._R8
    call material(fresh)
    t8 = 0._R8; t0 = 0._R8; tf = 0._R8
    call compute_source_(cond, t8, gas, RG, GAM, KG, MU, s8, grey,  6)
    call compute_source_(cond, t0, gas, RG, GAM, KG, MU, s0, off,   6)
    call compute_source_(cond, tf, gas, RG, GAM, KG, MU, sf, fresh, 6)

    q   = 0.8_R8*SIGMA*PI*D*D*n*(TG*TG*TG*TG - TP*TP*TP*TP)
    rel = abs(s8(5) - q)/abs(q)
    write(*,'(A,2ES24.16)') '   R1: energy slot (ICE, closed form) ', s8(5), q
    write(*,'(A,ES10.3)')   '   R1: relative difference ', rel
    call check(rel <= 1.e-14_R8, 'R1: emissivity 0.8 gives eps sigma pi d^2 n (Tg^4 - Tp^4) to 1e-14', ok)
    call check(s0(5) == 0._R8 .and. all(s0(1:4) == s8(1:4)) .and. s0(6) == s8(6) .and. t0 == t8, &
               'R2: emissivity 0 gives energy 0 and every other slot unchanged', ok)
    write(*,'(A,ES24.16)') '   R3: energy slot, emissivity never set ', sf(5)
    call check(sf(5) == 0._R8, 'R3: a material whose emissivity is never set radiates nothing', ok)
  end subroutine legs


  !> Constant density and cp, no table, no evaporation; every other component at its default.
  subroutine material(mat)
    use ICE_Config_Types_m, only: condensed_phase_t
    type(condensed_phase_t), intent(out) :: mat
    mat%evapSelect = 0; mat%intfSelect = 0
    mat%rho_al = RHO_L; mat%cs_al = CP_L
  end subroutine material

end program test_radiation
