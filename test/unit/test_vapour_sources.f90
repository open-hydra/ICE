!> Unit test of the vapour columns in the cell source (ICE_Mod_Sources' compute_source_): one cell of
!> case F's gas (air at 800 K and one atmosphere, at rest) holding water droplets of 1e-4 m at rest at
!> T_p = 337.5 K, closure MK, drag and heat-transfer Stokes, a Psat column on 280..400 K. A material
!> with Lv (2.4e6 - 5000 (T - 300) J/kg) and Cpv (1800 + 2 (T - 300) J/(kg K)) tabulated and both
!> scalars poisoned (1e30) must give the source of the same material with the tables off and the
!> scalars set to the tables at T_p: the mass, momentum and number slots and the relaxation time
!> exactly, the energy slot to 1e-15. T_p lies half a kelvin past a node and both tables have integer
!> slopes, so the table value there is exact however the compiler contracts the interpolation.
!>   VS1 d2-law: the rate (B_T reads Lv) and the latent sink
!>   VS2 CEM: the latent sink
!>   VS3 ASM: the gas-side heat ASM computes with Cpv, and the latent sink
!>   VS4 TC: the rate (the vapour Lewis number reads Cpv), the gas-side heat and the latent sink
!> Every internal procedure takes what it reads through its arguments. Prints every assertion and
!> exits non-zero if any failed.
program test_vapour_sources
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  implicit none
  logical :: ok

  ok = .true.
  call setup_model()
  call leg('VS1 d2-law', 'd2-law', ok)
  call leg('VS2 CEM',    'CEM',    ok)
  call leg('VS3 ASM',    'ASM',    ok)
  call leg('VS4 TC',     'TC',     ok)
  if (.not. ok) then
    write(*,'(A)') 'test_vapour_sources: FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_vapour_sources: PASS'

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


  !> Family 1 is MK, drag and heat-transfer Stokes, no blowing correction.
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
    call assign_heat('Stokes', obj_time_scheme%heatSelect)
    obj_time_scheme%blowSelect = 0
  end subroutine setup_model


  !> One evaporation model: the tabulated material against the reference, in one cell.
  subroutine leg(name, model, ok)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: assign_evaporation
    use ICE_Mod_Sources,     only: compute_source_
    character(len=*), intent(in)    :: name, model
    logical,          intent(inout) :: ok
    real(R8), parameter :: Tp = 337.5_R8, TG = 800._R8, RG = 287.05_R8, GAM = 1.4_R8, KG = 0.026_R8, MU = 1.8e-5_R8
    real(R8), parameter :: PI = acos(-1._R8), RHO_L = 1000._R8, D0 = 1.e-4_R8, RHO_P = 1._R8
    type(condensed_phase_t) :: tab, ref
    real(R8) :: lvt(280:400), cvt(280:400), pst(280:400), cond(6), gas(5), sa(6), sb(6), ta, tb
    integer  :: i, sel

    do i = 280, 400
      lvt(i) = 2.4e6_R8 - 5000._R8*real(i - 300, R8)
      cvt(i) = 1800._R8 + 2._R8*real(i - 300, R8)
      pst(i) = 101325._R8*exp(-4896.8_R8*(1._R8/real(i, R8) - 1._R8/373.15_R8))
    enddo
    call assign_evaporation(model, sel)

    ! the tables on, the scalars poisoned
    call material(tab, sel, pst)
    tab%use_lv = .true.; tab%use_cpv = .true.
    allocate(tab%lv_tab(280:400), source=lvt)
    allocate(tab%cpv_tab(280:400), source=cvt)
    tab%lv_al = 1.e30_R8; tab%cpv = 1.e30_R8
    call pack(tab)
    ! the tables off, the scalars at T_p
    call material(ref, sel, pst)
    ref%lv_al = value_at(lvt, Tp); ref%cpv = value_at(cvt, Tp)
    call pack(ref)

    cond = [RHO_P, 0._R8, 0._R8, 0._R8, Tp, RHO_P/(RHO_L*4._R8/3._R8*PI*(0.5_R8*D0)**3)]
    gas  = [101325._R8/(RG*TG), 0._R8, 0._R8, 0._R8, TG]
    ta = 0._R8; tb = 0._R8
    call compute_source_(cond, ta, gas, RG, GAM, KG, MU, sa, tab, 6)
    call compute_source_(cond, tb, gas, RG, GAM, KG, MU, sb, ref, 6)

    write(*,'(A,2ES24.16)') '   '//name//': mass slot   (table, reference) ', sa(1), sb(1)
    write(*,'(A,2ES24.16)') '   '//name//': energy slot (table, reference) ', sa(5), sb(5)
    call check(sb(1) < 0._R8, name//': the droplet evaporates', ok)
    call check(sa(1) == sb(1), name//': mass slot exact', ok)
    call check(all(sa(2:4) == sb(2:4)) .and. sa(6) == sb(6) .and. ta == tb, &
               name//': momentum and number slots and the relaxation time exact', ok)
    call check(abs(sa(5) - sb(5)) <= 1.e-15_R8*abs(sb(5)), name//': energy slot to 1e-15', ok)
  end subroutine leg


  !> Case F's water, at rest in no table but the Psat column: constant density and cp from the INI values.
  subroutine material(mat, sel, pst)
    use ICE_Config_Types_m, only: condensed_phase_t
    type(condensed_phase_t), intent(out) :: mat
    integer,                 intent(in)  :: sel
    real(R8),                intent(in)  :: pst(280:400)
    mat%evapSelect = sel; mat%intfSelect = 0
    mat%rho_al = 1000._R8; mat%cs_al = 4180._R8; mat%emiss = 0._R8
    mat%Mv = 18.015_R8; mat%Tboil = 373.15_R8; mat%Le = 1._R8; mat%Yinf = 0._R8; mat%alphaE = 1._R8
    mat%T_min = 280; mat%T_max = 400
    allocate(mat%psat_tab(280:400), source=pst)
    mat%use_psat = .true.
  end subroutine material


  !> The packed vapour data, as Setup_Materials builds them.
  subroutine pack(mat)
    use ICE_Config_Types_m,  only: condensed_phase_t
    use ICE_Lib_Evaporation, only: Ru, iMv, iLv, icpv, iLe, iYinf, iLvMvOverRu, iinvTboil, ialphaE
    type(condensed_phase_t), intent(inout) :: mat
    mat%ep = 0._R8
    mat%ep(iMv)         = mat%Mv
    mat%ep(iLv)         = mat%lv_al
    mat%ep(icpv)        = mat%cpv
    mat%ep(iLe)         = mat%Le
    mat%ep(iYinf)       = mat%Yinf
    mat%ep(iLvMvOverRu) = mat%lv_al*mat%Mv/Ru
    mat%ep(iinvTboil)   = 1._R8/mat%Tboil
    mat%ep(ialphaE)     = mat%alphaE
  end subroutine pack


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

end program test_vapour_sources
