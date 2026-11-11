!> Unit test of combustion = direct in the cell source (ICE_Mod_Sources' compute_source_), on Gallier and Godfroy's
!> aluminium (JPP 25(2) 2009, Tables 1 and 2): rho_p 2000 kg/m^3, c_p 1177 J/(kg K), Q 9.53e6 J/kg, L_v 10.8e6 J/kg,
!> T_sat 2791 K, a gas of mu 9.1e-5 Pa s, Pr 0.4, gamma 1.16, cp 2000 J/(kg K) at 3450 K and 8.71 MPa, drops of 120 um
!> with a residue of 60 um. MK closure, heat-transfer Ranz-Marshall, drag Stokes, emissivity 0.
!>   C1 a drop at T_sat 5 m/s slower than the gas (Re ~ 60): the mass rate of their eq. 4, typed here,
!>      2 pi (mu/Pr) n D ln(1 + B) (1 + 0.3 Re^1/2 Pr^1/3), B of their eq. 5, to 1e-13
!>   C2 the drop keeps its temperature: the source's internal-energy part is zero, to round-off of the energy
!>   C3 the gas receives the heat of reaction of the burnt mass, n mdot Q, apart from the source
!>   C4 a drop below T_sat, a drop of the residue's diameter, and a gas so cold that B <= 0: the source of the same
!>      material without combustion, bit for bit, and no heat to the gas
!>   C5 a material with emissivity 0.5 burns at the same rate and still keeps its temperature
!> Prints every assertion and exits non-zero if any failed.
program test_combustion
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  implicit none
  real(R8), parameter :: PI = acos(-1._R8)
  real(R8), parameter :: RHO_M = 2000._R8, CP_P = 1177._R8, Q = 9.53e6_R8, LV = 10.8e6_R8, TSAT = 2791._R8
  real(R8), parameter :: MU = 9.1e-5_R8, PR = 0.4_R8, GAM = 1.16_R8, CPG = 2000._R8, TG = 3450._R8, PG = 8.71e6_R8
  real(R8), parameter :: D0 = 120.e-6_R8, DRES = 60.e-6_R8, SLIP = 5._R8, RHO_P = 1._R8
  real(R8), parameter :: RG = CPG*(GAM - 1._R8)/GAM, KG = MU*CPG/PR
  logical :: ok

  ok = .true.
  call setup_model()
  call burning(ok)
  call inert(ok)
  if (.not. ok) then
    write(*,'(A)') 'test_combustion: FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_combustion: PASS'

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


  !> Family 1 is MK, drag Stokes, heat-transfer Ranz-Marshall, no blowing correction.
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


  !> Gallier's aluminium, burning (comb = 1) or not.
  subroutine material(mat, comb, emiss)
    use ICE_Config_Types_m, only: condensed_phase_t
    type(condensed_phase_t), intent(out) :: mat
    integer,                 intent(in)  :: comb
    real(R8),                intent(in)  :: emiss
    mat%rho_al = RHO_M; mat%cs_al = CP_P; mat%lv_al = LV; mat%emiss = emiss; mat%Tboil = TSAT
    mat%combSelect = comb; mat%qComb = Q; mat%dRes = DRES
  end subroutine material


  !> The MK state of drops of diameter d at temperature T, at rest, bulk density RHO_P.
  pure function drops(d, T) result(c)
    real(R8), intent(in) :: d, T
    real(R8) :: c(6)
    c = [RHO_P, 0._R8, 0._R8, 0._R8, T, RHO_P/(RHO_M*PI/6._R8*d**3)]
  end function drops


  !> The rate of the internal energy per unit volume the source gives, rho de/dt = S5 - u.Sm + (u^2/2 - e) S1
  pure function internal(c, s) result(r)
    real(R8), intent(in) :: c(6), s(6)
    real(R8) :: r
    r = s(5) - dot_product(c(2:4), s(2:4)) + (0.5_R8*dot_product(c(2:4), c(2:4)) - CP_P*c(5))*s(1)
  end function internal


  subroutine burning(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    use ICE_Mod_Sources,    only: compute_source_
    logical, intent(inout) :: ok
    type(condensed_phase_t) :: mat
    real(R8) :: c(6), gas(5), s(6), tau, heat, rhog, Re, B, w, ref, scale
    integer  :: e

    rhog = PG/(RG*TG)
    gas  = [rhog, SLIP, 0._R8, 0._R8, TG]
    c    = drops(D0, TSAT)
    Re   = rhog*SLIP*D0/MU
    B    = (CPG*(TG - TSAT) + Q)/LV
    ref  = 2._R8*PI*(MU/PR)*c(6)*D0*log(1._R8 + B)*(1._R8 + 0.3_R8*sqrt(Re)*PR**(1._R8/3._R8))
    write(*,'(A,F8.3,A,F8.5,A,ES12.5,A)') '   Re ', Re, ', B ', B, ', eq. 4: ', ref, ' kg/(m^3 s)'

    do e = 1, 2
      call material(mat, 1, merge(0._R8, 0.5_R8, e == 1))
      tau = 0._R8
      call compute_source_(c, tau, gas, RG, GAM, KG, MU, s, mat, 6, gas_heat=heat)
      w = -s(1)
      scale = w*(CP_P*TSAT + LV)
      if (e == 1) then
        write(*,'(A,ES24.16,A,ES10.3)') '   C1: -S1 ', w, ', relative to eq. 4 ', (w - ref)/ref
        call check(abs(w - ref) <= 1.e-13_R8*ref, 'C1 the burn rate is Gallier''s eq. 4 with eq. 5, to 1e-13', ok)
        write(*,'(A,ES10.3,A)') '   C2: internal-energy rate ', internal(c, s)/scale, ' of n mdot (c_p T + L_v)'
        call check(abs(internal(c, s)) <= 1.e-14_R8*scale, 'C2 the burning drop keeps its temperature', ok)
        write(*,'(A,ES24.16)') '   C3: gas heat / (n mdot Q) - 1 ', heat/(w*Q) - 1._R8
        call check(abs(heat - w*Q) <= 1.e-15_R8*w*Q, 'C3 the gas receives n mdot Q', ok)
      else
        call check(abs(w - ref) <= 1.e-13_R8*ref .and. abs(internal(c, s)) <= 1.e-14_R8*scale .and. &
                   abs(heat - w*Q) <= 1.e-15_R8*w*Q, &
                   'C5 with emissivity 0.5: the same rate, the temperature kept, n mdot Q to the gas', ok)
      endif
    enddo
  end subroutine burning


  !> C4: where the drop does not burn, the source is the inert material's.
  subroutine inert(ok)
    use ICE_Config_Types_m, only: condensed_phase_t
    use ICE_Mod_Sources,    only: compute_source_
    logical, intent(inout) :: ok
    type(condensed_phase_t) :: burn, none
    real(R8) :: c(6), gas(5), sa(6), sb(6), ta, tb, heat, heat_none
    integer  :: leg
    character(len=40) :: what(3)
    what = [character(len=40) :: 'below T_sat (2700 K)', 'just below the residue diameter', 'in a gas where B <= 0 (300 K, Q 0)']

    call material(none, 0, 0._R8)
    do leg = 1, 3
      call material(burn, 1, 0._R8)
      gas = [PG/(RG*TG), SLIP, 0._R8, 0._R8, TG]
      select case (leg)
      case (1); c = drops(D0, 2700._R8)
      case (2); c = drops(0.999_R8*DRES, TSAT)
      case (3); c = drops(D0, TSAT); gas(5) = 300._R8; gas(1) = PG/(RG*300._R8); burn%qComb = 0._R8
      end select
      ta = 0._R8; tb = 0._R8; heat = 1._R8
      call compute_source_(c, ta, gas, RG, GAM, KG, MU, sa, burn, 6, gas_heat=heat)
      call compute_source_(c, tb, gas, RG, GAM, KG, MU, sb, none, 6, heat_none)
      call check(all(sa == sb) .and. ta == tb .and. heat == 0._R8 .and. heat_none == 0._R8, &
                 'C4 '//trim(what(leg))//': the inert source bit for bit, no heat to the gas', ok)
    enddo
  end subroutine inert

end program test_combustion
