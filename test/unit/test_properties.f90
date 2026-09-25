!> Unit test of the property-table reader (ICE_Load_Table): the header grammar, the node and
!> column checks, and Load_Table on small tables written here (Tmin > 1, permuted columns, both
!> datums). Prints every assertion and exits non-zero if any failed.
program test_properties
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use ICE_Load_Table
  use ICE_Config_Types_m, only: obj_condensed
  use ICE_Global_m,       only: ICE_phase_prefix
  implicit none
  integer :: nfail = 0

  call test_tokens()
  call test_nodes()
  call test_columns()
  call test_load()

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_properties: ', nfail, ' FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_properties: PASS'

contains

  subroutine check(ok, what)
    logical,          intent(in) :: ok
    character(len=*), intent(in) :: what
    if (ok) then
      write(*,'(A)') '   OK   '//what
    else
      write(*,'(A)') '   FAIL '//what
      nfail = nfail + 1
    endif
  end subroutine check


  subroutine test_tokens()
    character(len=64), allocatable :: tk(:)
    integer :: icp, irho, ih, ips, code
    logical :: rel

    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_OK .and. icp == 1 .and. irho == 2 .and. ih == 3 .and. ips == 0 .and. rel, &
               'tokens: the ATLAS order binds Cp, Density, Enthalpy to vars 1, 2, 3, relative')
    tk = [character(len=64) :: 'Temperature', 'Enthalpy', 'Density', 'Cp']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_OK .and. icp == 3 .and. irho == 2 .and. ih == 1, 'tokens: permuted columns bind by name')
    tk = [character(len=64) :: 'temperature', 'cp', 'density', 'enthalpy_abs', 'psat']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_OK .and. .not. rel .and. ips == 4, 'tokens: lower-case aliases, Enthalpy_abs is absolute, Psat found')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Hvap', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_OK .and. irho == 3 .and. ih == 4, 'tokens: an unknown name is skipped, not misbound')
    tk = [character(len=64) :: 'Temperature', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_NO_CP, 'tokens: no Cp is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_NO_DENSITY, 'tokens: no Density is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_NO_ENTHALPY, 'tokens: no enthalpy is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Enthalpy_abs']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_TWO_ENTHALPY, 'tokens: both enthalpies are refused')
    tk = [character(len=64) :: 'Cp', 'Temperature', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_NO_TEMPERATURE, 'tokens: Temperature must come first')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Cp', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_DUPLICATE, 'tokens: a column named twice is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code)
    call check(code == TAB_DUPLICATE, 'tokens: one enthalpy named twice is a duplicate, not two datums')
  end subroutine test_tokens


  subroutine test_nodes()
    real(R8) :: T(121)
    integer  :: i
    T = [(real(279+i, R8), i = 1, 121)]
    call check(check_table_nodes(T) == TAB_OK, 'nodes: 280..400 K accepted')
    call check(check_table_nodes(T(1:1)) == TAB_FEW_ROWS, 'nodes: one row is refused')
    call check(check_table_nodes(T + 0.5_R8) == TAB_OFF_NODE, 'nodes: a half-kelvin offset is refused')
    T(60) = T(60) + 1._R8
    call check(check_table_nodes(T) == TAB_OFF_NODE, 'nodes: a missing kelvin is refused')
    T(60) = ieee_value(1._R8, ieee_quiet_nan)
    call check(check_table_nodes(T) == TAB_NONFINITE, 'nodes: a NaN temperature is refused')
  end subroutine test_nodes


  subroutine test_columns()
    real(R8) :: T(121), cp(121), rho(121), h(121)
    integer  :: i
    T = [(real(279+i, R8), i = 1, 121)]
    cp = 2000._R8; rho = 1500._R8; h = cp*T
    call check(check_table_columns(T, cp, rho, h, .true.) == TAB_OK, 'columns: constant cp, h = cp*T, relative')
    call check(abs(table_h_offset(T, cp, h)) == 0._R8, 'columns: hOff = 0 for h = cp*T')
    call check(check_table_columns(T, cp, rho, h - 1.5e7_R8, .false.) == TAB_OK, 'columns: absolute datum with an offset')
    call check(table_h_offset(T, cp, h - 1.5e7_R8) == -1.5e7_R8, 'columns: hOff = the absolute offset')
    call check(check_table_columns(T, cp, rho, h + 1.e5_R8, .true.) == TAB_DATUM_MISMATCH, &
               'columns: a relative header with an offset is refused')
    h(60) = h(60) + 100._R8
    call check(check_table_columns(T, cp, rho, h, .true.) == TAB_H_CP_MISMATCH, 'columns: a kink in h at constant cp')
    cp = 1000._R8 + 2._R8*T; h = 1000._R8*T + T*T - 1.58e7_R8
    h(60:) = h(60:) + 9._R8
    call check(check_table_columns(T, cp, rho, h, .false.) == TAB_OK, &
               'columns: a 9 J/kg seam in an absolute h (1e-6 of |h|) is accepted')
    h(60:) = h(60:) + 1.e5_R8
    call check(check_table_columns(T, cp, rho, h, .false.) == TAB_H_CP_MISMATCH, 'columns: a latent jump in h is refused')
    cp = 1000._R8 + 2._R8*T; h = 1000._R8*T + T*T
    call check(check_table_columns(T, cp, rho, h, .true.) == TAB_OK, 'columns: cp = 1000 + 2T with its exact integral')
    call check(table_h_offset(T, cp, h) == h(1) - T(1)*(h(2) - h(1)), 'columns: hOff from the first segment when cp varies')
    call check(check_table_columns(T, cp, rho, 1.01_R8*h, .true.) == TAB_H_CP_MISMATCH, 'columns: h 1 % off its cp')
    rho(10) = 0._R8
    call check(check_table_columns(T, cp, rho, h, .true.) == TAB_RHO_NONPOSITIVE, 'columns: zero density refused')
    rho = 1500._R8; cp(10) = -1._R8
    call check(check_table_columns(T, cp, rho, h, .true.) == TAB_CP_NONPOSITIVE, 'columns: negative cp refused')
    cp = 1000._R8 + 2._R8*T; h(20) = h(19)
    call check(check_table_columns(T, cp, rho, h, .true.) == TAB_H_NONMONOTONE, 'columns: flat h refused')
    h = 1000._R8*T + T*T; rho(5) = ieee_value(1._R8, ieee_quiet_nan)
    call check(check_table_columns(T, cp, rho, h, .true.) == TAB_NONFINITE, 'columns: NaN density refused')
  end subroutine test_columns


  !> Load_Table on tables written here, one per prefix.
  subroutine test_load()

    call execute_command_line('mkdir -p INPUT')

    call write_table('u1-', '"Temperature", "Cp", "Density", "Enthalpy"', 280, 400, 1)
    call load('u1-')
    call check(obj_condensed%use_table .and. obj_condensed%T_min == 280 .and. obj_condensed%T_max == 400, &
               'load: Tmin = 280, Tmax = 400 from the rows')
    call check(lbound(obj_condensed%rho_tab, 1) == 280 .and. ubound(obj_condensed%rho_tab, 1) == 400, &
               'load: the tables are indexed by temperature')
    call check(get_rho_al(250._R8) == 1500._R8 .and. get_rho_al(300.4_R8) == 1500._R8 .and. &
               get_cs_al(450._R8) == 2000._R8, 'load: constant properties inside and outside the range')
    call check(obj_condensed%h_datum == 'relative' .and. obj_condensed%h_off == 0._R8 .and. &
               .not. obj_condensed%rho_varies .and. .not. obj_condensed%cs_varies, 'load: relative datum, hOff = 0')

    call write_table('u2-', '"Temperature", "Enthalpy", "Density", "Cp"', 280, 400, 2)
    call load('u2-')
    call check(obj_condensed%rho_tab(300) == 1700._R8 .and. obj_condensed%rho_tab(301) == 1699._R8 .and. &
               all(obj_condensed%cs_tab == 2000._R8) .and. obj_condensed%rho_varies, &
               'load: permuted columns give the right density and cp')
    call check(get_rho_al(300.4_R8) == 1700._R8 .and. get_rho_al(100._R8) == 1720._R8, &
               'load: nearest kelvin inside, end value outside')

    call write_table('u3-', '"Temperature", "Cp", "Density", "Enthalpy_abs"', 1, 50, 3)
    call load('u3-')
    call check(obj_condensed%h_datum == 'absolute' .and. obj_condensed%h_off == -1.5e7_R8 .and. &
               obj_condensed%T_min == 1 .and. size(obj_condensed%h_tab) == 50, 'load: absolute datum keeps its offset')

    call write_table('u4-', '"Temperature", "Cp", "Density", "Hvap", "Enthalpy"', 280, 400, 4)
    call load('u4-')
    call check(obj_condensed%h_tab(300) == 6.e5_R8 .and. obj_condensed%h_tab(400) == 8.e5_R8 .and. &
               all(obj_condensed%rho_tab == 1500._R8), 'load: a fifth column after an unknown one binds by name')
  end subroutine test_load


  subroutine load(prefix)
    character(len=*), intent(in) :: prefix
    if (allocated(obj_condensed%rho_tab)) deallocate(obj_condensed%rho_tab)
    if (allocated(obj_condensed%cs_tab))  deallocate(obj_condensed%cs_tab)
    if (allocated(obj_condensed%h_tab))   deallocate(obj_condensed%h_tab)
    obj_condensed%use_table = .false.
    ICE_phase_prefix = prefix
    call Load_Table()
  end subroutine load


  !> kind 1: constant cp 2000, rho 1500, h = cp*T; kind 2: the same in the order T, h, rho, cp with
  !  rho = 2000 - T; kind 3: h = cp*T - 1.5e7 (absolute); kind 4: kind 1 with an unknown column before h.
  subroutine write_table(prefix, header, Tmin, Tmax, kind)
    character(len=*), intent(in) :: prefix, header
    integer,          intent(in) :: Tmin, Tmax, kind
    integer  :: unit, T
    real(R8) :: cp, rho, h
    open(newunit=unit, file='INPUT/'//prefix//'properties.dat', status='replace', action='write')
    write(unit,'(A)') 'TITLE = "unit test"'
    write(unit,'(A)') 'VARIABLES = '//header
    write(unit,'(A)') 'ZONE T="A"'
    write(unit,'(A,I0,A)') 'I=', Tmax - Tmin + 1, ', F=POINT'
    do T = Tmin, Tmax
      cp = 2000._R8; rho = 1500._R8; h = cp*T
      select case (kind)
      case (2)
        rho = 2000._R8 - T
        write(unit,'(F8.1,3(1X,ES24.16))') real(T, R8), h, rho, cp
        cycle
      case (3)
        h = h - 1.5e7_R8
      case (4)
        write(unit,'(F8.1,4(1X,ES24.16))') real(T, R8), cp, rho, 3.3e5_R8, h
        cycle
      end select
      write(unit,'(F8.1,3(1X,ES24.16))') real(T, R8), cp, rho, h
    enddo
    close(unit)
  end subroutine write_table

end program test_properties
