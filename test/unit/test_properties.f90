!> Unit test of the property-table reader (ICE_Load_Table): the header grammar, the node and
!> column checks, Load_Table on small tables written here (Tmin > 1, permuted columns, both
!> datums), the lookup and energy inversion of ICE_Lib_Properties, and the Psat column. Prints
!> every assertion and exits non-zero if any failed.
program test_properties
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use ICE_Load_Table
  use ICE_Config_Types_m, only: obj_condensed, condensed_phase_t, obj_time_scheme
  use ICE_Global_m,       only: ICE_phase_prefix
  implicit none
  integer :: nfail = 0

  call test_tokens()
  call test_nodes()
  call test_columns()
  call test_load()
  call test_lookup()
  call test_psat()

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
    call check(check_table_nodes(T - 300._R8) == TAB_NEGATIVE_T .and. check_table_nodes(T - 280._R8) == TAB_OK, &
               'nodes: a row below 0 K is refused, a table from 0 K is not')
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
    use ICE_Lib_Properties, only: mat_rho, mat_cp

    call execute_command_line('mkdir -p INPUT')

    call write_table('u1-', '"Temperature", "Cp", "Density", "Enthalpy"', 280, 400, 1)
    call load('u1-')
    call check(obj_condensed%use_table .and. obj_condensed%T_min == 280 .and. obj_condensed%T_max == 400, &
               'load: Tmin = 280, Tmax = 400 from the rows')
    call check(lbound(obj_condensed%rho_tab, 1) == 280 .and. ubound(obj_condensed%rho_tab, 1) == 400, &
               'load: the tables are indexed by temperature')
    call check(mat_rho(obj_condensed, 250._R8) == 1500._R8 .and. mat_rho(obj_condensed, 300.4_R8) == 1500._R8 .and. &
               mat_cp(obj_condensed, 450._R8) == 2000._R8, 'load: constant properties inside and outside the range')
    call check(obj_condensed%h_datum == 'relative' .and. obj_condensed%h_off == 0._R8 .and. &
               .not. obj_condensed%rho_varies .and. .not. obj_condensed%cs_varies, 'load: relative datum, hOff = 0')

    call write_table('u2-', '"Temperature", "Enthalpy", "Density", "Cp"', 280, 400, 2)
    call load('u2-')
    call check(obj_condensed%rho_tab(300) == 1700._R8 .and. obj_condensed%rho_tab(301) == 1699._R8 .and. &
               all(obj_condensed%cs_tab == 2000._R8) .and. obj_condensed%rho_varies, &
               'load: permuted columns give the right density and cp')
    call check(abs(mat_rho(obj_condensed, 300.4_R8) - 1699.6_R8) <= 1.e-12_R8*1699.6_R8 .and. &
               mat_rho(obj_condensed, 100._R8) == 1720._R8, &
               'load: linear inside, end value outside')

    call write_table('u3-', '"Temperature", "Cp", "Density", "Enthalpy_abs"', 1, 50, 3)
    call load('u3-')
    call check(obj_condensed%h_datum == 'absolute' .and. obj_condensed%h_off == -1.5e7_R8 .and. &
               obj_condensed%T_min == 1 .and. size(obj_condensed%h_tab) == 50, 'load: absolute datum keeps its offset')
    call check(lbound(obj_condensed%e_tab, 1) == 1 .and. obj_condensed%e_tab(1) == 2000._R8 .and. &
               obj_condensed%e_tab(50) == 1.e5_R8, 'load: the energy table is h - hOff (cp*T here)')

    call write_table('u4-', '"Temperature", "Cp", "Density", "Hvap", "Enthalpy"', 280, 400, 4)
    call load('u4-')
    call check(obj_condensed%h_tab(300) == 6.e5_R8 .and. obj_condensed%h_tab(400) == 8.e5_R8 .and. &
               all(obj_condensed%rho_tab == 1500._R8), 'load: a fifth column after an unknown one binds by name')
  end subroutine test_load


  !> P2: linear lookup with saturation, and the energy inversion with both extensions.
  subroutine test_lookup()
    use ICE_Lib_Properties, only: mat_rho, mat_cp, mat_e, mat_T_from_e
    type(condensed_phase_t) :: m
    real(R8) :: T, e, worst, Ts(9)
    integer  :: i, k

    ! no table: the INI constants
    m%rho_al = 2700._R8; m%cs_al = 1598._R8
    call check(mat_rho(m, 300.4_R8) == 2700._R8 .and. mat_cp(m, 1.e4_R8) == 1598._R8, 'lookup: no table gives the constants')

    ! a varying density on 280..400 K: 2000 up to 300 K, 1000 from 301 K
    m%use_table = .true.; m%T_min = 280; m%T_max = 400
    allocate(m%rho_tab(280:400), m%cs_tab(280:400), m%h_tab(280:400), m%e_tab(280:400))
    m%rho_tab = 1000._R8; m%rho_tab(280:300) = 2000._R8; m%rho_varies = .true.
    m%cs_tab = [(1000._R8 + 2._R8*i, i = 280, 400)]; m%cs_varies = .true.
    m%h_tab  = [(1000._R8*i + real(i, R8)**2, i = 280, 400)]
    m%h_off  = m%h_tab(280) - 280._R8*(m%h_tab(281) - m%h_tab(280))
    m%e_tab  = m%h_tab - m%h_off
    call check(abs(mat_rho(m, 300.4_R8) - 1600._R8) <= 1.e-9_R8*1600._R8, 'lookup: linear between the nodes (1600 at 300.4 K)')
    call check(mat_rho(m, 300._R8) == 2000._R8 .and. mat_rho(m, 301._R8) == 1000._R8, 'lookup: exact at the nodes')
    call check(mat_rho(m, 250._R8) == 2000._R8 .and. mat_rho(m, 450._R8) == 1000._R8, 'lookup: saturates outside the table')
    T = ieee_value(1._R8, ieee_quiet_nan)
    call check(mat_rho(m, T) /= mat_rho(m, T), 'lookup: NaN in, NaN out')
    call check(abs(mat_cp(m, 300.5_R8) - 1601._R8) <= 1.e-12_R8*1601._R8, 'lookup: cp linear too')
    m%rho_tab(399) = 1._R8; m%rho_tab(400) = 9007199254740994._R8
    call check(mat_rho(m, 400._R8) == 9007199254740994._R8 .and. mat_rho(m, 450._R8) == 9007199254740994._R8, &
               'lookup: the end value exactly at and beyond Tmax, whatever the last step')
    m%rho_varies = .false.; m%rho_tab = 1500._R8
    call check(mat_rho(m, 333.3_R8) == 1500._R8, 'lookup: a constant column returns its value')
    call check(mat_rho(m, T) == 1500._R8, 'lookup: a constant column ignores T, even a NaN one')
    m%cs_varies = .false.; m%cs_tab = 2000._R8
    call check(mat_cp(m, T) == 2000._R8 .and. mat_cp(m, 1.e5_R8) == 2000._R8, 'lookup: a constant cp ignores T, even a NaN one')

    ! energy: e(0) = 0, continuous at both ends, inverse to round-off below, inside and above
    call check(mat_e(m, 0._R8) == 0._R8, 'energy: e(0) = 0 along the first segment')
    call check(mat_e(m, 280._R8) == m%e_tab(280) .and. mat_e(m, 350._R8) == m%e_tab(350) .and. &
               mat_e(m, 400._R8) == m%e_tab(400), 'energy: the nodes return eTab exactly')
    Ts = [0._R8, 1.e-3_R8, 150._R8, 279.999_R8, 280._R8, 300.4_R8, 399.5_R8, 400._R8, 520._R8]
    worst = 0._R8
    do k = 1, size(Ts)
      e = mat_e(m, Ts(k))
      worst = max(worst, abs(mat_T_from_e(m, e) - Ts(k))/max(1._R8, Ts(k)))
    enddo
    call check(worst <= 4._R8*epsilon(1._R8), 'energy: T(e(T)) = T to round-off, both extensions included')
    call check(mat_T_from_e(m, 0._R8) == 0._R8, 'energy: a vacuum cell (e = 0) gives T = 0')
    m%h_tab = [(4184._R8*i - 1.5e7_R8 + merge(15._R8, -15._R8, mod(i, 2) == 0), i = 280, 400)]
    m%h_off = m%h_tab(280) - 4184._R8*280._R8
    m%e_tab = m%h_tab - m%h_off
    e = m%e_tab(280)
    call check(abs(mat_e(m, 280._R8 - 1.e-9_R8) - e) <= 1.e-10_R8*e .and. &
               mat_T_from_e(m, e*(1._R8 - 1.e-15_R8)) <= 280._R8, &
               'energy: below Tmin the extension meets eTab(Tmin), noisy h at constant cp included')
    m%h_tab = [(1000._R8*i + real(i, R8)**2, i = 280, 400)]
    m%h_off = m%h_tab(280) - 280._R8*(m%h_tab(281) - m%h_tab(280))
    m%e_tab = m%h_tab - m%h_off
    call check(mat_e(m, 300.2_R8) < mat_e(m, 300.3_R8) .and. mat_e(m, 410._R8) > mat_e(m, 400._R8), &
               'energy: increasing inside and above')
  end subroutine test_lookup


  !> P3: the Psat checks in IGLOO's order, the lookup, and Load_Table with and without evaporation.
  subroutine test_psat()
    use ICE_Lib_Properties, only: mat_psat
    real(R8), parameter :: P0 = 101325._R8
    real(R8) :: ps(300:400), bad(300:400)
    type(condensed_phase_t) :: m
    integer  :: i

    ps = [(psat_curve(real(i, R8)), i = 300, 400)]
    call check(validate_psat_column(ps, 300, 400, 373.15_R8, P0) == PSAT_OK, 'psat: a Clausius-Clapeyron column passes')
    bad = ps; bad(350) = ieee_value(1._R8, ieee_quiet_nan)
    call check(validate_psat_column(bad, 300, 400, 373.15_R8, P0) == PSAT_NONFINITE, 'psat: a NaN is refused')
    bad = 0._R8
    call check(validate_psat_column(bad, 300, 400, 373.15_R8, P0) == PSAT_ABSENT, 'psat: all zero counts as absent')
    bad = ps; bad(300) = -1._R8
    call check(validate_psat_column(bad, 300, 400, 373.15_R8, P0) == PSAT_NEGATIVE, 'psat: a negative value is refused')
    bad = ps(400:300:-1)
    call check(validate_psat_column(bad, 300, 400, 373.15_R8, P0) == PSAT_NONMONOTONE, 'psat: a decreasing column is refused')
    bad(300) = -1._R8
    call check(validate_psat_column(bad, 300, 400, 373.15_R8, P0) == PSAT_NEGATIVE, 'psat: negative is reported before decreasing')
    bad = 1000._R8
    call check(validate_psat_column(bad, 300, 400, 373.15_R8, P0) == PSAT_CONSTANT, 'psat: a constant column is refused')
    call check(validate_psat_column(ps, 300, 400, 400._R8, P0) == PSAT_TBOIL_RANGE .and. &
               validate_psat_column(ps, 300, 400, 299.9_R8, P0) == PSAT_TBOIL_RANGE .and. &
               validate_psat_column(ps, 300, 400, ieee_value(1._R8, ieee_quiet_nan), P0) == PSAT_TBOIL_RANGE .and. &
               validate_psat_column(0.5_R8*ps, 300, 400, 399._R8, P0) == PSAT_OK .and. &
               validate_psat_column(ps, 300, 400, 300._R8, P0) /= PSAT_TBOIL_RANGE, &
               'psat: the boiling temperature must lie in [Tmin, Tmax-1]')
    call check(validate_psat_column(0.49_R8*ps, 300, 400, 373.15_R8, P0) == PSAT_TBOIL_MISMATCH .and. &
               validate_psat_column(2.01_R8*ps, 300, 400, 373.15_R8, P0) == PSAT_TBOIL_MISMATCH .and. &
               validate_psat_column(0.51_R8*ps, 300, 400, 373.15_R8, P0) == PSAT_OK .and. &
               validate_psat_column(1.99_R8*ps, 300, 400, 373.15_R8, P0) == PSAT_OK, &
               'psat: 0.5 to 2 atm at the boiling temperature')

    m%T_min = 300; m%T_max = 400
    allocate(m%psat_tab(300:400), source=ps)
    call check(abs(mat_psat(m, 350.25_R8) - (0.75_R8*ps(350) + 0.25_R8*ps(351))) <= 1.e-12_R8*ps(351) .and. &
               mat_psat(m, 250._R8) == ps(300) .and. mat_psat(m, 450._R8) == ps(400) .and. mat_psat(m, 373._R8) == ps(373), &
               'psat: linear between the nodes, exact on them, end values outside')

    call write_table('u5-', '"Temperature", "Cp", "Density", "Enthalpy", "Psat"', 300, 400, 5)
    obj_condensed%Tboil = 373.15_R8
    obj_time_scheme%evapSelect = 0
    call load('u5-')
    call check(.not. obj_condensed%use_psat .and. .not. allocated(obj_condensed%psat_tab), &
               'psat: without evaporation the column is not read')
    obj_time_scheme%evapSelect = 2
    call load('u5-')
    call check(obj_condensed%use_psat .and. lbound(obj_condensed%psat_tab, 1) == 300 .and. &
               obj_condensed%psat_tab(373) == psat_curve(373._R8) .and. obj_condensed%psat_tab(400) == psat_curve(400._R8), &
               'psat: with evaporation the column is stored, indexed by temperature')
    obj_time_scheme%evapSelect = 0
  end subroutine test_psat


  !> A Clausius-Clapeyron curve through one atmosphere at 373.15 K (water-like).
  pure function psat_curve(T) result(p)
    real(R8), intent(in) :: T
    real(R8) :: p
    p = 101325._R8*exp(-4896.8_R8*(1._R8/T - 1._R8/373.15_R8))
  end function psat_curve


  subroutine load(prefix)
    character(len=*), intent(in) :: prefix
    if (allocated(obj_condensed%rho_tab)) deallocate(obj_condensed%rho_tab)
    if (allocated(obj_condensed%cs_tab))  deallocate(obj_condensed%cs_tab)
    if (allocated(obj_condensed%h_tab))   deallocate(obj_condensed%h_tab)
    if (allocated(obj_condensed%e_tab))   deallocate(obj_condensed%e_tab)
    if (allocated(obj_condensed%psat_tab)) deallocate(obj_condensed%psat_tab)
    obj_condensed%use_table = .false.
    obj_condensed%use_psat  = .false.
    ICE_phase_prefix = prefix
    call Load_Table()
  end subroutine load


  !> kind 1: constant cp 2000, rho 1500, h = cp*T; kind 2: the same in the order T, h, rho, cp with
  !  rho = 2000 - T; kind 3: h = cp*T - 1.5e7 (absolute); kind 4: kind 1 with an unknown column before h;
  !  kind 5: kind 1 with a Psat column, psat_curve(T).
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
      case (5)
        write(unit,'(F8.1,4(1X,ES24.16))') real(T, R8), cp, rho, h, psat_curve(real(T, R8))
        cycle
      end select
      write(unit,'(F8.1,3(1X,ES24.16))') real(T, R8), cp, rho, h
    enddo
    close(unit)
  end subroutine write_table

end program test_properties
