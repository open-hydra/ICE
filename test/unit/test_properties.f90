!> Unit test of the property-table reader (ICE_Load_Table): the header grammar, the node and
!> column checks, Load_Table on small tables written here (Tmin > 1, permuted columns, both
!> datums), the lookup and energy inversion of ICE_Lib_Properties, the Psat column, the Lv, Mv
!> and Cpv columns (their checks, and the ATLAS water table INPUT/aw-properties.dat), and the
!> energy the MK, IG and AG closures store, recover and lose with the evaporated mass. Prints
!> every assertion and exits non-zero if any failed.
program test_properties
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use ICE_Load_Table
  use ICE_Config_Types_m, only: obj_condensed, condensed_phase_t
  use ICE_Global_m,       only: ICE_phase_prefix, nmat
  implicit none
  integer :: nfail = 0

  nmat = 1
  allocate(obj_condensed(1))
  call test_tokens()
  call test_nodes()
  call test_columns()
  call test_load()
  call test_lookup()
  call test_psat()
  call test_vapour()
  call test_closures()

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
    integer :: icp, irho, ih, ips, code, ilv, imv, icpv, isig, imu
    logical :: rel

    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. icp == 1 .and. irho == 2 .and. ih == 3 .and. ips == 0 .and. rel, &
               'tokens: the ATLAS order binds Cp, Density, Enthalpy to vars 1, 2, 3, relative')
    tk = [character(len=64) :: 'Temperature', 'Enthalpy', 'Density', 'Cp']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. icp == 3 .and. irho == 2 .and. ih == 1, 'tokens: permuted columns bind by name')
    tk = [character(len=64) :: 'temperature', 'cp', 'density', 'enthalpy_abs', 'psat']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. .not. rel .and. ips == 4, 'tokens: lower-case aliases, Enthalpy_abs is absolute, Psat found')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Hvap', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. irho == 3 .and. ih == 4, 'tokens: an unknown name is skipped, not misbound')
    call check(ilv == 0 .and. imv == 0 .and. icpv == 0, 'tokens: Hvap is not bound as a vapour column')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy_abs', 'Psat', 'Lv', 'Mv', 'Cpv']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. ips == 4 .and. ilv == 5 .and. imv == 6 .and. icpv == 7 .and. .not. rel, &
               'tokens: the ATLAS order binds Psat, Lv, Mv, Cpv to vars 4, 5, 6, 7')
    tk = [character(len=64) :: 'temperature', 'CPV', 'mv', 'Density', 'LV', 'cp', 'enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. icpv == 1 .and. imv == 2 .and. ilv == 4 .and. icp == 5 .and. ips == 0, &
               'tokens: Lv, Mv, Cpv and their aliases bind by name in any order, without Psat')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Psat', 'Lv', 'lv']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_DUPLICATE, 'tokens: Lv named twice is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Psat', 'MV', 'Mv']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_DUPLICATE, 'tokens: Mv named twice is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Psat', 'cpv', 'Cpv']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_DUPLICATE, 'tokens: Cpv named twice is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'VISCOSITY', 'sigma']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. isig == 5 .and. imu == 4 .and. ips == 0 .and. ilv == 0 .and. icpv == 0, &
               'tokens: VISCOSITY and sigma bind their columns by name')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_OK .and. isig == 0 .and. imu == 0, 'tokens: no Sigma or Viscosity name: both 0')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Sigma', 'SIGMA']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_DUPLICATE, 'tokens: Sigma named twice is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'viscosity', 'Viscosity']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_DUPLICATE, 'tokens: Viscosity named twice is refused')
    tk = [character(len=64) :: 'Temperature', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_NO_CP, 'tokens: no Cp is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_NO_DENSITY, 'tokens: no Density is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_NO_ENTHALPY, 'tokens: no enthalpy is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Enthalpy_abs']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_TWO_ENTHALPY, 'tokens: both enthalpies are refused')
    tk = [character(len=64) :: 'Cp', 'Temperature', 'Density', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_NO_TEMPERATURE, 'tokens: Temperature must come first')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Cp', 'Enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
    call check(code == TAB_DUPLICATE, 'tokens: a column named twice is refused')
    tk = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'enthalpy']
    call classify_table_tokens(tk, icp, irho, ih, ips, rel, code, ilv, imv, icpv, isig, imu)
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
    call check(obj_condensed(1)%use_table .and. obj_condensed(1)%T_min == 280 .and. obj_condensed(1)%T_max == 400, &
               'load: Tmin = 280, Tmax = 400 from the rows')
    call check(lbound(obj_condensed(1)%rho_tab, 1) == 280 .and. ubound(obj_condensed(1)%rho_tab, 1) == 400, &
               'load: the tables are indexed by temperature')
    call check(mat_rho(obj_condensed(1), 250._R8) == 1500._R8 .and. mat_rho(obj_condensed(1), 300.4_R8) == 1500._R8 .and. &
               mat_cp(obj_condensed(1), 450._R8) == 2000._R8, 'load: constant properties inside and outside the range')
    call check(obj_condensed(1)%h_datum == 'relative' .and. obj_condensed(1)%h_off == 0._R8 .and. &
               .not. obj_condensed(1)%rho_varies .and. .not. obj_condensed(1)%cs_varies, 'load: relative datum, hOff = 0')

    call write_table('u2-', '"Temperature", "Enthalpy", "Density", "Cp"', 280, 400, 2)
    call load('u2-')
    call check(obj_condensed(1)%rho_tab(300) == 1700._R8 .and. obj_condensed(1)%rho_tab(301) == 1699._R8 .and. &
               all(obj_condensed(1)%cs_tab == 2000._R8) .and. obj_condensed(1)%rho_varies, &
               'load: permuted columns give the right density and cp')
    call check(abs(mat_rho(obj_condensed(1), 300.4_R8) - 1699.6_R8) <= 1.e-12_R8*1699.6_R8 .and. &
               mat_rho(obj_condensed(1), 100._R8) == 1720._R8, &
               'load: linear inside, end value outside')

    call write_table('u3-', '"Temperature", "Cp", "Density", "Enthalpy_abs"', 1, 50, 3)
    call load('u3-')
    call check(obj_condensed(1)%h_datum == 'absolute' .and. obj_condensed(1)%h_off == -1.5e7_R8 .and. &
               obj_condensed(1)%T_min == 1 .and. size(obj_condensed(1)%h_tab) == 50, 'load: absolute datum keeps its offset')
    call check(lbound(obj_condensed(1)%e_tab, 1) == 1 .and. obj_condensed(1)%e_tab(1) == 2000._R8 .and. &
               obj_condensed(1)%e_tab(50) == 1.e5_R8, 'load: the energy table is h - hOff (cp*T here)')

    call write_table('u4-', '"Temperature", "Cp", "Density", "Hvap", "Enthalpy"', 280, 400, 4)
    call load('u4-')
    call check(obj_condensed(1)%h_tab(300) == 6.e5_R8 .and. obj_condensed(1)%h_tab(400) == 8.e5_R8 .and. &
               all(obj_condensed(1)%rho_tab == 1500._R8), 'load: a fifth column after an unknown one binds by name')
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
    obj_condensed(1)%Tboil = 373.15_R8
    obj_condensed(1)%evapSelect = 0
    call load('u5-')
    call check(.not. obj_condensed(1)%use_psat .and. .not. allocated(obj_condensed(1)%psat_tab), &
               'psat: without evaporation the column is not read')
    obj_condensed(1)%evapSelect = 2
    call load('u5-')
    call check(obj_condensed(1)%use_psat .and. lbound(obj_condensed(1)%psat_tab, 1) == 300 .and. &
               obj_condensed(1)%psat_tab(373) == psat_curve(373._R8) .and. obj_condensed(1)%psat_tab(400) == psat_curve(400._R8), &
               'psat: with evaporation the column is stored, indexed by temperature')
    obj_condensed(1)%evapSelect = 0
  end subroutine test_psat


  !> The Lv, Mv, Cpv columns: the checks in IGLOO's order, the slope band, the ATLAS water table.
  subroutine test_vapour()
    call vapour_codes()
    call vapour_band()
    call vapour_atlas()
  end subroutine test_vapour


  !> Each code of validate_vapour_columns once, on the column sets of IGLOO's reader test: a Psat
  !  column with ln p = -5000/T + const, so Lv*Mv/Ru = 5000 K is its Clausius-Clapeyron slope.
  subroutine vapour_codes()
    use ICE_Lib_Evaporation, only: Ru
    real(R8) :: ps(100:110), Lv(100:110), Mv(100:110), cv(100:110), zero(100:110), bad(100:110)
    integer  :: k

    do k = 100, 110
      ps(k) = 101325._R8*exp(-5000._R8*(1._R8/real(k, R8) - 1._R8/105._R8))
    enddo
    Lv = 5000._R8*Ru/18._R8; Mv = 18._R8; cv = 1900._R8; zero = 0._R8
    call check(validate_vapour_columns(Lv, Mv, cv, ps, .true., 100, 110, 105._R8, Ru, 0._R8) == VAP_OK, &
               'vapour: the three columns on the Psat slope pass')
    call check(validate_vapour_columns(zero, zero, zero, ps, .true., 100, 110, 105._R8, Ru, 18._R8) == VAP_ABSENT, &
               'vapour: all three zero count as absent')
    bad = Lv; bad(103) = ieee_value(1._R8, ieee_quiet_nan)
    call check(validate_vapour_columns(bad, Mv, cv, ps, .true., 100, 110, 105._R8, Ru, 0._R8) == VAP_NONFINITE, &
               'vapour: a NaN in Lv is refused')
    bad = cv; bad(100:102) = 0._R8
    call check(validate_vapour_columns(Lv, Mv, bad, ps, .true., 100, 110, 105._R8, Ru, 0._R8) == VAP_PARTIAL, &
               'vapour: Cpv zero on three rows is refused')
    call check(validate_vapour_columns(Lv, Mv, -cv, ps, .true., 100, 110, 105._R8, Ru, 0._R8) == VAP_NONPOSITIVE, &
               'vapour: a negative Cpv is refused')
    bad = Mv; bad(104) = 18.5_R8
    call check(validate_vapour_columns(Lv, bad, cv, ps, .true., 100, 110, 105._R8, Ru, 0._R8) == VAP_MV_VARIES, &
               'vapour: an Mv that varies is refused')
    call check(validate_vapour_columns(Lv, Mv, cv, zero, .false., 100, 110, 105._R8, Ru, 0._R8) == VAP_WITHOUT_PSAT, &
               'vapour: the columns without a Psat column are refused')
    call check(validate_vapour_columns(Lv/1000._R8, Mv, cv, ps, .true., 100, 110, 105._R8, Ru, 0._R8) == &
               VAP_CC_MISMATCH, 'vapour: Lv in kJ/kg is refused (not the slope of Psat)')
    call check(validate_vapour_columns(Lv, zero, cv, ps, .true., 100, 110, 105._R8, Ru, 18._R8) == VAP_OK .and. &
               validate_vapour_columns(Lv, zero, cv, ps, .true., 100, 110, 105._R8, Ru, 18._R8/2.5_R8) == &
               VAP_CC_MISMATCH, 'vapour: without an Mv column the slope takes the key''s Mv (18 passes, 18/2.5 not)')
  end subroutine vapour_codes


  !> The slope check's band on psat_curve, a constant Lv = f*4896.8*Ru/Mv: the ratio is 1/f.
  subroutine vapour_band()
    use ICE_Lib_Evaporation, only: Ru
    real(R8), parameter :: f(5) = [0.49_R8, 0.51_R8, 1.99_R8, 2.01_R8, 1.e-3_R8]
    integer,  parameter :: expect(5) = [VAP_CC_MISMATCH, VAP_OK, VAP_OK, VAP_CC_MISMATCH, VAP_CC_MISMATCH]
    real(R8) :: ps(300:400), Lv(300:400), Mv(300:400), cv(300:400)
    character(len=96) :: what
    integer  :: i, k

    ps = [(psat_curve(real(i, R8)), i = 300, 400)]
    Mv = 18.015_R8; cv = 1900._R8
    do k = 1, size(f)
      Lv = f(k)*4896.8_R8*Ru/18.015_R8
      write(what, '(A,ES8.2,A,ES11.5,A)') 'vapour: Lv x ', f(k), ' (slope ratio ', 1._R8/f(k), ')'
      call check(validate_vapour_columns(Lv, Mv, cv, ps, .true., 300, 400, 373.15_R8, Ru, 0._R8) == expect(k), &
                 trim(what)//trim(merge(' accepted', ' refused ', expect(k) == VAP_OK)))
    enddo
  end subroutine vapour_band


  !> Load_Table on the table ATLAS writes for water (Psat, Lv, Mv, Cpv from the H2O(L)/H2O pair):
  !  not read without evaporation; with it the columns are the file's, the scalars Lv(298.15 K),
  !  the column's Mv and Cpv(boiling-temperature), and mat_lv, mat_cpv at six temperatures (below,
  !  the datum, a node, between nodes, the boiling temperature, above) equal the file's numbers by
  !  the lookup's formula: exact on a node and outside, 4 eps between nodes.
  subroutine vapour_atlas()
    use ICE_Lib_Properties, only: mat_lv, mat_cpv
    real(R8), parameter :: Tq(6) = [250._R8, 298.15_R8, 300._R8, 337.5_R8, 373.15_R8, 400._R8]
    logical,  parameter :: onNode(6) = [.true., .false., .true., .false., .false., .true.]
    real(R8) :: T(200), Lvf(200), Mvf(200), Cpvf(200), row(8), eLv, eCpv, tol
    integer  :: u, ios, n, k, lo, hi
    character(len=512) :: line
    character(len=96)  :: what

    obj_condensed(1)%Tboil = 373.15_R8
    obj_condensed(1)%evapSelect = 0
    call load('aw-')
    call check(.not. obj_condensed(1)%use_lv .and. .not. obj_condensed(1)%use_cpv .and. &
               .not. allocated(obj_condensed(1)%lv_tab) .and. .not. allocated(obj_condensed(1)%cpv_tab) .and. &
               obj_condensed(1)%lv_al == 10.8e6_R8 .and. obj_condensed(1)%Mv == 26.98_R8, &
               'vapour: without evaporation the Lv, Mv, Cpv columns are not read')
    obj_condensed(1)%evapSelect = 2
    call load('aw-')
    call check(obj_condensed(1)%use_psat .and. obj_condensed(1)%use_lv .and. obj_condensed(1)%use_cpv, &
               'vapour: with evaporation the ATLAS water table gives Psat, Lv and Cpv')
    if (.not. (allocated(obj_condensed(1)%lv_tab) .and. allocated(obj_condensed(1)%cpv_tab))) then
      call check(.false., 'vapour: the Lv and Cpv tables are allocated')
      obj_condensed(1)%evapSelect = 0
      return
    endif

    open(newunit=u, file='INPUT/aw-properties.dat', status='old', action='read')
    n = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      read(line, *, iostat=ios) row
      if (ios /= 0) cycle
      n = n + 1
      T(n) = row(1); Lvf(n) = row(6); Mvf(n) = row(7); Cpvf(n) = row(8)
    enddo
    close(u)
    lo = obj_condensed(1)%T_min; hi = obj_condensed(1)%T_max
    call check(n == 101 .and. lo == 280 .and. hi == 380 .and. nint(T(1)) == lo, &
               'vapour: the table spans 280..380 K in 101 rows')
    call check(lbound(obj_condensed(1)%lv_tab, 1) == lo .and. ubound(obj_condensed(1)%lv_tab, 1) == hi .and. &
               all(obj_condensed(1)%lv_tab == Lvf(1:n)) .and. all(obj_condensed(1)%cpv_tab == Cpvf(1:n)), &
               'vapour: the Lv and Cpv tables are the file''s columns, indexed by temperature')
    call check(obj_condensed(1)%Mv == 18.015_R8 .and. all(Mvf(1:n) == 18.015_R8), 'vapour: Mv = 18.015 from its column')
    call check(rel(obj_condensed(1)%lv_al, value_at(Lvf(1:n), lo, hi, 298.15_R8)) <= 4._R8*epsilon(1._R8), &
               'vapour: the scalar latent heat is Lv(298.15 K)')
    call check(rel(obj_condensed(1)%cpv, value_at(Cpvf(1:n), lo, hi, 373.15_R8)) <= 4._R8*epsilon(1._R8), &
               'vapour: the scalar vapour cp is Cpv(boiling-temperature)')
    do k = 1, size(Tq)
      eLv  = rel(mat_lv(obj_condensed(1), Tq(k)),  value_at(Lvf(1:n), lo, hi, Tq(k)))
      eCpv = rel(mat_cpv(obj_condensed(1), Tq(k)), value_at(Cpvf(1:n), lo, hi, Tq(k)))
      tol  = merge(0._R8, 4._R8*epsilon(1._R8), onNode(k))
      write(what, '(A,F6.2,A,2ES9.2,A)') 'vapour: Lv and Cpv at ', Tq(k), ' K are the file''s (rel.', eLv, eCpv, ')'
      call check(eLv <= tol .and. eCpv <= tol, trim(what))
    enddo
    obj_condensed(1)%evapSelect = 0
  end subroutine vapour_atlas


  !> A column at T: linear between the nodes, the end value outside (the lookup's rule, written here).
  pure function value_at(tab, lo, hi, T) result(v)
    integer,  intent(in) :: lo, hi
    real(R8), intent(in) :: tab(lo:hi), T
    real(R8) :: v, Tc
    integer  :: i
    Tc = min(max(T, real(lo, R8)), real(hi, R8))
    if (Tc >= real(hi, R8)) then
      v = tab(hi)
    else
      i = int(Tc)
      v = tab(i) + (tab(i+1) - tab(i))*(Tc - real(i, R8))
    endif
  end function value_at


  pure function rel(a, b) result(e)
    real(R8), intent(in) :: a, b
    real(R8) :: e
    e = abs(a - b)/abs(b)
  end function rel


  !> P6: a tabulated cp stores rho*e(T), e = h - hOff, and recovers T to round-off, in the
  !  conserved energy, the flux and the source of every closure; an empty cell recovers the T a
  !  constant cp gives it; a constant cp gives the same bits with or without a table.
  subroutine test_closures()
    use ICE_Lib_MK, only: prim_2_cons_MK, cons_2_prim_MK, flux_make_MK, source_make_MK
    use ICE_Lib_IG, only: prim_2_cons_IG, cons_2_prim_IG, flux_make_IG, source_make_IG
    use ICE_Lib_AG, only: prim_2_cons_AG, cons_2_prim_AG, flux_make_AG, source_make_AG
    use ICE_Lib_Properties, only: mat_e
    type(condensed_phase_t) :: m, c, t
    real(R8), parameter :: Ts(6) = [1.e-3_R8, 150._R8, 280._R8, 333.3_R8, 400._R8, 520._R8]
    real(R8), parameter :: nrm(3) = [0.6_R8, 0.8_R8, 0._R8], lv = 2.3e6_R8
    real(R8) :: pk(6), pi7(7), pa(12), q(12), b(12), f(12), s(12), force(6), e, ke, un, tr, Tv(2)
    real(R8) :: wr, wc, wf, ws
    integer  :: i, k

    ! cp = 1000 + 2T on 280..400 K, h its exact integral, anchored at 0 K along the first segment
    m%use_table = .true.; m%cs_varies = .true.; m%T_min = 280; m%T_max = 400; m%lv_al = lv
    allocate(m%rho_tab(280:400), source=1500._R8)
    allocate(m%cs_tab(280:400), m%h_tab(280:400))
    m%cs_tab = [(1000._R8 + 2._R8*i, i = 280, 400)]
    m%h_tab  = [(1000._R8*i + real(i, R8)**2, i = 280, 400)]
    m%h_off  = m%h_tab(280) - 280._R8*(m%h_tab(281) - m%h_tab(280))
    allocate(m%e_tab(280:400), source=m%h_tab - m%h_off)
    force = [2.e-3_R8, 0._R8, 0._R8, 0._R8, 0._R8, 1._R8]

    wr = 0._R8; wc = 0._R8; wf = 0._R8; ws = 0._R8
    do k = 1, size(Ts)
      e = mat_e(m, Ts(k))
      ! MK: rho, u, v, w, T, n
      pk = [1.2_R8, 10._R8, -3._R8, 0.5_R8, Ts(k), 1.e9_R8]
      ke = 0.5_R8*pk(1)*(pk(2)**2 + pk(3)**2 + pk(4)**2); un = pk(2)*nrm(1) + pk(3)*nrm(2)
      q(1:6) = prim_2_cons_MK(pk, m)
      b(1:6) = cons_2_prim_MK(q(1:6), m)
      wc = max(wc, abs(q(5) - ke - pk(1)*e)/(pk(1)*e + ke))
      wr = max(wr, maxval(abs(b(1:6) - pk)/abs(pk)))
      f(1:6) = flux_make_MK(pk, nrm, m)
      wf = max(wf, abs(f(5) - f(1)*(ke/pk(1) + e))/abs(f(1)*(ke/pk(1) + e)))
      s(1:6) = source_make_MK(pk, force, m)
      ws = max(ws, abs(s(5) + force(1)*(ke/pk(1) + e + lv))/(force(1)*(ke/pk(1) + e + lv)))
      ! IG: rho, u, v, w, P, T, n
      pi7 = [1.2_R8, 10._R8, -3._R8, 0.5_R8, 40._R8, Ts(k), 1.e9_R8]
      q(1:7) = prim_2_cons_IG(pi7, m)
      b(1:7) = cons_2_prim_IG(q(1:7), m)
      wc = max(wc, abs(q(6) - 0.5_R8*q(5) - pi7(1)*e)/(pi7(1)*e + 0.5_R8*q(5)))
      wr = max(wr, maxval(abs(b(1:7) - pi7)/abs(pi7)))
      f(1:7) = flux_make_IG(pi7, nrm, m)
      wf = max(wf, abs(f(6) - 2.5_R8*pi7(5)*un - f(1)*(ke/pi7(1) + e))/abs(f(1)*(ke/pi7(1) + e)))
      s(1:7) = source_make_IG(pi7, force, m)
      ws = max(ws, abs(s(6) - 0.5_R8*s(5) + force(1)*(e + lv))/(force(1)*(e + lv)))
      ! AG: rho, u, v, w, P11, P12, P13, P22, P23, P33, T, n
      pa = [1.2_R8, 10._R8, -3._R8, 0.5_R8, 40._R8, 3._R8, -2._R8, 30._R8, 1._R8, 20._R8, Ts(k), 1.e9_R8]
      tr = 0.5_R8*(pa(5) + pa(8) + pa(10))
      q(1:12) = prim_2_cons_AG(pa, m)
      b(1:12) = cons_2_prim_AG(q(1:12), m)
      wc = max(wc, abs(q(11) - 0.5_R8*(q(5) + q(8) + q(10)) - pa(1)*e)/(pa(1)*e + ke + tr))
      wr = max(wr, maxval(abs(b(1:12) - pa)/abs(pa)))
      f(1:12) = flux_make_AG(pa, nrm, m)
      wf = max(wf, abs(f(11) - (pa(2)*(pa(5)*nrm(1) + pa(6)*nrm(2)) + pa(3)*(pa(6)*nrm(1) + pa(8)*nrm(2)) &
                                + pa(4)*(pa(7)*nrm(1) + pa(9)*nrm(2))) - f(1)*(ke/pa(1) + tr/pa(1) + e)) &
                   /abs(f(1)*(ke/pa(1) + tr/pa(1) + e)))
      s(1:12) = source_make_AG(pa, force, m)
      ws = max(ws, abs(s(11) - 0.5_R8*(s(5) + s(8) + s(10)) + force(1)*(e + lv))/(force(1)*(e + lv)))
    enddo
    call check(wc <= 1.e-14_R8, 'closures: a tabulated cp stores rho*e(T), MK/IG/AG, below, inside and above the table')
    call check(wr <= 1.e-13_R8, 'closures: the state comes back to round-off, MK/IG/AG, T from 1 mK to 520 K')
    call check(wf <= 1.e-14_R8, 'closures: the energy flux carries e(T), MK/IG/AG')
    call check(ws <= 1.e-14_R8, 'closures: the evaporated mass takes e(T) + lv away, MK/IG/AG')

    ! an empty cell (rho = 1e-20): a constant cp equal to the first segment's slope gives the same T
    c%cs_al = m%e_tab(281) - m%e_tab(280); c%lv_al = lv
    Tv = [150._R8, 1.e-20_R8]
    wr = 0._R8
    do k = 1, 2
      pk = [1.e-20_R8, 1.e-20_R8, 1.e-20_R8, 1.e-20_R8, Tv(k), 1.e-20_R8]
      b(1:6) = cons_2_prim_MK(prim_2_cons_MK(pk, m), m)
      f(1:6) = cons_2_prim_MK(prim_2_cons_MK(pk, c), c)
      wr = max(wr, abs(b(5) - f(5))/f(5))
    enddo
    call check(wr <= 1.e-12_R8, 'closures: an empty cell recovers the T a constant cp gives it')

    ! a constant cp: the INI value and a constant column give the same bits
    t%use_table = .true.; t%T_min = 280; t%T_max = 400; t%lv_al = lv
    allocate(t%rho_tab(280:400), source=1500._R8)
    allocate(t%cs_tab(280:400), source=2000._R8)
    allocate(t%h_tab(280:400)); t%h_tab = [(2000._R8*i, i = 280, 400)]
    allocate(t%e_tab(280:400), source=t%h_tab)
    c%cs_al = 2000._R8
    pk = [1.2_R8, 10._R8, -3._R8, 0.5_R8, 333.3_R8, 1.e9_R8]
    pi7 = [1.2_R8, 10._R8, -3._R8, 0.5_R8, 40._R8, 333.3_R8, 1.e9_R8]
    pa = [1.2_R8, 10._R8, -3._R8, 0.5_R8, 40._R8, 3._R8, -2._R8, 30._R8, 1._R8, 20._R8, 333.3_R8, 1.e9_R8]
    call check(all(prim_2_cons_MK(pk, t) == prim_2_cons_MK(pk, c)) .and. &
               all(cons_2_prim_MK(prim_2_cons_MK(pk, t), t) == cons_2_prim_MK(prim_2_cons_MK(pk, c), c)) .and. &
               all(flux_make_MK(pk, nrm, t) == flux_make_MK(pk, nrm, c)) .and. &
               all(source_make_MK(pk, force, t) == source_make_MK(pk, force, c)) .and. &
               all(prim_2_cons_IG(pi7, t) == prim_2_cons_IG(pi7, c)) .and. &
               all(cons_2_prim_IG(prim_2_cons_IG(pi7, t), t) == cons_2_prim_IG(prim_2_cons_IG(pi7, c), c)) .and. &
               all(flux_make_IG(pi7, nrm, t) == flux_make_IG(pi7, nrm, c)) .and. &
               all(prim_2_cons_AG(pa, t) == prim_2_cons_AG(pa, c)) .and. &
               all(cons_2_prim_AG(prim_2_cons_AG(pa, t), t) == cons_2_prim_AG(prim_2_cons_AG(pa, c), c)) .and. &
               all(flux_make_AG(pa, nrm, t) == flux_make_AG(pa, nrm, c)), &
               'closures: a constant cp gives the same bits from the INI and from a constant column')
    q(1:6) = prim_2_cons_MK(pk, c)
    b(1:6) = cons_2_prim_MK(q(1:6), c)
    ke = 0.5_R8*pk(1)*(pk(2)**2 + pk(3)**2 + pk(4)**2)
    call check(abs(q(5) - (pk(1)*2000._R8*pk(5) + ke)) <= 4._R8*epsilon(1._R8)*q(5) .and. &
               all(abs(b(1:6) - pk) <= 4._R8*epsilon(1._R8)*abs(pk)), &
               'closures: a constant cp stores rho*cp*T and recovers T in one pass')

    call sink_legs(m, t, force)
  end subroutine test_closures


  !> The latent sink at the closure's temperature: a material with an Lv column (2.4e6 - 5000 (T - 300)
  !  J/kg on 280..400 K) and its scalar poisoned gives the sources of the same material without the
  !  column and the scalar set to the column at T, below, on, between and above the nodes, with a
  !  varying (mv) and a constant (mc) cp. Every slot exact but the energy, 1e-15 relative.
  subroutine sink_legs(mv, mc, force)
    use ICE_Lib_MK, only: source_make_MK
    use ICE_Lib_IG, only: source_make_IG
    use ICE_Lib_AG, only: source_make_AG
    type(condensed_phase_t), intent(in) :: mv, mc
    real(R8),                intent(in) :: force(6)
    real(R8), parameter :: Ts(6) = [1.e-3_R8, 280._R8, 333.3_R8, 337.5_R8, 400._R8, 520._R8]
    type(condensed_phase_t) :: tb(2), rf(2)
    real(R8) :: lvt(280:400), pk(6), pi7(7), pa(12), a(12), b(12)
    logical  :: ok(3)
    integer  :: i, j, k

    lvt = [(2.4e6_R8 - 5000._R8*real(i - 300, R8), i = 280, 400)]
    tb(1) = mv; tb(2) = mc; rf(1) = mv; rf(2) = mc
    do j = 1, 2
      tb(j)%use_lv = .true.; tb(j)%lv_al = 1.e30_R8
      allocate(tb(j)%lv_tab(280:400), source=lvt)
    enddo
    ok = .true.
    do j = 1, 2
      do k = 1, size(Ts)
        rf(j)%lv_al = value_at(lvt, 280, 400, Ts(k))
        pk = [1.2_R8, 10._R8, -3._R8, 0.5_R8, Ts(k), 1.e9_R8]
        a(1:6) = source_make_MK(pk, force, tb(j)); b(1:6) = source_make_MK(pk, force, rf(j))
        ok(1) = ok(1) .and. all(a([1,2,3,4,6]) == b([1,2,3,4,6])) .and. abs(a(5) - b(5)) <= 1.e-15_R8*abs(b(5))
        pi7 = [1.2_R8, 10._R8, -3._R8, 0.5_R8, 40._R8, Ts(k), 1.e9_R8]
        a(1:7) = source_make_IG(pi7, force, tb(j)); b(1:7) = source_make_IG(pi7, force, rf(j))
        ok(2) = ok(2) .and. all(a([1,2,3,4,5,7]) == b([1,2,3,4,5,7])) .and. abs(a(6) - b(6)) <= 1.e-15_R8*abs(b(6))
        pa = [1.2_R8, 10._R8, -3._R8, 0.5_R8, 40._R8, 3._R8, -2._R8, 30._R8, 1._R8, 20._R8, Ts(k), 1.e9_R8]
        a(1:12) = source_make_AG(pa, force, tb(j)); b(1:12) = source_make_AG(pa, force, rf(j))
        ok(3) = ok(3) .and. all(a(1:10) == b(1:10)) .and. a(12) == b(12) .and. abs(a(11) - b(11)) <= 1.e-15_R8*abs(b(11))
      enddo
    enddo
    call check(ok(1), 'closures: MK''s latent sink takes the Lv column at T, varying and constant cp, 1 mK to 520 K')
    call check(ok(2), 'closures: IG''s latent sink takes the Lv column at T, varying and constant cp, 1 mK to 520 K')
    call check(ok(3), 'closures: AG''s latent sink takes the Lv column at T, varying and constant cp, 1 mK to 520 K')
  end subroutine sink_legs


  !> A Clausius-Clapeyron curve through one atmosphere at 373.15 K (water-like).
  pure function psat_curve(T) result(p)
    real(R8), intent(in) :: T
    real(R8) :: p
    p = 101325._R8*exp(-4896.8_R8*(1._R8/T - 1._R8/373.15_R8))
  end function psat_curve


  subroutine load(prefix)
    character(len=*), intent(in) :: prefix
    if (allocated(obj_condensed(1)%rho_tab)) deallocate(obj_condensed(1)%rho_tab)
    if (allocated(obj_condensed(1)%cs_tab))  deallocate(obj_condensed(1)%cs_tab)
    if (allocated(obj_condensed(1)%h_tab))   deallocate(obj_condensed(1)%h_tab)
    if (allocated(obj_condensed(1)%e_tab))   deallocate(obj_condensed(1)%e_tab)
    if (allocated(obj_condensed(1)%psat_tab)) deallocate(obj_condensed(1)%psat_tab)
    if (allocated(obj_condensed(1)%lv_tab))  deallocate(obj_condensed(1)%lv_tab)
    if (allocated(obj_condensed(1)%cpv_tab)) deallocate(obj_condensed(1)%cpv_tab)
    obj_condensed(1)%use_table = .false.
    obj_condensed(1)%use_psat  = .false.
    obj_condensed(1)%use_lv    = .false.
    obj_condensed(1)%use_cpv   = .false.
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
