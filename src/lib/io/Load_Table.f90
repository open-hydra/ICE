!>@brief The condensed-material property table INPUT/<prefix>properties.dat.
!> The header names the columns; the rows sit on consecutive integer kelvins from any
!> Tmin. The header, node and column checks return a code, so the unit test can assert
!> them; Load_Table turns a failed code, or what the row scan finds, into a refusal.
module ICE_Load_Table
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: Load_Table
  public :: classify_table_tokens, check_table_nodes, check_table_columns, table_h_offset, table_reason
  public :: validate_psat_column, psat_reason

  integer, parameter :: ntok_max = 32, tok_len = 64, nzone_max = 64

  !> Codes of the header, node and column checks; table_reason gives the text.
  integer, parameter, public :: TAB_OK = 0, TAB_NO_TEMPERATURE = 1, TAB_NO_CP = 2, TAB_NO_DENSITY = 3, &
                                TAB_NO_ENTHALPY = 4, TAB_TWO_ENTHALPY = 5, TAB_DUPLICATE = 6, &
                                TAB_FEW_ROWS = 7, TAB_OFF_NODE = 8, TAB_NONFINITE = 9, &
                                TAB_RHO_NONPOSITIVE = 10, TAB_CP_NONPOSITIVE = 11, TAB_H_NONMONOTONE = 12, &
                                TAB_H_CP_MISMATCH = 13, TAB_DATUM_MISMATCH = 14, TAB_NEGATIVE_T = 15

  !> Codes of the Psat column checks (IGLOO's order); psat_reason gives the text.
  integer, parameter, public :: PSAT_OK = 0, PSAT_ABSENT = 1, PSAT_NONFINITE = 2, PSAT_NEGATIVE = 3, &
                                PSAT_NONMONOTONE = 4, PSAT_CONSTANT = 5, PSAT_TBOIL_RANGE = 6, &
                                PSAT_TBOIL_MISMATCH = 7

  character(len=*), parameter :: grammar = &
    'expected: VARIABLES = "Temperature", "Cp", "Density", "Enthalpy" (or "Enthalpy_abs")[, "Psat"], '// &
    'Temperature first and the others in any order, one zone, rows on consecutive integer kelvins'

contains

  !> Reads the optional table, one zone per material, and fills the table fields of obj_condensed(:).
  !  Absent => the constant [ICE-Physics] density/specific-heat of a single material (several
  !  materials need the table), and the log says which of the two applies.
  subroutine Load_Table()
    use Lib_ORION_data,     only: orion_data
    use Lib_Tecplot,        only: tec_read_points_multivars
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Global_m,       only: ICE_phase_prefix, nmat
    use ICE_Mod_MPI,        only: mpi_is_root
    use ICE_Lib_Evaporation, only: Patm
    use ICE_Lib_Properties, only: mat_psat
    use ICE_Parameters_m,   only: codename
    implicit none
    type(orion_data)        :: orion
    character(len=512)      :: tablefile
    character(len=tok_len)  :: tokens(ntok_max)
    integer                 :: ntok, icp, irho, ih, ips, code, err, Ni, Tmin, m
    integer                 :: nzone, nsize, ndata, ntrail, badline
    integer                 :: nannounced(nzone_max), nrows(nzone_max)
    logical                 :: exists, relative, found
    character(len=32)       :: expected
    real(R8), allocatable   :: T(:), cp(:), rho(:), h(:)

    tablefile = 'INPUT/'//trim(ICE_phase_prefix)//'properties.dat'

    inquire(file=trim(tablefile), exist=exists)
    if (.not. exists) then
      if (nmat > 1) call refuse(tablefile, 'absent, and '//trim(itoa(nmat))//' materials need it, one zone each')
      if (mpi_is_root) write(*,'(A)') ' [ICE] condensed properties: no '//trim(tablefile)// &
                     ' -- using the constant [ICE-Physics] density/specific-heat'
      return
    endif

    ! Header: the quoted names of the first VARIABLES line
    call read_variables_line(tablefile, tokens, ntok, found)
    if (.not. found) call refuse(tablefile, 'no VARIABLES line')
    call classify_table_tokens(tokens(1:ntok), icp, irho, ih, ips, relative, code)
    if (code /= TAB_OK) call refuse(tablefile, table_reason(code))

    ! Rows: ORION reads them, after a scan that guards what ORION cannot see
    call scan_rows(tablefile, ntok, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
    if (ntrail > 0) call refuse(tablefile, 'text after the last data row (remove it)')
    if (nzone /= nmat .or. nsize /= nmat .or. ndata /= nmat) then
      if (nmat == 1) then
        expected = 'one zone expected'
      else
        expected = trim(itoa(nmat))//' zones expected'
      endif
      call refuse(tablefile, trim(expected)//' (one per material), found '//trim(itoa(nzone))//' ZONE lines, '// &
                  trim(itoa(nsize))//' I= lines and '//trim(itoa(ndata))//' blocks of rows')
    endif
    if (badline > 0) call refuse(tablefile, 'unreadable rows: line '//trim(itoa(badline))//' does not hold '// &
                                 trim(itoa(ntok))//' numbers')
    do m = 1, nmat
      if (nannounced(m) /= nrows(m)) call refuse(tablefile, trim(zone_label(m))//'the zone announces '// &
                                                 trim(itoa(nannounced(m)))//' rows but holds '//trim(itoa(nrows(m))))
    enddo
    err = tec_read_points_multivars(orion, ntok-1, trim(tablefile))
    if (err /= 0) call refuse(tablefile, 'unreadable rows')
    Ni = orion%block(1)%Ni
    do m = 1, nmat
      if (orion%block(m)%Ni /= nrows(m)) call refuse(tablefile, trim(zone_label(m))//'the zone announces '// &
                                                     trim(itoa(orion%block(m)%Ni))//' rows but holds '//trim(itoa(nrows(m))))
      if (orion%block(m)%Ni /= Ni .or. nint(orion%block(m)%mesh(1,1,1,1)) /= nint(orion%block(1)%mesh(1,1,1,1))) &
        call refuse(tablefile, trim(zone_label(m))//'spans '//trim(itoa(orion%block(m)%Ni))//' rows from T = '// &
                    trim(itoa(nint(orion%block(m)%mesh(1,1,1,1))))//', zone 1 '//trim(itoa(Ni))//' from T = '// &
                    trim(itoa(nint(orion%block(1)%mesh(1,1,1,1))))//' (every zone covers the same temperatures)')
    enddo

    allocate(T(Ni), cp(Ni), rho(Ni), h(Ni))
    do m = 1, nmat
      T   = orion%block(m)%mesh(1,1:Ni,1,1)
      cp  = orion%block(m)%vars(icp,1:Ni,1,1)
      rho = orion%block(m)%vars(irho,1:Ni,1,1)
      h   = orion%block(m)%vars(ih,1:Ni,1,1)

      code = check_table_nodes(T)
      if (code /= TAB_OK) call refuse(tablefile, trim(zone_label(m))//table_reason(code))
      code = check_table_columns(T, cp, rho, h, relative)
      if (code /= TAB_OK) call refuse(tablefile, trim(zone_label(m))//table_reason(code))
      call check_ini_value(tablefile, 'density',       obj_condensed(m)%rho_al, rho)
      call check_ini_value(tablefile, 'specific-heat', obj_condensed(m)%cs_al,  cp)

      Tmin = nint(T(1))
      associate (mat => obj_condensed(m))
        mat%T_min = Tmin
        mat%T_max = Tmin + Ni - 1
        allocate(mat%rho_tab(Tmin:Tmin+Ni-1), source=rho)
        allocate(mat%cs_tab(Tmin:Tmin+Ni-1),  source=cp)
        allocate(mat%h_tab(Tmin:Tmin+Ni-1),   source=h)
        mat%rho_varies = any(rho /= rho(1))
        mat%cs_varies  = any(cp /= cp(1))
        mat%h_datum    = merge('relative', 'absolute', relative)
        mat%h_off      = table_h_offset(T, cp, h)
        allocate(mat%e_tab(Tmin:Tmin+Ni-1),   source=h - mat%h_off)
        mat%use_table  = .true.
        mat%description = 'Table-based rho_al(T) and cs_al(T) from '//trim(tablefile)
      end associate

      if (mpi_is_root) then
        if (m == 1) write(*,'(A,I0,A,I0,A,I0,A)') ' [ICE] condensed properties: '//trim(tablefile)//', ', Ni, &
          ' rows on T = ', Tmin, '..', Tmin+Ni-1, ' K (linear between the nodes, end values outside)'
        write(*,'(A,I0,A)') '   material ', m, ': density '//trim(range_text(rho))//', cp '//trim(range_text(cp))// &
          ', enthalpy '//trim(obj_condensed(m)%h_datum)//' (hOff = '//trim(rtoa(obj_condensed(m)%h_off))//' J/kg)'
      endif
      if (ips > 0) call load_psat(orion%block(m)%vars(ips,1:Ni,1,1), m)
    enddo

  contains

    !> The Psat column serves the evaporation models only: a material that does not evaporate does
    !  not use it, a column of zeros is absent, and otherwise it must pass validate_psat_column.
    subroutine load_psat(ps, m)
      real(R8), intent(in) :: ps(:)
      integer,  intent(in) :: m
      integer :: pcode
      if (obj_condensed(m)%evapSelect == 0) then
        if (mpi_is_root) write(*,'(A)') '   Psat column not used: this material does not evaporate'
        return
      endif
      pcode = validate_psat_column(ps, Tmin, Tmin+Ni-1, obj_condensed(m)%Tboil, Patm)
      if (pcode == PSAT_ABSENT) then
        if (mpi_is_root) write(*,'(A)') '   Psat column all zero: Clausius-Clapeyron is used'
        return
      endif
      if (pcode /= PSAT_OK) call refuse(tablefile, trim(zone_label(m))//'Psat column: '//psat_reason(pcode))
      allocate(obj_condensed(m)%psat_tab(Tmin:Tmin+Ni-1), source=ps)
      obj_condensed(m)%use_psat = .true.
      if (mpi_is_root) write(*,'(A)') '   Psat from the table, psat(boiling-temperature) = '// &
        trim(rtoa(mat_psat(obj_condensed(m), obj_condensed(m)%Tboil)/Patm))//' atm'
    end subroutine load_psat


    !> With a table, a set INI value must equal a constant column; against a varying one it is refused.
    subroutine check_ini_value(file, key, ini, col)
      use ICE_Input_Registry, only: param_is_set
      character(len=*), intent(in) :: file, key
      real(R8),         intent(in) :: ini, col(:)
      if (.not. param_is_set(trim(codename)//'-Physics', key)) return
      if (any(col /= col(1))) then
        call refuse(file, trim(zone_label(m))//'['//trim(codename)//'-Physics] '//key// &
                    ' is set but the table''s column varies with T (remove the key: the table owns the property)')
      elseif (abs(ini - col(1)) > 1.e-9_R8*abs(col(1))) then
        call refuse(file, trim(zone_label(m))//'['//trim(codename)//'-Physics] '//key//' = '//trim(rtoa(ini))// &
                    ' differs from the table''s constant '//trim(rtoa(col(1)))//' (remove the key or make them equal)')
      endif
    end subroutine check_ini_value


    !> Refusals name the zone when there are several materials.
    function zone_label(m) result(txt)
      integer, intent(in) :: m
      character(len=96) :: txt
      txt = ''
      if (nmat > 1) txt = 'zone '//trim(itoa(m))//' ('//trim(obj_condensed(m)%name)//'): '
    end function zone_label

  end subroutine Load_Table


  !> Maps the header names to column indices of vars (column 1 is Temperature, read as the mesh).
  pure subroutine classify_table_tokens(tokens, icp, irho, ih, ips, relative, code)
    character(len=*), intent(in)  :: tokens(:)
    integer,          intent(out) :: icp, irho, ih, ips, code
    logical,          intent(out) :: relative
    integer :: t

    icp = 0; irho = 0; ih = 0; ips = 0; relative = .true.
    code = TAB_NO_TEMPERATURE
    if (size(tokens) < 1) return
    if (tokens(1) /= 'Temperature' .and. tokens(1) /= 'temperature') return
    code = TAB_OK
    do t = 2, size(tokens)
      select case (trim(tokens(t)))
      case ('Cp', 'cp')
        if (icp > 0) code = TAB_DUPLICATE
        icp = t - 1
      case ('Density', 'density')
        if (irho > 0) code = TAB_DUPLICATE
        irho = t - 1
      case ('Enthalpy', 'enthalpy', 'Enthalpy_abs', 'enthalpy_abs')
        if (ih > 0) code = merge(TAB_DUPLICATE, TAB_TWO_ENTHALPY, &
                                 relative .eqv. (tokens(t) == 'Enthalpy' .or. tokens(t) == 'enthalpy'))
        ih = t - 1
        relative = (tokens(t) == 'Enthalpy' .or. tokens(t) == 'enthalpy')
      case ('Psat', 'psat', 'PSAT')
        if (ips > 0) code = TAB_DUPLICATE
        ips = t - 1
      end select
      if (code /= TAB_OK) return
    enddo
    if (icp == 0) code = TAB_NO_CP
    if (irho == 0) code = TAB_NO_DENSITY
    if (ih == 0) code = TAB_NO_ENTHALPY
  end subroutine classify_table_tokens


  !> At least two rows, none below 0 K, each within 1e-6 K of the integer node Tmin+i-1.
  pure function check_table_nodes(T) result(code)
    real(R8), intent(in) :: T(:)
    integer :: code, i, Tmin

    code = TAB_FEW_ROWS
    if (size(T) < 2) return
    code = TAB_NONFINITE
    if (.not. all(ieee_is_finite(T))) return
    code = TAB_NEGATIVE_T
    Tmin = nint(T(1))
    if (Tmin < 0) return
    code = TAB_OFF_NODE
    do i = 1, size(T)
      if (abs(T(i) - real(Tmin+i-1, R8)) > 1.e-6_R8) return
    enddo
    code = TAB_OK
  end function check_table_nodes


  !> Finite, positive density and cp, strictly increasing enthalpy consistent with cp (a
  !  varying cp: the trapezoid within 1e-3 of the step or 1e-6 of |h|, which absorbs the seams
  !  of a polynomial fit), and a relative enthalpy without offset when cp is constant.
  pure function check_table_columns(T, cp, rho, h, relative) result(code)
    real(R8), intent(in) :: T(:), cp(:), rho(:), h(:)
    logical,  intent(in) :: relative
    integer  :: code, i, n
    real(R8) :: hOff, dh

    n = size(T)
    code = TAB_NONFINITE
    if (.not. (all(ieee_is_finite(cp)) .and. all(ieee_is_finite(rho)) .and. all(ieee_is_finite(h)))) return
    code = TAB_RHO_NONPOSITIVE
    if (any(rho <= 0._R8)) return
    code = TAB_CP_NONPOSITIVE
    if (any(cp <= 0._R8)) return
    code = TAB_H_NONMONOTONE
    if (any(h(2:n) <= h(1:n-1))) return
    code = TAB_H_CP_MISMATCH
    if (all(cp == cp(1))) then
      hOff = table_h_offset(T, cp, h)
      do i = 1, n
        if (abs(h(i) - cp(1)*T(i) - hOff) > 1.e-6_R8*max(abs(h(i)), cp(1)*T(i))) return
      enddo
      code = TAB_DATUM_MISMATCH
      if (relative .and. abs(hOff) > 1.e-6_R8*cp(1)*T(1)) return
    else
      do i = 1, n-1
        dh = h(i+1) - h(i)
        if (abs(dh - 0.5_R8*(cp(i) + cp(i+1))*(T(i+1) - T(i))) > &
            max(1.e-3_R8*dh, 1.e-6_R8*max(abs(h(i)), abs(h(i+1))))) return
      enddo
    endif
    code = TAB_OK
  end function check_table_columns


  !> Enthalpy at 0 K along the first segment: h(Tmin) - cp*Tmin for a constant cp.
  pure function table_h_offset(T, cp, h) result(hOff)
    real(R8), intent(in) :: T(:), cp(:), h(:)
    real(R8) :: hOff
    if (all(cp == cp(1))) then
      hOff = h(1) - cp(1)*T(1)
    else
      hOff = h(1) - T(1)*(h(2) - h(1))
    endif
  end function table_h_offset


  !> Finite; all zero means absent; non-negative; non-decreasing; not constant; the boiling
  !  temperature inside [Tmin, Tmax-1] and psat there within a factor 2 of one atmosphere.
  pure function validate_psat_column(tab, lo, hi, Tboil, Patm) result(code)
    integer,  intent(in) :: lo, hi
    real(R8), intent(in) :: tab(lo:hi), Tboil, Patm
    integer  :: code, i
    real(R8) :: ratio
    code = PSAT_NONFINITE
    if (.not. all(ieee_is_finite(tab))) return
    code = PSAT_ABSENT
    if (all(tab == 0._R8)) return
    code = PSAT_NEGATIVE
    if (any(tab < 0._R8)) return
    code = PSAT_NONMONOTONE
    if (any(tab(lo+1:hi) < tab(lo:hi-1))) return
    code = PSAT_CONSTANT
    if (all(tab == tab(lo))) return
    code = PSAT_TBOIL_RANGE
    if (.not. (Tboil >= real(lo, R8) .and. Tboil <= real(hi-1, R8))) return
    i = int(Tboil)
    ratio = (tab(i) + (tab(i+1) - tab(i))*(Tboil - real(i, R8)))/Patm
    code = PSAT_TBOIL_MISMATCH
    if (ratio < 0.5_R8 .or. ratio > 2._R8) return
    code = PSAT_OK
  end function validate_psat_column


  pure function psat_reason(code) result(txt)
    integer, intent(in) :: code
    character(len=:), allocatable :: txt
    select case (code)
    case (PSAT_NONFINITE);      txt = 'a value that is not finite'
    case (PSAT_NEGATIVE);       txt = 'a negative pressure'
    case (PSAT_NONMONOTONE);    txt = 'a pressure that decreases with T'
    case (PSAT_CONSTANT);       txt = 'a constant pressure'
    case (PSAT_TBOIL_RANGE);    txt = 'boiling-temperature outside [Tmin, Tmax-1] of the table'
    case (PSAT_TBOIL_MISMATCH); txt = 'psat(boiling-temperature) is not within a factor 2 of one atmosphere'
    case default;               txt = 'unknown Psat error'
    end select
  end function psat_reason


  pure function table_reason(code) result(txt)
    integer, intent(in) :: code
    character(len=:), allocatable :: txt
    select case (code)
    case (TAB_NO_TEMPERATURE); txt = 'the first column is not "Temperature"'
    case (TAB_NO_CP);          txt = 'no "Cp" column'
    case (TAB_NO_DENSITY);     txt = 'no "Density" column'
    case (TAB_NO_ENTHALPY);    txt = 'no "Enthalpy" or "Enthalpy_abs" column'
    case (TAB_TWO_ENTHALPY);   txt = 'both "Enthalpy" and "Enthalpy_abs" (the enthalpy datum is ambiguous)'
    case (TAB_DUPLICATE);      txt = 'a column named twice'
    case (TAB_FEW_ROWS);       txt = 'fewer than two rows'
    case (TAB_OFF_NODE);       txt = 'rows not on consecutive integer kelvins'
    case (TAB_NEGATIVE_T);     txt = 'a temperature below 0 K'
    case (TAB_NONFINITE);      txt = 'a value that is not finite'
    case (TAB_RHO_NONPOSITIVE); txt = 'a density that is not positive'
    case (TAB_CP_NONPOSITIVE); txt = 'a cp that is not positive'
    case (TAB_H_NONMONOTONE);  txt = 'an enthalpy that does not increase with T'
    case (TAB_H_CP_MISMATCH);  txt = 'an enthalpy that disagrees with cp'
    case (TAB_DATUM_MISMATCH); txt = 'a relative "Enthalpy" column with an offset (h(Tmin) /= cp*Tmin: tag it "Enthalpy_abs")'
    case default;              txt = 'unknown table error'
    end select
  end function table_reason


  !> The quoted names of the first line that contains VARIABLES.
  subroutine read_variables_line(file, tokens, ntok, found)
    character(len=*),       intent(in)  :: file
    character(len=tok_len), intent(out) :: tokens(ntok_max)
    integer,                intent(out) :: ntok
    logical,                intent(out) :: found
    character(len=4096) :: line
    integer :: unit, ios, i, j, k

    ntok = 0; found = .false.; tokens = ''
    open(newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, 'VARIABLES') == 0) cycle
      found = .true.
      k = 1
      do while (ntok < ntok_max .and. k < len(line))
        i = index(line(k:), '"')
        if (i == 0) exit
        i = k + i - 1
        j = index(line(i+1:), '"')
        if (j == 0) exit
        j = i + j
        ntok = ntok + 1
        tokens(ntok) = line(i+1:j-1)
        k = j + 1
      enddo
      exit
    enddo
    close(unit)
  end subroutine read_variables_line


  !> What ORION's point reader does not check, line by line: the ZONE lines (its zone count), the
  !  lines holding I= (it gives each one a block), the rows each of them announces and each block
  !  of rows holds, every data row holding ntok numbers (ORION keeps only the last row's status),
  !  and text after the last row.
  subroutine scan_rows(file, ntok, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
    character(len=*), intent(in)  :: file
    integer,          intent(in)  :: ntok
    integer,          intent(out) :: nzone, nsize, ndata, nannounced(nzone_max), nrows(nzone_max), ntrail, badline
    character(len=4096) :: line
    real(R8) :: x, v(ntok_max)
    integer  :: unit, ios, rd, iline
    logical  :: numeric, prev_numeric

    nzone = 0; nsize = 0; ndata = 0; nannounced = -1; nrows = 0; ntrail = 0; badline = 0
    prev_numeric = .false.; iline = 0
    open(newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      iline = iline + 1
      if ((index(line, 'ZONE') > 0 .and. index(line, 'ZONETYPE') == 0) .or. index(line, 'Zone') > 0 .or. &
          index(line, 'ZONE T') > 0) nzone = nzone + 1     ! ORION's own rule
      if (index(line, 'I=') > 0) then
        nsize = nsize + 1
        if (nsize <= nzone_max) nannounced(nsize) = announced_rows(line)
      endif
      read(line, *, iostat=rd) x
      numeric = (rd == 0 .and. index(line, 'DATA') == 0)
      if (numeric) then
        if (.not. prev_numeric) ndata = ndata + 1
        if (ndata <= nzone_max) nrows(ndata) = nrows(ndata) + 1
        ntrail = 0
        if (badline == 0) then
          read(line, *, iostat=rd) v(1:ntok)
          if (rd /= 0 .or. index(line, '/') > 0) badline = iline
        endif
      else if (ndata > 0) then
        ntrail = ntrail + 1
      endif
      prev_numeric = numeric
    enddo
    close(unit)
  end subroutine scan_rows


  !> The row count of a size line as ORION takes it: I= in one of its first two comma-separated fields.
  pure function announced_rows(line) result(n)
    character(len=*), intent(in) :: line
    integer :: n, f, c, k, ios
    character(len=len(line)) :: rest, field
    n = -1
    rest = line
    do f = 1, 2
      c = index(rest, ',')
      if (c > 0) then
        field = rest(:c-1)
        rest  = rest(c+1:)
      else
        field = rest
        rest  = ''
      endif
      k = index(field, 'I=')
      if (k > 0) then
        read(field(k+2:), *, iostat=ios) n
        if (ios /= 0) n = -1
        return
      endif
    enddo
  end function announced_rows


  !> Every rank prints, as the BC reader does: a rank that stops first cannot swallow the reason.
  subroutine refuse(file, reason)
    character(len=*), intent(in) :: file, reason
    write(*,'(A)') ' [ERROR] [ICE::Load_Table] '//trim(file)//': '//reason
    write(*,'(A)') '         '//grammar
    error stop 1
  end subroutine refuse


  pure function range_text(col) result(txt)
    real(R8), intent(in) :: col(:)
    character(len=:), allocatable :: txt
    if (all(col == col(1))) then
      txt = trim(rtoa(col(1)))//' (constant)'
    else
      txt = trim(rtoa(minval(col)))//'..'//trim(rtoa(maxval(col)))//' (varies)'
    endif
  end function range_text


  pure function rtoa(x) result(txt)
    real(R8), intent(in) :: x
    character(len=24) :: txt
    write(txt, '(ES12.5)') x
    txt = adjustl(txt)
  end function rtoa


  pure function itoa(i) result(txt)
    integer, intent(in) :: i
    character(len=16) :: txt
    write(txt, '(I0)') i
  end function itoa

end module ICE_Load_Table
