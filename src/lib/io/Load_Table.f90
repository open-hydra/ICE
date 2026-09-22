module ICE_Load_Table
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_Table, get_rho_al, get_cs_al

contains

  !> Optional table of rho(T) and cs(T) for the condensed phase, read from
  !  INPUT/<phase-prefix>properties.dat. Absent => get_rho_al/get_cs_al fall back to the constant
  !  [ICE-Physics] rho/cs. Both outcomes are REPORTED: the fallback is a legitimate configuration,
  !  but a silent one is indistinguishable from "the table was there and I failed to find it", which
  !  is exactly how a mis-named file turns into a plausible wrong answer.
  !
  !  The path MUST be built from ICE_phase_prefix, like IO_BC.f90 and IO_Solution.f90. It used to be
  !  the literal 'INPUT/part-properties.dat', which agreed with the prefix only because the default
  !  is 'part-' (Global_m.f90) -- so any case whose phase file was not part-phase.txt silently lost
  !  its table. Callers: Wrap_Setup.f90 calls this unconditionally, AFTER the prefix is assigned.
  subroutine Load_Table()
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Global_m,       only: ICE_phase_prefix
    use ICE_Mod_MPI,        only: mpi_is_root
    implicit none
    integer  :: ios, unitFile, npts, i, iT
    real(R8) :: T, cp, rho, h
    character(len=512) :: line, token
    character(len=512) :: tablefile

    tablefile = 'INPUT/'//trim(ICE_phase_prefix)//'properties.dat'

    open(newunit=unitFile, file=trim(tablefile), status='old', iostat=ios)
    if (ios /= 0) then
      if (mpi_is_root) write(*,'(A)') ' [ICE] condensed properties: no '//trim(tablefile)// &
                     ' -- using the constant [ICE-Physics] rho/cs'
      return
    endif

    ! Skip TITLE line
    read(unitFile, '(A)', iostat=ios) line
    ! Skip VARIABLES line
    read(unitFile, '(A)', iostat=ios) line
    ! Skip ZONE line
    read(unitFile, '(A)', iostat=ios) line
    ! Read "I=N, F=POINT" line and extract N
    read(unitFile, '(A)', iostat=ios) line
    call extract_integer_after(line, 'I=', npts)
    if (npts <= 0) then
      close(unitFile)
      if (mpi_is_root) write(*,'(A)') ' [ICE] condensed properties: '//trim(tablefile)// &
                     ' has no readable "I=<n>" point count -- using the constant [ICE-Physics] rho/cs'
      return
    end if

    obj_condensed%T_min = 1
    obj_condensed%T_max = npts
    allocate(obj_condensed%rho_tab(1:npts))
    allocate(obj_condensed%cs_tab(1:npts))

    do i = 1, npts
      read(unitFile, *, iostat=ios) T, cp, rho, h
      if (ios /= 0) exit
      iT = nint(T)
      if (iT >= 1 .and. iT <= npts) then
        obj_condensed%cs_tab(iT)  = cp
        obj_condensed%rho_tab(iT) = rho
      end if
    end do
    close(unitFile)

    obj_condensed%use_table   = .true.
    obj_condensed%description = 'Table-based rho_al(T) and cs_al(T) from '//trim(tablefile)

    if (mpi_is_root) write(*,'(A,I0,A)') ' [ICE] condensed properties: table from '//trim(tablefile)//' (', npts, &
                        ' points, indexed by integer T)'

  end subroutine Load_Table


  function get_rho_al(Tp) result(rho)
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    real(R8), intent(in) :: Tp
    real(R8)             :: rho
    integer              :: iT
    if (obj_condensed%use_table) then
      iT  = min(max(nint(Tp), obj_condensed%T_min), obj_condensed%T_max)
      rho = obj_condensed%rho_tab(iT)
    else
      rho = obj_condensed%rho_al
    end if
  end function get_rho_al


  function get_cs_al(Tp) result(cs)
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    real(R8), intent(in) :: Tp
    real(R8)             :: cs
    integer              :: iT
    if (obj_condensed%use_table) then
      iT = min(max(nint(Tp), obj_condensed%T_min), obj_condensed%T_max)
      cs = obj_condensed%cs_tab(iT)
    else
      cs = obj_condensed%cs_al
    end if
  end function get_cs_al


  subroutine extract_integer_after(line, key, val)
    implicit none
    character(len=*), intent(in)  :: line, key
    integer,          intent(out) :: val
    integer :: pos, ios
    val = 0
    pos = index(line, key)
    if (pos == 0) return
    read(line(pos+len(key):), *, iostat=ios) val
  end subroutine extract_integer_after

end module ICE_Load_Table
