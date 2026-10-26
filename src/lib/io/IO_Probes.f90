module ICE_IO_Probes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m
  use ICE_Base_Types_m,     only: ICE_vector_3D_type
  use ICE_Config_Types_m,   only: obj_io_probes, obj_sim_param
  use ICE_Advanced_Types_m, only: ICE_domain_type

  implicit none
  private
  public :: Setup_Probes, Write_Probes_Data, Close_Probes, Restart_Probes, Drop_Probe_Rows, Drop_Rows_After

  type :: real_ptr
    real(R8), pointer :: p => null()
  end type real_ptr

  type :: obj_probe
    character(len=clen)               :: name   !< its section
    character(len=llen)               :: file   !< OUTPUT/<name>.txt
    type(ICE_vector_3D_type)          :: location
    integer, dimension(4)             :: ilocation
    integer                           :: nvar
    character(len=clen), allocatable  :: names(:)
    integer, allocatable              :: family(:), column(:)   !< where each variable is in prim
    type(real_ptr), allocatable       :: variables(:)
    real(R8)                          :: dtime
    integer                           :: ntime = 1   !< next multiple of dtime with a time-based row
    integer                           :: diter
    integer                           :: unit
  contains
    private
    procedure, pass(self) :: Place
  end type obj_probe

  integer, public :: nprobes = 0
  type(obj_probe), allocatable, target :: probe(:)

contains


  subroutine Setup_Probes(domain)
    use ICE_Config_Types_m, only: obj_io_probes, obj_sim_param
    use IR_Precision,       only: str
    use strings,            only: parse
    use ICE_Mod_MPI,        only: mpi_abort_all
#ifdef USE_MPI
    use ICE_Mod_MPI, only: is_local_block
#endif
    implicit none
    type(ICE_domain_type), intent(in) :: domain
    integer             :: i, ii, error
    character(len=clen) :: string(16)

    if (.not. allocated(obj_io_probes)) return
    nprobes = size(obj_io_probes)
    if (nprobes == 0) return
    allocate(probe(1:nprobes))

    do i = 1, nprobes

      ! The key names the probe's section, and the section names the file. The
      ! path is built apart: the key holds clen characters.
      probe(i)%name = obj_io_probes(i)%file
      probe(i)%file = 'OUTPUT/'//trim(probe(i)%name)//'.txt'

      string = ''
      probe(i)%nvar = 0
      call parse(obj_io_probes(i)%varnames, ' ', string)
      do ii = 1, size(string)
        if (trim(string(ii)) /= '') probe(i)%nvar = probe(i)%nvar + 1
      end do
      allocate(probe(i)%names(1:probe(i)%nvar))
      do ii = 1, probe(i)%nvar
        probe(i)%names(ii) = trim(string(ii))
      end do

      probe(i)%dtime = obj_io_probes(i)%dtime
      probe(i)%diter = obj_io_probes(i)%diter

      probe(i)%location%c = obj_io_probes(i)%loc
      probe(i)%ilocation  = obj_io_probes(i)%iloc

      if (sum(obj_io_probes(i)%iloc) == 0) then
        call probe(i)%Place(domain)
      else
        call Check_Index(probe(i), domain)
      end if

      ! Every rank holds the whole domain and the same input, so every rank
      ! takes the names and stops on the same error.
      call Resolve_Variables(probe(i))

#ifdef USE_MPI
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
#endif

      call Assign_Variables(probe(i), domain)

      ! A restart appends to the probe file, creating it if missing (a probe
      ! added, a file moved away); a new run replaces it. Only the rank that
      ! owns the probe gets here, so a failure aborts all ranks.
      if (.not. obj_sim_param%newrun) then
        open(newunit=probe(i)%unit, file=trim(probe(i)%file), &
             status='UNKNOWN', position='APPEND', action='WRITE', iostat=error)
      else
        open(newunit=probe(i)%unit, file=trim(probe(i)%file), &
             status='REPLACE', iostat=error)
      end if
      if (error /= 0) call mpi_abort_all('cannot open probe file '//trim(probe(i)%file))

    end do

  end subroutine Setup_Probes


  !> An index-position names a cell of the mesh: a block of the domain and
  !> i, j, k inside it (ghost cells excluded).
  subroutine Check_Index(probe, domain)
    implicit none
    type(obj_probe),       intent(in) :: probe
    type(ICE_domain_type), intent(in) :: domain
    integer :: b

    b = probe%ilocation(1)
    if (b >= 1 .and. b <= domain%nb) then
      if (all(probe%ilocation(2:4) >= 1 .and. probe%ilocation(2:4) <= domain%blk(b)%dim)) return
    end if
    write(*,'(A,4(1X,I0),A,I0,A)') ' [ERROR] [ICE::Setup_Probes] probe '//trim(probe%name)//': index-position', &
      probe%ilocation, ' is not a cell of the mesh (block 1..', domain%nb, ', i j k inside the block)'
    error stop 1
  end subroutine Check_Index


  !> Each variable name to its family and its column of prim. A family p
  !> answers to the names of the solution files, rho_p<p>, u_p<p>, v_p<p>,
  !> w_p<p>, T_p<p>, n_p<p> (and f_p<p>, chi_p<p> for a solidifying family),
  !> and to the same names without the p (rho_<p>, ...). Any other name stops
  !> the run.
  subroutine Resolve_Variables(probe)
    use IR_Precision, only: str
    use ICE_Global_m, only: ngroups, ncond, nbase, solid_of
    implicit none
    type(obj_probe),  intent(inout) :: probe
    character(len=3), parameter :: base(8) = ['rho', 'u  ', 'v  ', 'w  ', 'T  ', 'n  ', 'f  ', 'chi']
    integer :: v, p, q, col(8)

    allocate(probe%family(probe%nvar), probe%column(probe%nvar))
    probe%family = 0
    probe%column = 0

    do v = 1, probe%nvar
      do p = 1, ngroups
        col = [1, 2, 3, 4, nbase(p)-1, nbase(p), nbase(p)+1, ncond(p)]
        do q = 1, merge(8, 6, solid_of(p))
          if (probe%names(v) == trim(base(q))//'_p'//trim(str(.true.,p)) .or. &
              probe%names(v) == trim(base(q))//'_'//trim(str(.true.,p))) then
            probe%family(v) = p
            probe%column(v) = col(q)
          end if
        end do
      end do
      if (probe%family(v) == 0) then
        write(*,'(A)') ' [ERROR] [ICE::Setup_Probes] probe '//trim(probe%name)//': unknown variable "'// &
          trim(probe%names(v))//'"'
        write(*,'(A,I0)') '   accepted: rho_p<n> u_p<n> v_p<n> w_p<n> T_p<n> n_p<n> (f_p<n> chi_p<n> for a '// &
          'solidifying family), or the same without the p, for a family n = 1..', ngroups
        error stop 1
      end if
    end do

  end subroutine Resolve_Variables


  subroutine Place(self, domain)
    implicit none
    class(obj_probe),                  intent(inout) :: self
    type(ICE_domain_type), intent(in)                :: domain
    integer :: i, j, k, b
    real(R8) :: d0, d

    d0 = huge(1._R8)
    do b = 1, domain%nb
      do k = 0, domain%blk(b)%dim(3)
      do j = 0, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        d = norm2(self%location%c - domain%blk(b)%node(i,j,k)%c)
        if (d < d0) then
          d0 = d
          self%ilocation = [b, i, j, k]
        end if
      end do; end do; end do
    end do

  end subroutine Place


  subroutine Assign_Variables(probe, domain)
    implicit none
    type(obj_probe),       intent(inout), target :: probe
    type(ICE_domain_type), intent(in),    target :: domain
    integer :: v, b, i, j, k

    b = probe%ilocation(1)
    i = probe%ilocation(2)
    j = probe%ilocation(3)
    k = probe%ilocation(4)

    allocate(probe%variables(1:probe%nvar))

    do v = 1, probe%nvar
      probe%variables(v)%p => domain%blk(b)%cond_phase(probe%family(v))%prim(probe%column(v),i,j,k)
    end do

  end subroutine Assign_Variables


  !> The time counter of the time-based rows at a restart. Write_Probes_Data
  !> writes a time-based row when time >= dtime*ntime and then counts it, so the
  !> rows fall on the multiples of dtime, from the first (ntime starts at 1). A
  !> restart that goes on in time (resume) takes for ntime the first multiple
  !> its time has not reached, computed with the same product as
  !> Write_Probes_Data so that an exact multiple counts as written; otherwise,
  !> and in a steady run (time < 0), the rows start again from the first.
  subroutine Restart_Probes(time, resume)
    implicit none
    real(R8), intent(in) :: time
    logical,  intent(in) :: resume
    integer :: i

    do i = 1, nprobes
      probe(i)%ntime = 1
      if (.not. resume .or. time < 0._R8) cycle
      if (time / probe(i)%dtime >= 0.5_R8 * real(huge(1), R8)) cycle
      probe(i)%ntime = max(0, int(time / probe(i)%dtime) - 1)
      do while (time >= probe(i)%dtime * probe(i)%ntime)
        probe(i)%ntime = probe(i)%ntime + 1
      end do
    end do

  end subroutine Restart_Probes


  !> At a restart, the rows of the probe files written after its solution
  !> (Drop_Rows_After: the first value of a row, an iteration or a time, against
  !> last), each file on the rank that owns its probe: closed, cut, and opened
  !> again at its end, as Setup_Probes opens it.
  subroutine Drop_Probe_Rows(last, times)
    use ICE_Mod_MPI, only: mpi_abort_all
#ifdef USE_MPI
    use ICE_Mod_MPI, only: is_local_block
#endif
    implicit none
    real(R8), intent(in) :: last
    logical,  intent(in) :: times
    integer :: i, error

    do i = 1, nprobes
#ifdef USE_MPI
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
#endif
      close(probe(i)%unit)
      call Drop_Rows_After(trim(probe(i)%file), 1 + probe(i)%nvar, 1, last, times)
      open(newunit=probe(i)%unit, file=trim(probe(i)%file), &
           status='UNKNOWN', position='APPEND', action='WRITE', iostat=error)
      if (error /= 0) call mpi_abort_all('cannot open probe file '//trim(probe(i)%file))
    end do

  end subroutine Drop_Probe_Rows


  !> Drop the rows of a history file written after the solution a restart
  !> reads, which the restart writes again. A row holds ncol values and starts
  !> a line; written list-directed it may take more than one line, so the values
  !> are counted to find where each row starts. Value col of a row places it: a
  !> time (written with a decimal point or an exponent) when times, else an
  !> iteration. The first row whose value is of that kind and greater than last
  !> is dropped, with every row after it; a row whose value is of the other kind
  !> is not compared (the iteration rows of a probe in a time-accurate run, whose
  !> iteration count a restart starts again). A last row with fewer than ncol
  !> values (a file cut in the middle of a row) is dropped too, and a line that
  !> starts with '#' is kept and holds no values. The file is truncated and ends
  !> with a newline (End_With_Newline); a missing file is left alone.
  subroutine Drop_Rows_After(file, ncol, col, last, times)
    implicit none
    character(len=*), intent(in) :: file
    integer,          intent(in) :: ncol, col
    real(R8),         intent(in) :: last
    logical,          intent(in) :: times
    character(len=32768) :: line
    integer  :: unit, ios, nkeep, nrow, nval, k1, k2
    real(R8) :: v
    logical  :: later

    open(newunit=unit, file=file, status='OLD', action='READWRITE', iostat=ios)
    if (ios /= 0) return
    nkeep = 0      ! lines kept
    nrow  = 0      ! lines kept before the row being read
    nval  = 0      ! values counted
    later = .false.
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) then
        if (line(verify(line, ' '):verify(line, ' ')) == '#') then
          nkeep = nkeep + 1
          cycle
        end if
        if (mod(nval, ncol) == 0) nrow = nkeep
      end if
      ! The blank-separated values of the line, line(k1:k2)
      k2 = 0
      do
        k1 = verify(line(k2+1:), ' ')
        if (k1 == 0) exit
        k1 = k2 + k1
        k2 = k1 + scan(line(k1:), ' ') - 2
        if (k2 < k1) k2 = len(line)
        nval = nval + 1
        if (mod(nval - 1, ncol) + 1 == col .and. ((scan(line(k1:k2), '.eEdD') > 0) .eqv. times)) then
          read(line(k1:k2), *, iostat=ios) v
          if (ios == 0) later = v > last
          if (later) exit
        end if
        if (k2 >= len(line)) exit
      end do
      if (later) exit
      nkeep = nkeep + 1
    end do
    ! The row later than the solution, or a last row cut short: dropped from its first line
    if (later .or. mod(nval, ncol) /= 0) nkeep = nrow
    rewind(unit)
    do k1 = 1, nkeep
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
    end do
    endfile(unit)
    close(unit)
    call End_With_Newline(file)

  end subroutine Drop_Rows_After


  !> End a file with a newline if its last line has none: the rows a restart
  !> appends would go on that line. Read and written as a stream of bytes.
  subroutine End_With_Newline(file)
    use iso_fortran_env, only: int64
    implicit none
    character(len=*), intent(in) :: file
    integer        :: unit, ios
    integer(int64) :: nbytes
    character      :: last

    open(newunit=unit, file=file, access='STREAM', form='UNFORMATTED', status='OLD', action='READWRITE', iostat=ios)
    if (ios /= 0) return
    inquire(unit=unit, size=nbytes)
    if (nbytes > 0) then
      read(unit, pos=nbytes, iostat=ios) last
      if (ios == 0 .and. last /= new_line('a')) write(unit, pos=nbytes+1, iostat=ios) new_line('a')
    end if
    close(unit)

  end subroutine End_With_Newline


  !> The probe rows due at this iteration, each passed to the file system when
  !> it is written: a run stopped before its end (killed, at the time limit of
  !> a queue) keeps every row it wrote, which the buffer of the unit would lose.
  subroutine Write_Probes_Data(iter, time)
#ifdef USE_MPI
    use ICE_Mod_MPI, only: is_local_block
#endif
    implicit none
    integer, intent(in) :: iter
    real(R8), intent(in) :: time
    integer :: i, v

    do i = 1, nprobes
#ifdef USE_MPI
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
#endif
      if (mod(iter, probe(i)%diter) == 0) then
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p, v=1,probe(i)%nvar)
        flush(probe(i)%unit)
      else if (time >= probe(i)%dtime * probe(i)%ntime) then
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p, v=1,probe(i)%nvar)
        flush(probe(i)%unit)
        probe(i)%ntime = probe(i)%ntime + 1
      end if
    end do

  end subroutine Write_Probes_Data


  !> Close the probe files at the end of the run, each on the rank that owns
  !> its probe, before MPI is finalized: a unit left open is closed by the
  !> run-time library only when the process ends.
  subroutine Close_Probes()
#ifdef USE_MPI
    use ICE_Mod_MPI, only: is_local_block
#endif
    implicit none
    integer :: i

    do i = 1, nprobes
#ifdef USE_MPI
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
#endif
      close(probe(i)%unit)
    end do

  end subroutine Close_Probes


end module ICE_IO_Probes
