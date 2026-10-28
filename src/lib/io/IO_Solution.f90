module ICE_IO_Solution
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IR_Precision
  use Lib_ORION_data
  use ICE_Global_m,      only: ICE_phase_prefix, ngroups, ncond, solid_of
  use ICE_Config_Types_m, only: obj_time_scheme
  use ICE_Parameters_m,  only: llen, hlen
  use ICE_Finite_m,      only: Not_Finite

  implicit none

  integer(kind=I4) :: error = 0, error_gas = 0   !< what the readers return, of the condensed and of the gas field
  real(kind=R8)    :: IOtime

  character(llen)  :: condinit, gasinit

  !> A restart: whether the solution read holds its time, and whether it is a
  !> steady one, whose solution time is then its global iteration count
  !> (restart_time)
  logical :: restart_time_found = .true., restart_steady = .false.

  !> The counters stored with every solution next to CYCLE, the global
  !> iteration count (write_vtk_tec), and read back by a restart
  !> (restart_counters), in this order: see restart_counters_t
  character(len=20), parameter, public :: restart_counter_names(5) = [character(len=20) :: &
    'ICE_ITERATION', 'ICE_OUTPUT_NUMBER', 'ICE_GRID_LEVEL', 'ICE_LAST_OUTPUT_TIME', 'ICE_RES_REFERENCE']

  !> A counter of the restart solution that is not what the writer writes (not
  !> finite, or a count out of range): check_ic stops the run with it
  character(len=llen) :: restart_bad_counter = ''

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(r_solution_if), pointer, public :: read_ic
  procedure(w_solution_if), pointer, public :: write_solution
  public :: io_extension

  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
  subroutine r_solution_if (IOfield_cond, IOfield_gas)
    use Lib_ORION_data
    implicit none
    type(orion_data), intent(inout)           :: IOfield_cond
    type(orion_data), intent(inout), optional :: IOfield_gas
  end subroutine r_solution_if

  subroutine w_solution_if (grid, IOfield, file, format)
    use Lib_ORION_data
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout)        :: grid
    type(orion_data),      intent(inout)        :: IOfield
    character(len=*),      intent(in)           :: file
    character(len=*),      intent(in), optional :: format(2)
  end subroutine w_solution_if
  end interface


contains


  !> File extension for a `<writer> <mode>` format pair. Tecplot binary must be
  !  `.szplt`: TecIO refuses to write a binary file under any other name, and ORION
  !  reads back `.szplt` but not `.plt`, so it is the only binary name that survives
  !  a write-then-restart round trip.
  pure function io_extension(fmt) result(ext)
    implicit none
    character(len=*), intent(in)  :: fmt(2)
    character(len=:), allocatable :: ext

    select case (trim(fmt(1)))
    case ('vtk')
      ext = '.vtm'
    case default
      if (trim(fmt(2)) == 'binary') then
        ext = '.szplt'
      else
        ext = '.tec'
      end if
    end select

  end function io_extension


  !> Initial setup: wires procedure pointers and resolves initial condition paths.
  subroutine input_solution_setup()
    use ICE_Config_Types_m, only: obj_io, obj_sim_param
    implicit none
    logical      :: present
    integer      :: i
    character(llen) :: try

    read_ic => read_vtk_tec

    !> A new run reads INPUT/part-ic.* in ic-format; a restart reads back a solution
    !  ICE wrote itself, so it is sol-format that names it.
    if (.not. obj_sim_param%newrun) obj_io%ini_fmt = obj_io%sol_fmt
    obj_io%extension = io_extension(obj_io%ini_fmt)

    !> The gas field follows `ic-format` like every other input.
    if (obj_sim_param%owcoupled) &
      gasinit = trim(obj_io%gaspath)//'gas'//trim(obj_io%extension)

    if (obj_sim_param%newrun) then
      condinit = 'INPUT/'//trim(ICE_phase_prefix)//'ic'//trim(obj_io%extension)
    else
      condinit = 'OUTPUT/'//trim(ICE_phase_prefix)//'field'//trim(obj_io%extension)
      inquire(file=condinit, exist=present)
      !> Without it (a run stopped before its last write), the last of the
      !  numbered solutions of sol-overwrite = false, named as ICE_postprocess
      !  names them: field1, field2, ...
      if (.not. present) then
        i = 0
        do
          i = i + 1
          try = 'OUTPUT/'//trim(ICE_phase_prefix)//'field'//trim(str(.true.,i))//trim(obj_io%extension)
          inquire(file=try, exist=present)
          if (present) then
            condinit = try
          else
            exit
          end if
        end do
      end if
    end if

  end subroutine input_solution_setup


  !> Output setup
  subroutine output_solution_setup(IOfield)
    use ICE_Config_Types_m, only: obj_io
    implicit none
    type(orion_data), intent(inout) :: IOfield
    integer(kind=I4) :: p, b
    integer(kind=I4) :: Onvar

    !> Condensed phase variables name
    obj_io%Ovarnames = ' '
    Onvar = 0
    do p = 1, ngroups
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "rho_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "u_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "v_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "w_p'//trim(str(.true.,p))//'"'
      select case (trim(obj_time_scheme%model(p)))
      case ('MK')
      case ('IG')
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "P_p'//trim(str(.true.,p))//'"'
      case ('AG')
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "P11_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "P12_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "P13_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "P22_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "P23_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "P33_p'//trim(str(.true.,p))//'"'
      end select
      Onvar = Onvar + ncond(p)
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "T_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "n_p'//trim(str(.true.,p))//'"'
      if (solid_of(p)) then
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "f_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "chi_p'//trim(str(.true.,p))//'"'
      end if
    end do

    !> VTK specification
    IOfield%vtk%format = 'raw'
    IOfield%vtk%node   = .false.

    !> Tecplot specification
    IOfield%tec%node   = .false.
    IOfield%tec%format = 'ascii'
    IOfield%tec%bc     = .false.

    !> Concretize the sol subroutine
    write_solution => write_vtk_tec

    if (Onvar /= Size(IOfield%block(1)%vars, 1)) then
      do b = 1, Size(IOfield%block)
        IOfield%block(b)%name = 'Block'//trim(str(.true.,b))
        deallocate(IOfield%block(b)%vars)
        allocate(IOfield%block(b)%vars(1:Onvar, &
                                       1:IOfield%block(b)%Ni, &
                                       1:IOfield%block(b)%Nj, &
                                       1:IOfield%block(b)%Nk))
      end do
    end if

  end subroutine output_solution_setup


  !> Update the orion-field data and write accordingly to the chosen format
  subroutine write_vtk_tec(grid, IOfield, file, format)
    use IR_Precision, only: str
    use Lib_VTK
    use Lib_Tecplot
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_io, obj_sim_param, obj_multigrid
    use ICE_Mod_MPI,           only: mpi_is_root, mpi_abort_all
    use ICE_Mod_GhostExchange, only: gather_prim_to_root, mpi_io_barrier
    implicit none
    type(ICE_domain_type), intent(inout)        :: grid
    type(orion_data),      intent(inout)        :: IOfield
    character(len=*),      intent(in)           :: file
    character(len=*),      intent(in), optional :: format(2)
    character(len=llen) :: path, localpath_vtk, written
    character(len=hlen) :: what
    character(len=len(obj_io%Ovarnames)) :: vtknames
    integer(kind=I4)    :: b, p, c, nstart, nend, E_IO
    real(R8P)           :: counters(size(restart_counter_names))

    !> Collect the blocks owned by other ranks on root, which writes the file
    call gather_prim_to_root(grid)

    if (mpi_is_root) then

      path = 'OUTPUT/'

      do b = 1, size(IOfield%block)
        IOfield%block(b)%name = 'Block'//trim(str(.true.,b))
        nend = 0
        do p = 1, ngroups
          nstart = nend + 1 ; nend = nend + ncond(p)
          IOfield%block(b)%vars(nstart:nend, &
                                 1:IOfield%block(b)%Ni, &
                                 1:IOfield%block(b)%Nj, &
                                 1:IOfield%block(b)%Nk) = &
            grid%blk(b)%cond_phase(p)%prim(1:ncond(p), &
                                            1:grid%blk(b)%dim(1), &
                                            1:grid%blk(b)%dim(2), &
                                            1:grid%blk(b)%dim(3))
        end do
      end do

      !> The counters that a restart reads back (read_vtk_tec) to go on as the
      !  uninterrupted run would, the same in every format, next to CYCLE, the
      !  global iteration count. grid is the finest level, whose iteration count
      !  stays 0 while a coarse one is solved.
      counters = [real(grid%iter, R8P), real(obj_io%output_number, R8P), real(obj_multigrid%MG_level, R8P), &
                  obj_sim_param%time_from_call, obj_sim_param%res_reference]

      E_IO = 0
      select case (format(1))
      case ('vtk')
        IOfield%vtk%format = format(2)
        !> The names become XML attributes (Name="..."): the quotes of the
        !  Tecplot list go, ORION splits the names on blanks.
        vtknames = obj_io%Ovarnames
        do c = 1, len_trim(vtknames)
          if (vtknames(c:c) == '"') vtknames(c:c) = ' '
        end do
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        written = trim(path)//trim(file)//'.vtm'
        !> TIME is the solution time of the Tecplot files (ICE_postprocess):
        !  minus the global iteration count in a steady run. CYCLE and the
        !  counters are field data of every block file.
        E_IO = vtk_write_structured_multiblock(orion=IOfield, &
                 vtspath=trim(localpath_vtk)//trim(file), &
                 vtmpath=trim(path)//trim(file), varnames=vtknames, time=IOfield%solutiontime, &
                 cycle=int(obj_sim_param%iter_general, I4P), fldnames=restart_counter_names, fldvalues=counters)
      case ('tecplot')
        IOfield%tec%format = format(2)
        written = trim(path)//trim(file)//io_extension(format)
        !> CYCLE and the counters as data set auxiliary data
        E_IO = tec_write_structured_multiblock(orion=IOfield, varnames=obj_io%Ovarnames, &
                 filename=trim(written), cycle=int(obj_sim_param%iter_general, I4P), &
                 fldnames=restart_counter_names, fldvalues=counters)
      end select

      !> A write that the library refuses or that fails returns a non-zero
      !  value, and may leave no file: the run stops on it, with the file. The
      !  other ranks wait in the barrier below; MPI_ABORT ends them.
      if (E_IO /= 0) then
        write(what,'(A,I0,A)') 'cannot write '//trim(written)//' (error ', E_IO, ' of the writer)'
        call mpi_abort_all(trim(what))
      end if

    end if

    call mpi_io_barrier()

  end subroutine write_vtk_tec


  subroutine read_vtk_tec(IOfield_cond, IOfield_gas)
    use Lib_Tecplot
    use Lib_VTK
    use ICE_Config_Types_m, only: obj_io, obj_sim_param
    implicit none
    type(orion_data), intent(inout)           :: IOfield_cond
    type(orion_data), intent(inout), optional :: IOfield_gas
    character(len=llen) :: vtspath
    integer(I4P)        :: ncycle
    real(R8P)           :: counters(size(restart_counter_names))
    logical             :: found(size(restart_counter_names))

    !> A vtk field is a .vtm and the block files it lists, read from the
    !  prefix vtm_block_path finds, in the format they declare (the reader's
    !  own default is base64).
    if (present(IOfield_gas)) then
      select case (trim(obj_io%ini_fmt(1)))
      case ('tecplot')
        IOfield_gas%tec%format = obj_io%ini_fmt(2)
        error_gas = tec_read_structured_multiblock(orion=IOfield_gas, filename=trim(gasinit))
      case ('vtk')
        vtspath = vtm_block_path(gasinit)
        IOfield_gas%vtk%format = vtk_block_format(gasinit, vtspath, IOfield_gas%vtk%format)
        error_gas = vtk_read_structured_multiblock(orion=IOfield_gas, &
                      vtmpath=gasinit(1:len(trim(gasinit))-4), &
                      vtspath=trim(vtspath), time=IOtime)
      end select
    end if

    !> A restart reads the counters stored with its solution too (write_vtk_tec)
    ncycle   = 0
    counters = 0._R8P
    found    = .false.
    select case (trim(obj_io%ini_fmt(1)))
    case ('tecplot')
      IOfield_cond%tec%format = obj_io%ini_fmt(2)
      if (obj_sim_param%newrun) then
        error = tec_read_structured_multiblock(orion=IOfield_cond, filename=trim(condinit))
      else
        error = tec_read_structured_multiblock(orion=IOfield_cond, filename=trim(condinit), &
                  cycle=ncycle, fldnames=restart_counter_names, fldvalues=counters, fldfound=found)
      end if
    case ('vtk')
      vtspath = vtm_block_path(condinit)
      IOfield_cond%vtk%format = vtk_block_format(condinit, vtspath, IOfield_cond%vtk%format)
      if (obj_sim_param%newrun) then
        error = vtk_read_structured_multiblock(orion=IOfield_cond, &
                  vtmpath=condinit(1:len(trim(condinit))-4), &
                  vtspath=trim(vtspath), time=IOtime)
      else
        error = vtk_read_structured_multiblock(orion=IOfield_cond, &
                  vtmpath=condinit(1:len(trim(condinit))-4), &
                  vtspath=trim(vtspath), time=IOtime, &
                  cycle=ncycle, fldnames=restart_counter_names, fldvalues=counters, fldfound=found)
      end if
      !> The time of the files (their TIME, 0 when they hold none) where the
      !  Tecplot reader puts the time of the zone headers
      IOfield_cond%solutiontime = IOtime
    case default
      !> Neither reader matched, so the field would be left EMPTY and every later size()
      !  silently returns zero -- allocate_data then makes a zero-length bc array and the run
      !  segfaults far from here, with nothing in the log pointing back. Note the registered
      !  hits this path. Fail where the cause is still visible.
      write(*,'(A)') '  [ERROR] ICE cannot read the initial condition.'
      write(*,'(A)') '  The format is "'//trim(obj_io%ini_fmt(1))//'", which is neither "tecplot"'
      write(*,'(A)') '  nor "vtk", so no reader was selected. Set it explicitly, e.g.'
      write(*,'(A)') '    ic-format = tecplot ascii'
      error stop 1
    end select

    if (.not. obj_sim_param%newrun) then
      call restart_time(IOfield_cond)
      call restart_counters(ncycle, counters, found)
    end if

  end subroutine read_vtk_tec


  !> The time of the solution that a restart reads, as Wrap_Setup takes it.
  !> ORION's Tecplot readers give the time of the zone headers: a steady
  !> solution holds minus its global iteration count (the sign written since
  !> ORION a978930) or the count with STRANDID = 0 (the files of earlier
  !> ORION), a time-accurate one its time, >= 0. The count of a steady solution
  !> is kept positive. A file whose zone header has no SOLUTIONTIME holds no
  !> time, and check_ic stops the run: ORION reads it as -10 (8cda645) or -huge
  !> (a978930), and -10 is also the time a978930 reads for a steady solution at
  !> iteration 10, so an ASCII file is asked itself (tec_has_time). A vtk
  !> solution gives its TIME, negative for a steady one (0 when the files hold
  !> none, which the reader does not tell). A time that is not finite is left
  !> to check_ic.
  subroutine restart_time(IOfield)
    use ICE_Config_Types_m, only: obj_io
    implicit none
    type(orion_data), intent(inout) :: IOfield

    restart_time_found = .true.
    restart_steady     = .false.
    if (Not_Finite(IOfield%solutiontime)) return

    select case (trim(obj_io%ini_fmt(1)))
    case ('tecplot')
      restart_time_found = IOfield%solutiontime > -0.5_R8 * huge(1._R8)
      if (trim(obj_io%ini_fmt(2)) == 'ascii') restart_time_found = restart_time_found .and. tec_has_time(condinit)
      restart_steady = restart_time_found .and. (IOfield%solutiontime < 0._R8 .or. IOfield%strandid == 0)
    case ('vtk')
      restart_steady = IOfield%solutiontime < 0._R8
    end select
    if (restart_steady) IOfield%solutiontime = abs(IOfield%solutiontime)

  end subroutine restart_time


  !> The counters that a restart reads with its solution (read_vtk_tec), kept
  !> in obj_io%restart when the solution holds them all and each is what the
  !> writer writes: finite, and the first three counts. One that is not stops
  !> the run (restart_bad_counter, check_ic). A solution without them (written
  !> by another program, or by ICE before it wrote them) leaves found false,
  !> and the run starts its counters as a new one (Wrap_Setup).
  subroutine restart_counters(ncycle, counters, found)
    use ICE_Config_Types_m, only: obj_io
    implicit none
    integer(I4P), intent(in) :: ncycle
    real(R8P),    intent(in) :: counters(:)
    logical,      intent(in) :: found(:)
    character(len=32) :: value
    integer :: k

    obj_io%restart%found = .false.
    restart_bad_counter  = ''
    do k = 1, size(counters)
      if (.not. found(k)) cycle
      if (Not_Finite(counters(k))) then
        restart_bad_counter = trim(restart_counter_names(k))//' is not finite (NaN or infinity)'
        return
      end if
      ! ICE_ITERATION, ICE_OUTPUT_NUMBER and ICE_GRID_LEVEL count
      if (k <= 3 .and. (counters(k) < 0._R8P .or. counters(k) >= real(huge(1), R8P))) then
        write(value, '(ES24.16E3)') counters(k)
        restart_bad_counter = trim(restart_counter_names(k))//' = '//trim(adjustl(value))//' is not a count'
        return
      end if
    end do
    if (.not. all(found)) return

    obj_io%restart%found            = .true.
    obj_io%restart%cycle            = ncycle
    obj_io%restart%iteration        = nint(counters(1))
    obj_io%restart%output_number    = nint(counters(2))
    obj_io%restart%grid_level       = nint(counters(3))
    obj_io%restart%last_output_time = counters(4)
    obj_io%restart%res_reference    = counters(5)

  end subroutine restart_counters


  !> The first block file that a .vtm lists (its file attribute without
  !> .vts), '' when the .vtm cannot be read or lists none.
  function vtm_first_entry(vtmfile) result(block)
    implicit none
    character(len=*), intent(in) :: vtmfile
    character(len=llen) :: block
    character(len=llen) :: line
    integer :: unit, ios, i1, i2

    block = ''
    open(newunit=unit, file=vtmfile, status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      i1 = index(line, 'file="')
      if (i1 == 0) cycle
      i2 = index(line, '.vts"')
      if (i2 > i1 + 6) block = line(i1+6:i2-1)
      exit
    end do
    close(unit)

  end function vtm_first_entry


  !> The prefix that the vtk reader puts before the block files a .vtm lists:
  !> the folder of the .vtm, which ORION's writer lists them from, or
  !> INPUT/vtk/field, where earlier versions of ICE took the blocks of a field,
  !> when the first block file is there and not in the folder.
  function vtm_block_path(vtmfile) result(vtspath)
    implicit none
    character(len=*), intent(in) :: vtmfile
    character(len=llen) :: vtspath
    character(len=llen) :: block
    logical :: there

    vtspath = vtmfile(1:index(vtmfile, '/', back=.true.))
    block = vtm_first_entry(vtmfile)
    if (len_trim(block) == 0) return
    inquire(file=trim(vtspath)//trim(block)//'.vts', exist=there)
    if (there) return
    inquire(file='INPUT/vtk/field'//trim(block)//'.vts', exist=there)
    if (there) vtspath = 'INPUT/vtk/field'

  end function vtm_block_path


  !> The vtk format of the block files of a .vtm, 'ascii', 'binary' (base64)
  !> or 'raw' (appended data), as its first block file declares it in its
  !> first DataArray (the blocks are written in one call, in one format), or
  !> fallback when the file cannot tell.
  function vtk_block_format(vtmfile, vtspath, fallback) result(fmt)
    implicit none
    character(len=*), intent(in) :: vtmfile, vtspath, fallback
    character(len=6) :: fmt
    character(len=llen) :: line, block
    integer :: unit, ios

    fmt = fallback
    block = vtm_first_entry(vtmfile)
    if (len_trim(block) == 0) return
    open(newunit=unit, file=trim(vtspath)//trim(block)//'.vts', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, '<DataArray') == 0) cycle
      if (index(line, 'format="appended"') > 0) then
        fmt = 'raw'
      else if (index(line, 'format="binary"') > 0) then
        fmt = 'binary'
      else if (index(line, 'format="ascii"') > 0) then
        fmt = 'ascii'
      end if
      exit
    end do
    close(unit)

  end function vtk_block_format


  !> Whether the first zone header of a Tecplot ASCII file names SOLUTIONTIME:
  !> the lines from the one that starts with ZONE to the first that starts with
  !> a number, in any case. A file that cannot be opened names none.
  logical function tec_has_time(file)
    implicit none
    character(len=*), intent(in) :: file
    character(len=32768) :: line
    integer :: unit, ios, i, c
    logical :: zone

    tec_has_time = .false.
    zone = .false.
    open(newunit=unit, file=file, status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      do i = 1, len_trim(line)
        c = iachar(line(i:i))
        if (c >= iachar('a') .and. c <= iachar('z')) line(i:i) = achar(c - 32)
      end do
      if (zone .and. scan(line(1:1), '0123456789+-.') > 0) exit
      if (line(1:4) == 'ZONE') zone = .true.
      if (zone .and. index(line, 'SOLUTIONTIME') > 0) then
        tec_has_time = .true.
        exit
      end if
    end do
    close(unit)

  end function tec_has_time


  !> An error of the read of the initial condition (or of the solution a
  !> restart reads, or of the gas field it was read with), a restart solution
  !> that holds no time (restart_time), or the time of a restart solution or a
  !> value the condensed field gives the domain that is not finite (NaN,
  !> infinity) stops the run at set-up, with the file and the
  !> place: the node coordinates of every block and the cell variables that
  !> setup_cond reads, the first sum(ncond) of each cell. Every rank reads the
  !> files: the ranks decide together (the largest of a flag), the first rank
  !> that found an error writes it, and every rank stops.
  subroutine check_ic(IOfield_cond, gas_read)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: output_unit
    use ICE_Config_Types_m,    only: obj_sim_param
    use ICE_Mod_MPI,           only: mpi_rank_, mpi_size_, mpi_allreduce_max_r8, mpi_allreduce_min_r8
    use ICE_Mod_GhostExchange, only: mpi_io_barrier
    implicit none
    type(orion_data), intent(in) :: IOfield_cond
    logical,          intent(in) :: gas_read   !< the gas field was read with it (ICE alone, 1-way coupled)
    character(len=1), parameter :: axis(3) = ['x', 'y', 'z']
    character(len=hlen) :: what
    character(len=llen) :: name
    integer(kind=I4)    :: b, i, j, k, d, nvar, off
    real(kind=R8)       :: flag, flag_all, first

    what = ''
    if (error /= 0) then
      write(what,'(A,I0,A)') 'cannot read '//trim(condinit)//' (error ', error, ' of the reader)'
    elseif (gas_read .and. error_gas /= 0) then
      write(what,'(A,I0,A)') 'cannot read '//trim(gasinit)//' (error ', error_gas, ' of the reader)'
    elseif (.not. allocated(IOfield_cond%block)) then
      what = 'cannot read '//trim(condinit)//': it gives no block'
    elseif (.not. obj_sim_param%newrun .and. Not_Finite(IOfield_cond%solutiontime)) then
      what = trim(condinit)//': the solution time is not finite (NaN or infinity)'
    elseif (.not. obj_sim_param%newrun .and. .not. restart_time_found) then
      what = trim(condinit)//': the solution holds no time (its zone header has no SOLUTIONTIME), and a restart '// &
             'goes on from the time of its solution; to start from it, use it as the initial condition (newrun = true)'
    elseif (.not. obj_sim_param%newrun .and. len_trim(restart_bad_counter) > 0) then
      what = trim(condinit)//': '//trim(restart_bad_counter)
    else
      blocks: do b = 1, size(IOfield_cond%block)
        associate (mesh => IOfield_cond%block(b)%mesh)
          do k = lbound(mesh, 4), ubound(mesh, 4)
          do j = lbound(mesh, 3), ubound(mesh, 3)
          do i = lbound(mesh, 2), ubound(mesh, 2)
          do d = 1, size(mesh, 1)
            if (.not. ieee_is_finite(mesh(d,i,j,k))) then
              write(what,'(3A,I0,A,2(I0,A),I0)') 'coordinate ', axis(min(d,3)), ' of block ', b, &
                ', node i, j, k = ', i, ', ', j, ', ', k
              exit blocks
            end if
          end do; end do; end do; end do
        end associate
        if (.not. allocated(IOfield_cond%block(b)%vars)) cycle
        associate (vars => IOfield_cond%block(b)%vars)
          nvar = min(size(vars, 1), sum(ncond))
          do k = lbound(vars, 4), ubound(vars, 4)
          do j = lbound(vars, 3), ubound(vars, 3)
          do i = lbound(vars, 2), ubound(vars, 2)
          do d = 1, nvar
            if (.not. ieee_is_finite(vars(d,i,j,k))) then
              ! The name, when the file gives the names (after the coordinates)
              name = ''
              if (allocated(IOfield_cond%varnames)) then
                off = size(IOfield_cond%varnames) - size(vars, 1)
                if (off >= 0) name = ' ('//trim(adjustl(IOfield_cond%varnames(off+d)))//')'
              end if
              write(what,'(A,I0,2A,I0,A,2(I0,A),I0)') 'variable ', d, trim(name), ' of block ', b, &
                ', cell i, j, k = ', i, ', ', j, ', ', k
              exit blocks
            end if
          end do; end do; end do; end do
        end associate
      end do blocks
      if (len_trim(what) > 0) what = trim(condinit)//' holds a value that is not finite (NaN or infinity): '//trim(what)
    end if

    flag = 0._R8
    if (len_trim(what) > 0) flag = 1._R8
    call mpi_allreduce_max_r8(flag, flag_all)
    if (flag_all > 0._R8) then
      ! Written once, by the first rank that found it, before any rank stops
      call mpi_allreduce_min_r8(merge(real(mpi_rank_, R8), huge(1._R8), flag > 0._R8), first)
      if (flag > 0._R8 .and. nint(first) == mpi_rank_) then
        if (mpi_size_ > 1) then
          write(*,'(A,I0,A)') ' [ERROR] [ICE::read_ic] '//trim(what)//' (rank ', mpi_rank_, ')'
        else
          write(*,'(A)') ' [ERROR] [ICE::read_ic] '//trim(what)
        end if
        flush(output_unit)
      end if
      call mpi_io_barrier()
      error stop 1
    end if

  end subroutine check_ic


end module ICE_IO_Solution
