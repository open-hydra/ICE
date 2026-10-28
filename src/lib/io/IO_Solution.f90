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
    character(5) :: n
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
      if (.not. present) then
        i = 0
        do
          i = i + 1
          write(n, '(I5.5)') i
          try = 'OUTPUT/'//trim(ICE_phase_prefix)//'field'//trim(adjustl(n))//trim(obj_io%extension)
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
    use ICE_Config_Types_m, only: obj_io
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
        E_IO = vtk_write_structured_multiblock(orion=IOfield, &
                 vtspath=trim(localpath_vtk)//trim(file), &
                 vtmpath=trim(path)//trim(file), varnames=vtknames, time=grid%time)
      case ('tecplot')
        IOfield%tec%format = format(2)
        written = trim(path)//trim(file)//io_extension(format)
        E_IO = tec_write_structured_multiblock(orion=IOfield, varnames=obj_io%Ovarnames, &
                 filename=trim(written))
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

    if (present(IOfield_gas)) then
      select case (trim(obj_io%ini_fmt(1)))
      case ('tecplot')
        IOfield_gas%tec%format = obj_io%ini_fmt(2)
        error_gas = tec_read_structured_multiblock(orion=IOfield_gas, filename=trim(gasinit))
      case ('vtk')
        IOfield_gas%tec%format = obj_io%ini_fmt(2)
        error_gas = vtk_read_structured_multiblock(orion=IOfield_gas, &
                      vtmpath=gasinit(1:len(trim(gasinit))-4), &
                      vtspath='INPUT/vtk/field', time=IOtime)
      end select
    end if

    select case (trim(obj_io%ini_fmt(1)))
    case ('tecplot')
      IOfield_cond%tec%format = obj_io%ini_fmt(2)
      error = tec_read_structured_multiblock(orion=IOfield_cond, filename=trim(condinit))
    case ('vtk')
      IOfield_cond%tec%format = obj_io%ini_fmt(2)
      error = vtk_read_structured_multiblock(orion=IOfield_cond, &
                vtmpath=condinit(1:len(trim(condinit))-4), &
                vtspath='INPUT/vtk/field', time=IOtime)
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

    if (.not. obj_sim_param%newrun) call restart_time(IOfield_cond)

  end subroutine read_vtk_tec


  !> The time of the solution that a restart reads, as Wrap_Setup takes it.
  !> ORION's Tecplot readers give the time of the zone headers: a steady
  !> solution holds minus its global iteration count (the sign written since
  !> ORION a978930) or the count with STRANDID = 0 (the files of earlier
  !> ORION), a time-accurate one its time, >= 0. The count of a steady solution
  !> is kept positive. A file whose zone header has no SOLUTIONTIME holds no
  !> time, and check_ic stops the run: ORION reads it as -10 (8cda645) or -huge
  !> (a978930), and -10 is also the time a978930 reads for a steady solution at
  !> iteration 10, so an ASCII file is asked itself (tec_has_time). A time that
  !> is not finite is left to check_ic.
  subroutine restart_time(IOfield)
    use ICE_Config_Types_m, only: obj_io
    implicit none
    type(orion_data), intent(inout) :: IOfield

    restart_time_found = .true.
    restart_steady     = .false.
    if (Not_Finite(IOfield%solutiontime)) return
    if (trim(obj_io%ini_fmt(1)) /= 'tecplot') return

    restart_time_found = IOfield%solutiontime > -0.5_R8 * huge(1._R8)
    if (trim(obj_io%ini_fmt(2)) == 'ascii') restart_time_found = restart_time_found .and. tec_has_time(condinit)
    restart_steady = restart_time_found .and. (IOfield%solutiontime < 0._R8 .or. IOfield%strandid == 0)
    if (restart_steady) IOfield%solutiontime = abs(IOfield%solutiontime)

  end subroutine restart_time


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
