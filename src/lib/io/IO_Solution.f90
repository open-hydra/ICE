module ICE_IO_Solution
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IR_Precision
  use Lib_ORION_data
  use ICE_Global_m,      only: ICE_phase_prefix, ngroups, ncond
  use ICE_Config_Types_m, only: obj_time_scheme
  use ICE_Parameters_m,  only: llen

  implicit none

  integer(kind=I4) :: error
  real(kind=R8)    :: IOtime

  character(llen)  :: condinit, gasinit

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(r_solution_if), pointer, public :: read_bck
  procedure(w_solution_if), pointer, public :: write_bck, write_solution

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
    type(ICE_domain_type), intent(in)           :: grid
    type(orion_data),      intent(inout)        :: IOfield
    character(len=*),      intent(in)           :: file
    character(len=*),      intent(in), optional :: format(2)
  end subroutine w_solution_if
  end interface


contains


  !> Initial setup: wires procedure pointers and resolves initial condition paths.
  subroutine input_solution_setup()
    use ICE_Config_Types_m, only: obj_io, obj_sim_param
    implicit none
    logical      :: present
    integer      :: i
    character(5) :: n
    character(llen) :: try

    select case (trim(obj_io%bck_fmt(1)))
    case ('vtk')
      read_bck  => read_vtk_tec
      write_bck => write_vtk_tec
    case default
      read_bck  => read_vtk_tec
      write_bck => write_vtk_tec
    end select

    if (obj_sim_param%owcoupled) gasinit = 'INPUT/gas.tec'

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
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"rho_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"u_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"v_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"w_p'//trim(str(.true.,p))//'"'
      select case (trim(obj_time_scheme%model(p)))
      case ('MK')
        Onvar = Onvar + 6
      case ('IG')
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"P_p'//trim(str(.true.,p))//'"'
        Onvar = Onvar + 7
      case ('AG')
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"P11_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"P12_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"P13_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"P22_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"P23_p'//trim(str(.true.,p))//'"'
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"P33_p'//trim(str(.true.,p))//'"'
        Onvar = Onvar + 12
      end select
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"T_p'//trim(str(.true.,p))//'"'
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"n_p'//trim(str(.true.,p))//'"'
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
    implicit none
    type(ICE_domain_type), intent(in)           :: grid
    type(orion_data),      intent(inout)        :: IOfield
    character(len=*),      intent(in)           :: file
    character(len=*),      intent(in), optional :: format(2)
    character(len=llen) :: path, localpath_vtk
    integer(kind=I4)    :: b, p, nstart, nend, E_IO

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

    select case (format(1))
    case ('vtk')
      IOfield%vtk%format = format(2)
      localpath_vtk = trim(path)//'vtk/'
      call execute_command_line('mkdir -p '//trim(localpath_vtk))
      E_IO = vtk_write_structured_multiblock(orion=IOfield, &
               vtspath=trim(localpath_vtk)//trim(file), &
               vtmpath=trim(path)//trim(file), varnames=obj_io%Ovarnames, time=grid%time)
    case ('tecplot')
      IOfield%tec%format = format(2)
      E_IO = tec_write_structured_multiblock(orion=IOfield, varnames=obj_io%Ovarnames, &
               filename=trim(path)//trim(file)//'.tec')
    end select

  end subroutine write_vtk_tec


  subroutine read_vtk_tec(IOfield_cond, IOfield_gas)
    use Lib_Tecplot
    use Lib_VTK
    use ICE_Config_Types_m, only: obj_io
    implicit none
    type(orion_data), intent(inout)           :: IOfield_cond
    type(orion_data), intent(inout), optional :: IOfield_gas

    if (present(IOfield_gas)) then
      select case (trim(obj_io%bck_fmt(1)))
      case ('tecplot')
        IOfield_gas%tec%format = 'ascii'
        error = tec_read_structured_multiblock(orion=IOfield_gas, filename=trim(gasinit))
      case ('vtk')
        IOfield_gas%tec%format = obj_io%bck_fmt(2)
        error = vtk_read_structured_multiblock(orion=IOfield_gas, &
                  vtmpath=gasinit(1:len(trim(gasinit))-4), &
                  vtspath='INPUT/vtk/field', time=IOtime)
      end select
    end if

    select case (trim(obj_io%bck_fmt(1)))
    case ('tecplot')
      IOfield_cond%tec%format = 'ascii'
      error = tec_read_structured_multiblock(orion=IOfield_cond, filename=trim(condinit))
    case ('vtk')
      IOfield_cond%tec%format = obj_io%bck_fmt(2)
      error = vtk_read_structured_multiblock(orion=IOfield_cond, &
                vtmpath=condinit(1:len(trim(condinit))-4), &
                vtspath='INPUT/vtk/field', time=IOtime)
    end select

  end subroutine read_vtk_tec


end module ICE_IO_Solution
