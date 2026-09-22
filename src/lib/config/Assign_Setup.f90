module ICE_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private
  public :: Assign_Setup

contains

  subroutine Assign_Setup()
    use ICE_Config_Types_m
    use ICE_Global_m,      only: ngroups, nrk, npop, ncond
    use strings,           only: parse
    use ICE_Lib_Limiters,  only: assign_limiter
    use ICE_Lib_Drag,      only: assign_drag
    use ICE_Lib_Heat,      only: assign_heat
    implicit none
    integer :: p
    logical :: gas_present

    ! --- Detect one-way coupling from gas file presence ---
    inquire(file='INPUT/gas.tec', exist=gas_present)
    obj_sim_param%owcoupled = gas_present

    ! --- Parse format strings into arrays ---
    call parse(obj_io%sol_format, ' ', obj_io%sol_fmt)
    call parse(obj_io%bck_format, ' ', obj_io%bck_fmt)

    ! --- Set file extension from backup format ---
    select case (trim(obj_io%bck_fmt(1)))
    case ('vtk')
      obj_io%extension = '.vtm'
    case default
      obj_io%extension = '.tec'
    end select

    ! --- Derive SD flag and handle missing limiter ---
    obj_space_scheme%SD = (index(obj_space_scheme%space_reconstruction, 'SD') > 0)
    if (index(obj_space_scheme%space_reconstruction, 'MUSCL') > 0) then
      if (trim(obj_space_scheme%flux_limiter) == 'none') then
        obj_space_scheme%flux_limiter = 'VANLEER'
        write(*, '(A)') ' [WARNING] MUSCL without flux-limiter. Van Leer by default.'
      end if
      call assign_limiter(obj_space_scheme%flux_limiter)
    end if

    ! --- Assign global drag/heat (only used when coupled) ---
    if (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) then
      call assign_drag(obj_time_scheme%drag, obj_time_scheme%dragSelect)
      call assign_heat(obj_time_scheme%heat)
    end if

    ! --- Validate models and compute ncond/npop ---
    allocate(npop(1:3)); npop = 0
    allocate(ncond(1:ngroups)); ncond = 0
    do p = 1, ngroups
      select case (trim(obj_time_scheme%model(p)))
      case ('MK')
        npop(1)  = npop(1) + 1
        ncond(p) = 6
      case ('IG')
        npop(2)  = npop(2) + 1
        ncond(p) = 7
      case ('AG')
        npop(3)  = npop(3) + 1
        ncond(p) = 12
      end select
    end do

    ! --- nrk global (re-read at first Explicit_Step call) ---
    nrk = 1

  end subroutine Assign_Setup

end module ICE_Assign_Setup
