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
    use ICE_IO_Solution,   only: io_extension
    use ICE_Lib_Evaporation, only: assign_evaporation, assign_interface, assign_blowing, &
                                   Ru, iMv, iLv, icpv, iLe, iYinf, iLvMvOverRu,     &
                                   iinvTboil, ialphaE
    implicit none
    integer :: p
    logical :: gas_present

    ! --- Parse format strings into arrays ---
    call parse(obj_io%ini_format, ' ', obj_io%ini_fmt)
    call parse(obj_io%sol_format, ' ', obj_io%sol_fmt)

    ! --- Detect one-way coupling from gas file presence ---
    ! The name follows `gas-path` and `ic-format`, as the reader does.
    inquire(file=trim(obj_io%gaspath)//'gas'//io_extension(obj_io%ini_fmt), &
            exist=gas_present)
    obj_sim_param%owcoupled = gas_present

    ! --- Shock detector ---
    obj_space_scheme%SD = (trim(obj_space_scheme%shock_detector) == 'Jameson')

    ! --- Space reconstruction and its limiter ---
    if (trim(obj_space_scheme%space_reconstruction) == 'MUSCL') then
      if (trim(obj_space_scheme%flux_limiter) == 'none') then
        obj_space_scheme%flux_limiter = 'vanleer'
        write(*, '(A)') ' [WARNING] MUSCL without flux-limiter. Van Leer by default.'
      end if
      call assign_limiter(obj_space_scheme%flux_limiter)
    else
      call assign_limiter('IORD')
    end if

    ! --- Assign the global exchange models (only consulted when coupled) ---
    if (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) then
      call assign_drag(obj_time_scheme%drag, obj_time_scheme%dragSelect)
      call assign_heat(obj_time_scheme%heat, obj_time_scheme%heatSelect)
      call assign_evaporation(obj_time_scheme%evaporation, obj_time_scheme%evapSelect)
      call assign_interface(obj_time_scheme%interface_model, obj_time_scheme%intfSelect)
      call assign_blowing(obj_time_scheme%blowing, obj_time_scheme%blowSelect)
    end if

    ! --- Pack the vapour properties into the array Lib_Evaporation indexes ---
    !  The two derived entries are what the saturation pressure is actually built
    !  from, so they are formed once here rather than per cell and per stage.
    obj_condensed%ep = 0._R8
    obj_condensed%ep(iMv)         = obj_condensed%Mv
    obj_condensed%ep(iLv)         = obj_condensed%lv_al
    obj_condensed%ep(icpv)        = obj_condensed%cpv
    obj_condensed%ep(iLe)         = obj_condensed%Le
    obj_condensed%ep(iYinf)       = obj_condensed%Yinf
    obj_condensed%ep(iLvMvOverRu) = obj_condensed%lv_al*obj_condensed%Mv/Ru
    obj_condensed%ep(iinvTboil)   = 1._R8/obj_condensed%Tboil
    obj_condensed%ep(ialphaE)     = obj_condensed%alphaE

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
