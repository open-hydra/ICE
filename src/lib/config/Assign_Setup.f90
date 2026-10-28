module ICE_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private
  public :: Assign_Setup, Check_Riemann_Solver

contains

  subroutine Assign_Setup()
    use ICE_Config_Types_m
    use ICE_Global_m,      only: ngroups, nrk, npop, ncond, nbase, solid_of, ICE_phase_prefix
    use strings,           only: parse
    use ICE_Lib_Limiters,  only: assign_limiter
    use ICE_Lib_Drag,      only: assign_drag
    use ICE_Lib_Heat,      only: assign_heat
    use ICE_IO_Solution,   only: io_extension
    use ICE_Lib_Evaporation, only: assign_evaporation, assign_interface, assign_blowing
    use ICE_Mod_MPI,       only: mpi_is_root
    implicit none
    integer :: p
    logical :: gas_present

    ! --- Condensed phase to read: [ICE-Parameters] phase = <ATLAS phase name> ---
    if (len_trim(adjustl(obj_sim_param%phase)) > 0) then
      obj_sim_param%phase = adjustl(obj_sim_param%phase)
      if (index(trim(obj_sim_param%phase), '-') > 0) then
        write(*,'(A)') ' [ERROR] [ICE-Parameters] phase must not contain "-" (it is the ATLAS phase name)'
        error stop 1
      endif
      ICE_phase_prefix = trim(obj_sim_param%phase)//'-'
    endif

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
        if (mpi_is_root) write(*, '(A)') ' [WARNING] MUSCL without flux-limiter. Van Leer by default.'
      end if
      call assign_limiter(obj_space_scheme%flux_limiter)
    else
      if (trim(obj_space_scheme%flux_limiter) /= 'none' .and. mpi_is_root) &
        write(*, '(A)') ' [WARNING] flux-limiter = '//trim(obj_space_scheme%flux_limiter)// &
                        ' ignored: space-reconstruction is first-order.'
      call assign_limiter('IORD')
    end if

    ! --- Riemann solver: Saurel for MK families only ---
    call Check_Riemann_Solver()

    ! --- Implicit residual smoothing: steady runs only ---
    ! The smoothed residual is not the residual of the space discretisation, so a
    ! time-accurate run with it does not follow the transient of the scheme.
    if (obj_irs%enabled .and. obj_time_scheme%time_accurate) then
      write(*,'(A)') ' [ERROR] [ICE-Numerics] irs = true needs time-accurate = false: implicit residual '// &
                     'smoothing accelerates steady runs and changes the transient of a time-accurate one'
      error stop 1
    end if

    ! --- Assign the global exchange models (only consulted when coupled) ---
    if (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) then
      call assign_drag(obj_time_scheme%drag, obj_time_scheme%dragSelect)
      call assign_heat(obj_time_scheme%heat, obj_time_scheme%heatSelect)
      call assign_evaporation(obj_time_scheme%evaporation, obj_time_scheme%evapSelect)
      call assign_interface(obj_time_scheme%interface_model, obj_time_scheme%intfSelect)
      call assign_blowing(obj_time_scheme%blowing, obj_time_scheme%blowSelect)
    end if

    ! --- Validate models and compute ncond/npop ---
    allocate(npop(1:3)); npop = 0
    allocate(ncond(1:ngroups)); ncond = 0
    allocate(nbase(1:ngroups)); nbase = 0
    allocate(solid_of(1:ngroups)); solid_of = .false.
    do p = 1, ngroups
      select case (trim(obj_time_scheme%model(p)))
      case ('MK')
        npop(1)  = npop(1) + 1
        nbase(p) = 6
      case ('IG')
        npop(2)  = npop(2) + 1
        nbase(p) = 7
      case ('AG')
        npop(3)  = npop(3) + 1
        nbase(p) = 12
      end select
    end do
    ncond = nbase

    ! --- nrk global (re-read at first Explicit_Step call) ---
    nrk = 1

  end subroutine Assign_Setup


  !> Saurel upwinds on the mean normal velocity and assembles the flux from the MK variable
  !> layout, so it is meaningless for the Gaussian closures. riemann-solver is resolved at
  !> every step (Mod_Explicit), so the choice is checked at set-up and after each re-read
  !> of input.ini. With err present (a re-read with ini-error = warn) err returns the error,
  !> blank when there is none, instead of the stop.
  subroutine Check_Riemann_Solver(err)
    use ICE_Config_Types_m, only: obj_time_scheme
    use ICE_Global_m,       only: ngroups
    implicit none
    character(len=*), intent(out), optional :: err
    character(len=256) :: msg
    integer :: p

    if (present(err)) err = ''
    if (trim(obj_time_scheme%riemann) /= 'Saurel') return
    do p = 1, ngroups
      if (trim(obj_time_scheme%model(p)) /= 'MK') then
        write(msg,'(A,I0,A)') ' [ERROR] [ICE-Numerics] riemann-solver = Saurel is specific to the MK closure, '// &
                              'and family ', p, ' is '//trim(obj_time_scheme%model(p))//'; use Rusanov or HLLE '// &
                              'with IG and AG.'
        if (present(err)) then
          err = msg
          return
        end if
        write(*,'(A)') trim(msg)
        error stop 1
      end if
    end do

  end subroutine Check_Riemann_Solver

end module ICE_Assign_Setup
