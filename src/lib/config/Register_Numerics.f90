module ICE_Read_Numerics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Config_Types_m, only: obj_time_scheme, obj_space_scheme, obj_multigrid
  use ICE_Input_Registry
  implicit none
  private
  public :: Register_Numerics, Register_Families

contains

  subroutine Register_Numerics(nmgl)
    use ICE_Parameters_m, only: codename
    use ICE_Config_Types_m, only: obj_irs
    implicit none
    integer, intent(in) :: nmgl
    character(len=256)  :: section

    section = trim(codename)//'-Numerics'

    !! ------------------------------------------------------
    !! Time scheme ------------------------------------------
    !! ------------------------------------------------------
    obj_time_scheme%warning_message = 'none'
    obj_time_scheme%error_message   = 'none'
    obj_time_scheme%description     = 'none'
    obj_time_scheme%solver_type     = 'RK2'

    call reg%add(trim(section), 'time-scheme', obj_time_scheme%solver_type,            &
                 'RK2',    'Time integration solver', 'euler, RK2, RK3', .true.)

    ! Stability coefficients and related options
    call reg%add(trim(section), 'cfl', obj_time_scheme%cfl,                            &
                 '0.5',    'CFL number', '> 0', .true.)
    call reg%add(trim(section), 'dt-max', obj_time_scheme%dt_max,                      &
                 '1e-4',   'Ceiling on the time step [s], applied after the CFL factor', &
                 '> 0', .false.)
    call reg%add(trim(section), 'tau-factor', obj_time_scheme%tau_factor,              &
                 '1.0',    'Ceiling on the time step as a multiple of the particle relaxation time, '// &
                           'coupled runs only (0 = off)', '>= 0', .false.)
    call reg%add(trim(section), 'cfl-rise-threshold', obj_time_scheme%cfl_rampa_iter,  &
                 '0',      'CFL rise threshold', '>= 0', .false.)

    ! Time-accurate switch
    call reg%add(trim(section), 'time-accurate', obj_time_scheme%time_accurate,        &
                 '.true.', 'Time accurate switch', 'logical', .true.)

    ! Implicit residual smoothing --------------------------
    obj_irs%description     = 'none'
    obj_irs%warning_message = 'none'
    obj_irs%error_message   = 'none'
    call reg%add(trim(section), 'irs', obj_irs%enabled,                                &
                 '.false.', 'Implicit Residual Smoothing', 'logical', .false.)
    call reg%add(trim(section), 'irs-beta', obj_irs%beta,                              &
                 '0.5',     'IRS beta parameter', '>= 0', .false.)

    !! ------------------------------------------------------
    !! Space scheme -----------------------------------------
    !! ------------------------------------------------------
    obj_space_scheme%warning_message      = 'none'
    obj_space_scheme%error_message        = 'none'
    obj_space_scheme%description          = 'none'
    obj_space_scheme%space_reconstruction = 'first-order'
    obj_space_scheme%flux_limiter         = 'none'
    obj_space_scheme%shock_detector       = 'none'

    call reg%add(trim(section), 'space-reconstruction', obj_space_scheme%space_reconstruction, &
                 'first-order', 'Space reconstruction method', 'MUSCL, first-order', .true.)
    call reg%add(trim(section), 'flux-limiter', obj_space_scheme%flux_limiter,          &
                 'none',   'Flux limiter for space reconstruction',                     &
                 'minmod, vanalbada, vanleer, ospre, umist, osher, sweby, mc, koren, '// &
                 'superbee, none', .false.)

    ! Shock detector ---------------------------------------
    call reg%add(trim(section), 'shock-detector', obj_space_scheme%shock_detector,      &
                 'none',   'Shock detector method', 'Jameson, none', .false.)

    ! Riemann solver ---------------------------------------
    call reg%add(trim(section), 'riemann-solver', obj_time_scheme%riemann,              &
                 '',       'Riemann solver (empty: Saurel for MK, Rusanov for IG and AG)', &
                 'Saurel, Rusanov, HLLE', .false.)

    ! Multigrid levels -------------------------------------
    call Register_Multigrid_Levels(nmgl)

  end subroutine Register_Numerics


  subroutine Register_Multigrid_Levels(nmgl)
    use ICE_Parameters_m, only: codename
    use IR_Precision,     only: str
    implicit none
    integer, intent(in) :: nmgl
    character(len=256)  :: section
    integer :: m

    allocate(obj_multigrid%iter_threshold(nmgl))
    obj_multigrid%iter_threshold = 1000000000

    section = trim(codename)//'-Multigrid'

    call reg%add(trim(section), 'levels', obj_multigrid%MGL, &
                 '1', 'Number of grid levels; every block dimension must be '// &
                      'divisible by 2^(levels-1)', '> 0', .false.)
    obj_multigrid%MGL = nmgl

    do m = 1, nmgl
      call reg%add(trim(section), 'level'//trim(str(.true.,m))//'-iter', &
                   obj_multigrid%iter_threshold(m),                      &
                   '1000000000', 'Iterations for multigrid level '//trim(str(.true.,m)), &
                   '> 0', .false.)
    end do

  end subroutine Register_Multigrid_Levels


  subroutine Register_Families(ngroups)
    use ICE_Parameters_m, only: codename
    use IR_Precision,     only: str
    implicit none
    integer, intent(in) :: ngroups
    character(len=256)  :: section
    integer :: p

    allocate(obj_time_scheme%model(ngroups))
    obj_time_scheme%model = ''

    do p = 1, ngroups
      section = trim(codename)//'-Family'//trim(str(.true.,p))
      call reg%add(trim(section), 'closure', obj_time_scheme%model(p), &
                   '', 'Kinetic closure for this family', 'MK, IG, AG', .true.)
    end do

  end subroutine Register_Families

end module ICE_Read_Numerics
