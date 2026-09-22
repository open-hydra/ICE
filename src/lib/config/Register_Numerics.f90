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

    obj_time_scheme%warning_message = 'none'
    obj_time_scheme%error_message   = 'none'
    obj_time_scheme%description     = 'none'
    obj_time_scheme%drag            = 'None'
    obj_time_scheme%heat            = 'None'
    obj_time_scheme%solver_type     = '2'

    obj_space_scheme%warning_message     = 'none'
    obj_space_scheme%error_message       = 'none'
    obj_space_scheme%description         = 'none'
    obj_space_scheme%space_reconstruction = ''
    obj_space_scheme%flux_limiter        = 'none'

    section = trim(codename)//'-Parameters'

    call reg%add(trim(section), 'cfl',                obj_time_scheme%cfl,            &
                 '0.5',    'CFL stability parameter',            '> 0', .false.)
    call reg%add(trim(section), 'dt-max',             obj_time_scheme%dt_max,         &
                 '1e-4',   'Ceiling on the local time step [s], applied before the CFL '// &
                           'factor: the step never exceeds cfl * dt-max', '> 0', .false.)
    call reg%add(trim(section), 'cfl-rise-threshold', obj_time_scheme%cfl_rampa_iter, &
                 '0',      'CFL ramp start iteration',           '>= 0',.false.)
    call reg%add(trim(section), 'time-accurate',      obj_time_scheme%time_accurate,  &
                 '.true.', 'Time-accurate integration flag',     '',    .false.)
    call reg%add(trim(section), 'irs',      obj_irs%enabled,                         &
                 '.false.', 'Enable implicit residual smoothing', '',   .false.)
    call reg%add(trim(section), 'irs-beta', obj_irs%beta,                             &
                 '0.5',     'IRS Jacobi smoothing coefficient',   '',   .false.)

    section = trim(codename)//'-Scheme'

    call reg%add(trim(section), 'space-reconstruction', obj_space_scheme%space_reconstruction, &
                 '',      'Space reconstruction (MUSCL, MUSCL-SD, or empty)', '', .false.)
    call reg%add(trim(section), 'flux-limiter',          obj_space_scheme%flux_limiter,        &
                 'none',  'Flux limiter (VANLEER, MINMOD, MC, SUPERBEE)',      '', .false.)
    call reg%add(trim(section), 'time',                  obj_time_scheme%solver_type,          &
                 '2',     'Time integrator (1=Euler, 2=RK2, 3=RK3)',          '', .false.)
    call reg%add(trim(section), 'drag',                  obj_time_scheme%drag,                 &
                 'None',  'Drag model (global for all families)',               '', .false.)
    call reg%add(trim(section), 'heat',                  obj_time_scheme%heat,                 &
                 'None',  'Heat transfer model (global for all families)',      '', .false.)

    call Register_Multigrid_Levels(nmgl)

  end subroutine Register_Numerics


  subroutine Register_Multigrid_Levels(nmgl)
    use ICE_Parameters_m, only: codename
    use IR_Precision,     only: str
    implicit none
    integer, intent(in) :: nmgl
    character(len=256)  :: section
    integer :: m

    obj_multigrid%MGL = nmgl
    allocate(obj_multigrid%iter_threshold(nmgl))
    obj_multigrid%iter_threshold = 1000000000

    section = trim(codename)//'-Multigrid'

    do m = 1, nmgl
      call reg%add(trim(section), 'level'//trim(str(.true.,m))//'-iter', &
                   obj_multigrid%iter_threshold(m),                      &
                   '1000000000', 'Max iterations at multigrid level '//trim(str(.true.,m)), &
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
      call reg%add(trim(section), 'model', obj_time_scheme%model(p), &
                   '', 'Particle model (MK/IG/AG)', '', .true.)
    end do

  end subroutine Register_Families

end module ICE_Read_Numerics
