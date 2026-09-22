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
                 '0',      'Ramp the CFL number linearly over this many iterations '// &
                           '(0 = no ramp)',                     '>= 0',.false.)
    call reg%add(trim(section), 'time-accurate',      obj_time_scheme%time_accurate,  &
                 '.true.', 'Advance every cell with the global minimum step (true) or '// &
                           'with its own local step, for steady state (false)', '', .false.)
    call reg%add(trim(section), 'irs',      obj_irs%enabled,                         &
                 '.false.', 'Enable implicit residual smoothing', '',   .false.)
    call reg%add(trim(section), 'irs-beta', obj_irs%beta,                             &
                 '0.5',     'IRS Jacobi smoothing coefficient',   '',   .false.)

    section = trim(codename)//'-Scheme'

    call reg%add(trim(section), 'space-reconstruction', obj_space_scheme%space_reconstruction, &
                 '',      'Space reconstruction: MUSCL, MUSCL-SD (MUSCL with the '// &
                          'density shock detector), or empty for first order', '', .false.)
    call reg%add(trim(section), 'flux-limiter',          obj_space_scheme%flux_limiter,        &
                 'none',  'Flux limiter, used only with MUSCL: IORD, MINMOD, VANALBADA, '// &
                          'VANLEER, OSPRE, UMIST, OSHER, SWEBY, MC, KOREN, SUPERBEE',  '', .false.)
    call reg%add(trim(section), 'time',                  obj_time_scheme%solver_type,          &
                 '2',     'Time integrator: 1 = forward Euler, 2 = SSP-RK2, 3 = SSP-RK3', &
                 '', .false.)
    call reg%add(trim(section), 'drag',                  obj_time_scheme%drag,                 &
                 'None',  'Drag model, global for all families: Newton, Stokes, '// &
                          'Schlichting, Schiller-Naumann, Wen-Yu, Putnam, Clift-Gauvin, '// &
                          'Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, Hermsen', '', .false.)
    call reg%add(trim(section), 'heat',                  obj_time_scheme%heat,                 &
                 'None',  'Heat transfer model, global for all families: Stokes, JAXA1, '// &
                          'JAXA2, JAXA3, Chang, Ranz-Marshall, Kavanau-Drake',          '', .false.)

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
                 '1', 'Number of grid levels. Each coarse level halves every block '// &
                      'dimension, so every block must be divisible by 2^(levels-1)', &
                 '> 0', .false.)
    obj_multigrid%MGL = nmgl

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
                   '', 'Closure for this family: MK (monokinetic, 6 variables), '// &
                   'IG (isotropic Gaussian, 7) or AG (anisotropic Gaussian, 12)', '', .true.)
    end do

  end subroutine Register_Families

end module ICE_Read_Numerics
