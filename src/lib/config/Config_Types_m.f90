module ICE_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m

  implicit none
  private

  !! ------------------------------------------------------
  !! Simulation Parameters --------------------------------
  !! ------------------------------------------------------
  type :: simulation_parameters_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    logical   :: newrun          ! Restart flag (true = new run)
    real(R8)  :: res_threshold   ! Min residual to stop execution
    real(R8)  :: time_threshold  ! Max physical time to stop execution
    integer   :: iter_threshold  ! Max iterations to stop execution
    ! Computed at runtime (not registered)
    logical   :: owcoupled = .false.
    logical   :: twcoupled = .false.
    ! Useful variables
    integer   :: iter_general   = 0
    integer   :: iter_from_call = 0
    real(R8)  :: time_from_call = 0._R8
    integer   :: nthreads
    real(R8)  :: cputime(2)
    integer   :: TODO
    logical   :: HYDRA_time_accurate = .false.
    logical   :: HYDRA_postprocess   = .false.
    logical   :: HYDRA_MG            = .false.
    ! Global residual (L2 norm, updated each iteration)
    real(R8)  :: residuotot(5) = 0._R8
  end type simulation_parameters_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Input-Output -----------------------------------------
  !! ------------------------------------------------------
  type :: io_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    integer             :: sol_diter, bck_diter, shell_diter, ini_diter, res_diter
    real(R8)            :: sol_dtime, bck_dtime
    logical             :: sol_overwrite, bck_overwrite
    character(len=llen) :: sol_format        ! e.g. 'tecplot ascii'
    character(len=llen) :: bck_format        ! e.g. 'native binary'
    character(len=clen) :: sol_fmt(2)        ! Parsed: (writer, mode), set by Assign_Setup
    character(len=clen) :: bck_fmt(2)        ! Parsed: (writer, mode), set by Assign_Setup
    character(4)        :: extension
    character(len=hlen) :: gaspath           ! Gas-phase solution path (restart/coupling)
    integer             :: init              ! Initialisation flag
    ! Useful variables
    character(len=hlen) :: Ovarnames
    integer             :: Onvar
    real(R8)            :: IOtime
    ! Residuals file
    integer             :: unitRES
    character(len=llen) :: unitRES_format
  end type io_t

  type :: io_probes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: file
    character(len=hlen) :: varnames
    integer             :: iloc(4)
    real(R8)            :: loc(3)
  end type io_probes_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!! NUMERICAL SCHEME !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Time Scheme ------------------------------------------
  !! ------------------------------------------------------
  type :: time_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS (global)
    real(R8) :: cfl            ! CFL stability parameter
    real(R8) :: dt_max         ! Ceiling on the time step, before the CFL factor
    integer  :: cfl_rampa_iter ! Iteration at which CFL ramp starts
    logical  :: time_accurate  ! Time-accurate integration flag
    character(len=llen) :: solver_type   ! Time integrator: '1'=Euler, '2'=RK2, '3'=RK3
    character(len=llen) :: drag          ! Drag model (global, same for all families)
    integer  :: dragSelect     ! Drag model as the selector Lib_Drag dispatches on
    character(len=llen) :: heat          ! Heat transfer model (global, same for all families)
    ! Per-family (only model type differs across families)
    character(len=llen), allocatable :: model(:)
  end type time_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Implicit Residual Smoothing --------------------------
  !! ------------------------------------------------------
  type :: irs_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8) :: beta    = 0.5_R8   ! Jacobi smoothing coefficient
    ! Useful variables
    logical  :: enabled = .false.
  end type irs_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Space Scheme -----------------------------------------
  !! ------------------------------------------------------
  type :: space_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    character(len=llen) :: space_reconstruction  ! 'MUSCL', 'MUSCL-SD', or empty (1st order)
    character(len=llen) :: flux_limiter          ! 'VANLEER', 'MINMOD', etc.
    ! Useful variables
    logical :: SD = .false.
  end type space_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Multigrid --------------------------------------------
  !! ------------------------------------------------------
  type :: multigrid_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    integer              :: MGL = 1             ! Number of multigrid levels
    integer, allocatable :: iter_threshold(:)   ! Iterations for each level
    ! Useful variables
    integer :: MG_level  = 1      ! Current level being solved
    logical :: change_MG = .false.
  end type multigrid_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! BC ---------------------------------------------------
  !! ------------------------------------------------------
  type :: io_bc_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! Useful variables
    logical, allocatable :: viscous_flag(:,:)
    logical, allocatable :: coupling_flag(:,:)
  end type io_bc_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! PHYSICS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Condensed Phase (ICE-specific) -----------------------
  !! ------------------------------------------------------
  type :: condensed_phase_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8) :: rho_al = 2700._R8    ! Particle density            [kg/m^3]
    real(R8) :: cs_al  = 1598._R8    ! Particle specific heat      [J/(kg K)]
    real(R8) :: lv_al  = 10.8e6_R8   ! Particle latent heat        [J/kg]
    real(R8) :: q_al   = 9.53e6_R8   ! Particle combustion energy  [J/kg]
    real(R8) :: emiss  = 1._R8       ! Surface emissivity          [-]
    ! Table-based properties rho(T) and cs(T) (optional, loaded by Load_Table)
    logical                   :: use_table = .false.
    integer                   :: T_min = 0, T_max = 0
    real(R8), allocatable     :: rho_tab(:)   ! indexed T_min:T_max
    real(R8), allocatable     :: cs_tab(:)    ! indexed T_min:T_max
  end type condensed_phase_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! MODULE-LEVEL INSTANCES !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(simulation_parameters_t),  public :: obj_sim_param
  type(io_t),                      public :: obj_io
  type(io_probes_t), allocatable,  public :: obj_io_probes(:)
  type(io_bc_t),                   public :: obj_io_bc
  type(time_scheme_t),             public :: obj_time_scheme
  type(irs_t),                     public :: obj_irs
  type(space_scheme_t),            public :: obj_space_scheme
  type(multigrid_t),               public :: obj_multigrid
  type(condensed_phase_t),         public :: obj_condensed

end module ICE_Config_Types_m
