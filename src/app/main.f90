program ICE_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use ICE_Advanced_Types_m,  only: ICE_simulation_type
  use ICE_Config_Types_m,    only: obj_sim_param
  use ICE_Procedures_m,      only: ICE_type
  use ICE_Mod_MPI
#ifdef USE_MPI
  use ICE_Mod_GhostExchange, only: cleanup_ghost_schedule
#endif
  implicit none
  type(ICE_type)            :: ICE
  type(ICE_simulation_type) :: simulation

  ! Initialize MPI environment (no-op if USE_MPI is not defined)
  call mpi_init_env()

#if defined (_OPENMP)
  !$omp parallel
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  if (mpi_is_root) then
    write(*,'(A)') ' Parallel execution'
    write(*,'(A)') ' OpenMP:'
    write(*,'(A,I4)') ' -  Number of threads --> ', obj_sim_param%nthreads
  end if
#else
  if (mpi_is_root) write(*,'(A)') ' Serial execution'
  obj_sim_param%nthreads = 1
#endif

#ifdef USE_MPI
  if (mpi_is_root) then
    write(*,'(A)') ' MPI:'
    write(*,'(A,I4)') ' -  Number of ranks   --> ', mpi_size_
  end if
#endif

  ! Solving with ICE. The timers (when [ICE-IO] timers is on) start at the end
  ! of ICE%setup, time every ICE%solve, and report from ICE%postprocess: the
  ! same three calls hydra makes, so a coupled run is timed the same way.
  call ICE%setup(simulation)

  obj_sim_param%TODO = 1
  do while (obj_sim_param%TODO <= 2)
    call ICE%solve(simulation, Dummy_Function)
    if (obj_sim_param%TODO <= 2) call ICE%postprocess(simulation)
  enddo

  call ICE%postprocess(simulation)

  ! Free persistent MPI requests before finalizing
#ifdef USE_MPI
  call cleanup_ghost_schedule()
#endif

  ! Finalize MPI environment
  call mpi_finalize_env()

contains

  subroutine Dummy_Function
    ! Empty subroutine to be passed as an argument to ICE%solve.
    ! It can be used for user-defined operations during the solution process.
    ! It is used in HYDRA coupling procedures.
  end subroutine Dummy_Function

end program ICE_program