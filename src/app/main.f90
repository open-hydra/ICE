program ICE_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use ICE_Advanced_Types_m,  only: ICE_simulation_type
  use ICE_Config_Types_m,    only: obj_sim_param, obj_io
  use iso_fortran_env,       only: int64
  use ICE_Procedures_m,      only: ICE_type
  use ICE_Mod_MPI
  use ICE_Mod_Timers,        only: timer_run_begin, timer_iter_begin, &
                                   timer_iter_end, timer_report, timer_summary
#ifdef USE_MPI
  use ICE_Mod_GhostExchange, only: cleanup_ghost_schedule
#endif
  implicit none
  type(ICE_type)            :: ICE
  type(ICE_simulation_type) :: simulation
  integer(int64)            :: ncells_total

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

  ! Solving with ICE
  call ICE%setup(simulation)

  ! Timers start after set-up, so the OpenMP pool is at full size when the
  ! per-thread cycle counters are attached and the loop cost excludes the I/O
  ! and the partitioning.
  call count_cells(simulation, ncells_total)
  call timer_run_begin(obj_io%timers, ncells_total)

  obj_sim_param%TODO = 1
  do while (obj_sim_param%TODO <= 2)
    call timer_iter_begin()
    call ICE%solve(simulation, Dummy_Function)
    call timer_iter_end()
    if (obj_sim_param%TODO <= 2) call ICE%postprocess(simulation)
  enddo

  call timer_report(simulation%domain(1)%iter)

  call ICE%postprocess(simulation)
  call timer_summary()

  ! Free persistent MPI requests before finalizing
#ifdef USE_MPI
  call cleanup_ghost_schedule()
#endif

  ! Finalize MPI environment
  call mpi_finalize_env()

contains

  !> Cells this rank owns and cells in the whole domain, for the throughput
  !> figure in the timer summary. Ghost cells are excluded: they are work, but
  !> they are not part of the problem being solved.
  subroutine count_cells(sim, ntot)
    use ICE_Advanced_Types_m, only: ICE_simulation_type
    type(ICE_simulation_type), intent(in)  :: sim
    integer(int64),            intent(out) :: ntot
    integer :: b
    ntot = 0_int64
    do b = 1, sim%domain(1)%nb
      associate (d => sim%domain(1)%blk(b)%dim)
        ntot = ntot + int(d(1), int64) * int(d(2), int64) * int(d(3), int64)
      end associate
    end do
  end subroutine count_cells


  subroutine Dummy_Function
    ! Empty subroutine to be passed as an argument to ICE%solve.
    ! It can be used for user-defined operations during the solution process.
    ! It is used in HYDRA coupling procedures.
  end subroutine Dummy_Function

end program ICE_program