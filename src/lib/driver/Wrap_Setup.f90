module ICE_Wrap_Setup

  implicit none
  private
  public :: ICE_setup

contains

  subroutine ICE_setup(sim, IOgas)
    use Lib_ORION_data
    use ICE_Advanced_Types_m,  only: ICE_simulation_type
    use ICE_Config_Types_m
    use ICE_Read_Ini,          only: Read_Inifile
    use ICE_Assign_Setup,      only: Assign_Setup
    use ICE_IO_Solution
    use ICE_Global_m
    use ICE_Mod_Allocate_Data
    use ICE_Mod_Metrics
    use ICE_IO_BC,          only: Setup_BC
    use ICE_Setup_Materials, only: Setup_Materials
    use ICE_Lib_Heat,       only: heat_formula
    use ICE_IO_Probes,      only: Setup_Probes
    use ICE_Mod_BC_Fluxes
    use ICE_Mod_Phase
    use ICE_Mod_Multigrid,  only: Setup_Multigrid
    use ICE_Mod_MPI,        only: mpi_is_root, partition_blocks
    use ICE_Mod_GhostExchange, only: build_ghost_schedule, build_local_bc_index, report_bc_mix
    use ICE_Mod_Timers,     only: timer_run_begin
    use iso_fortran_env,    only: int64
    !$ use omp_lib,         only: omp_set_schedule, omp_sched_static, omp_sched_dynamic, omp_sched_guided
    implicit none
    type(ICE_simulation_type), intent(inout)  :: sim
    type(orion_data), intent(inout), optional :: IOgas
    logical :: coupled
    integer :: ios, b

    ! Print header
    if (mpi_is_root) call Print_Header()

    ! Read and validate input.ini
    call Read_Inifile()

    ! Post-read: assign solvers, compute ncond, detect coupling
    call Assign_Setup()

    ! The materials of the phase file, their models, and the property table
    call Setup_Materials()

    obj_sim_param%HYDRA_MG = (obj_multigrid%MGL > 1)

    coupled = obj_sim_param%owcoupled .or. obj_sim_param%twcoupled

    ! Wire IO procedure pointers and resolve IC file paths
    call input_solution_setup()

    ! Allocate domain and ODP arrays (one per multigrid level)
    allocate(sim%domain(obj_multigrid%MGL))
    allocate(sim%ODP(obj_multigrid%MGL))

    ! Read mesh and initial conditions into fine level
    if (coupled) then
      if (present(IOgas)) then
        call read_ic(sim%ODP(1))
        call copyORION(IOgas, sim%OCP)
      else
        call read_ic(sim%ODP(1), sim%OCP)
      end if
    else
      call read_ic(sim%ODP(1))
    end if

    ! Allocate data structure for fine level
    call allocate_data(sim%domain(1), coupled, sim%ODP(1))

    ! Setup metrics for fine level
    call setup_metrics(sim%domain(1), sim%ODP(1))

    ! Print simulation info onto the logfile/shell
    if (mpi_is_root) call print_simulation_info()

    ! Setup boundaries. Each grid level reads its own file: the fine level
    ! <prefix>bc.txt, a coarse level <prefix>bc<level>.txt (see Setup_Multigrid).
    call Setup_BC(sim%domain(1), 1)

    ! Distribute the blocks over the MPI ranks (every rank keeps the whole domain but
    ! updates only its own blocks) and build the halo exchange for the ghost fill
    call partition_blocks(sim%domain(1)%nb, &
                          [(product(sim%domain(1)%blk(b)%dim), b = 1, sim%domain(1)%nb)])
    call build_ghost_schedule(sim%domain(1), 1)
    call build_local_bc_index(sim%domain(1))

    ! Setup gaseous phase (fine level only)
    if (coupled) call setup_gas(sim%domain(1), sim%OCP)

    ! Setup condensed phase (fine level)
    call setup_cond(sim%domain(1), sim%ODP(1))

    ! With multigrid, build the coarse levels: grid, metrics, boundaries, halo
    ! schedule and the initial solution restricted onto each of them.
    if (obj_multigrid%MGL > 1) then
      call Setup_Multigrid(sim)
      ! Start solver on coarsest level; level cycling is handled in Wrap_Solve
      obj_multigrid%MG_level = obj_multigrid%MGL
    end if

    ! Fine-level itermax: iter-threshold governs for MGL=1; min with level1-iter when using MG
    if (obj_multigrid%MGL > 1) then
      sim%domain(1)%itermax = min(obj_sim_param%iter_threshold, obj_multigrid%iter_threshold(1))
    else
      sim%domain(1)%itermax = obj_sim_param%iter_threshold
    end if

    call output_solution_setup(sim%ODP(1))

    ! Setup probes (no-op if nprobes == 0)
    call Setup_Probes(sim%domain(1))

    ! Open residuals file (written by root only)
    if (mpi_is_root) &
      open(newunit=obj_io%unitRES, &
           file='OUTPUT/'//trim(ICE_phase_prefix)//'residual-history.dat', &
           status='unknown', form='formatted')

    if (.not. obj_sim_param%newrun) then
      if (mpi_is_root) then
        ios = 0
        do while (ios == 0)
          read(obj_io%unitRES, *, iostat=ios)
        end do
        backspace(obj_io%unitRES)
      end if

      if (obj_time_scheme%time_accurate) then
        sim%domain(1)%time           = sim%ODP(1)%solutiontime
        obj_sim_param%iter_general   = 0
      else
        obj_sim_param%iter_general   = int(sim%ODP(1)%solutiontime)
        sim%domain(1)%time           = -1._8
      end if
    else
      if (obj_time_scheme%time_accurate) then
        sim%domain(1)%time = 0._8
      else
        sim%domain(1)%time = -1.12358132135_8
      end if
      obj_sim_param%iter_general = 0
    end if

    obj_sim_param%time_from_call = sim%domain(1)%time
    sim%domain(1)%iter           = 0
    obj_sim_param%iter_from_call = 0

    write(obj_io%unitRES_format, '(A11,I0,A7)') '(I8,E20.10,', 5, 'E20.10)'

    call Cpu_Time(obj_sim_param%cputime(1))

    ! The timers start here, at the end of set-up, so that the loop cost
    ! excludes the I/O and the partitioning, the cycle counters attach to the
    ! OpenMP pool at its full size, and a coupled run -- hydra calls ICE%setup,
    ! never ICE's own main -- gets them too.
    !> The step's per-cell loops are SCHEDULE(RUNTIME) and ICE sets that schedule
    !> itself: static. OMP_SCHEDULE therefore has no effect on ICE (a dynamic
    !> schedule on these loops dispatches per cell and cost up to 25x at 80
    !> threads); ICE_OMP_SCHEDULE=<kind>[,<chunk>] is the measurement override.
    !> ICE_TILES_PER_THREAD sets the flux kernel's tiles per thread (default 4:
    !> its tile loops are scheduled dynamically). Neither changes a result: which
    !> thread sweeps a cell is not arithmetic.
    block
      character(len=32) :: env
      integer :: status, ios, chunk
      !$ call omp_set_schedule(omp_sched_static, 0)
      call get_environment_variable('ICE_OMP_SCHEDULE', env, status=status)
      if (status == 0) then
        chunk = 0
        if (index(env, ',') > 0) then
          read(env(index(env, ',')+1:), *, iostat=ios) chunk
          if (ios /= 0) chunk = 0
          env = env(1:index(env, ',')-1)
        end if
        select case (trim(adjustl(env)))
        !$ case ('dynamic'); call omp_set_schedule(omp_sched_dynamic, chunk)
        !$ case ('guided');  call omp_set_schedule(omp_sched_guided,  chunk)
        !$ case ('static');  call omp_set_schedule(omp_sched_static,  chunk)
        case default
          if (mpi_is_root) write(*,'(A)') '  OpenMP: ICE_OMP_SCHEDULE not understood, static kept'
        end select
        if (mpi_is_root) write(*,'(A)') '  OpenMP: cell loops scheduled '//trim(env)
      end if
      call get_environment_variable('ICE_TILES_PER_THREAD', env, status=status)
      if (status == 0) then
        read(env, *, iostat=ios) tiles_per_thread
        if (ios /= 0 .or. tiles_per_thread < 1) tiles_per_thread = 1
      end if
      if (status == 0 .and. mpi_is_root) &
        write(*,'(A,I0)') '  OpenMP: flux tiles per thread ', tiles_per_thread
    end block
    if (obj_io%timers) call report_bc_mix(sim%domain(1))
    call timer_run_begin(obj_io%timers, count_cells(sim))


  contains


    !> Cells in the whole domain, for the throughput figure in the timer
    !> summary. Ghost cells are excluded: they are work, but they are not part
    !> of the problem being solved.
    integer(int64) function count_cells(s) result(ntot)
      type(ICE_simulation_type), intent(in) :: s
      integer :: b
      ntot = 0_int64
      do b = 1, s%domain(1)%nb
        associate (d => s%domain(1)%blk(b)%dim)
          ntot = ntot + int(d(1), int64) * int(d(2), int64) * int(d(3), int64)
        end associate
      end do
    end function count_cells


    subroutine Print_Header()
      write(*,*)
      write(*,'(A)') ' ██  ███████  ███████'
      write(*,'(A)') ' ██  ██       ██     '
      write(*,'(A)') ' ██  ██       █████  '
      write(*,'(A)') ' ██  ██       ██     '
      write(*,'(A)') ' ██  ███████  ███████'
      write(*,*)
      write(*,'(A)') ' Integration of a Condensed phase via Eulerian methods'
      write(*,*)
    end subroutine Print_Header


    subroutine print_simulation_info()
      integer :: m

      write(*,*)
      write(*,'(A)') " Checking input file..."
      write(*,'(A)') " - Mesh                   ---> OK"
      write(*,'(A)') " - Initial conditions     ---> OK"
      write(*,'(A)') " - Boundary conditions    ---> OK"
      write(*,'(A)') " - Particle properties    ---> OK"

      write(*,*)
      if (obj_time_scheme%time_accurate) then
        write(*,'(A)', advance='no') " Global time step"
      else
        write(*,'(A)', advance='no') " Local time step"
      end if
      if (obj_sim_param%twcoupled) then
        write(*,'(A)') " 2-way coupled simulation"
      elseif (obj_sim_param%owcoupled) then
        write(*,'(A)') " 1-way coupled simulation"
      else
        write(*,'(A)') " 0-way coupled simulation"
      end if

      write(*,*)
      write(*,'(A)') " ICE phase model:"
      if (npop(1)>0) write(*,'(A,I4)') " - MK particle families  --> ", npop(1)
      if (npop(2)>0) write(*,'(A,I4)') " - IG particle families  --> ", npop(2)
      if (npop(3)>0) write(*,'(A,I4)') " - AG particle families  --> ", npop(3)

      write(*,*)
      write(*,'(A)') " ICE numerical scheme:"
      if (trim(obj_space_scheme%space_reconstruction) /= 'MUSCL') then
        write(*,'(A)') " - Space   --> first order"
      else
        write(*,'(A)') " - Space   --> MUSCL with "//trim(obj_space_scheme%flux_limiter)//" flux limiter"
      end if
      if (obj_space_scheme%SD) &
        write(*,'(A)') " - Shock   --> "//trim(obj_space_scheme%shock_detector)//" detector"
      if (trim(obj_time_scheme%solver_type) == 'euler') then
        write(*,'(A)') " - Time    --> Explicit Euler"
      else
        write(*,'(A)') " - Time    --> Explicit "//trim(obj_time_scheme%solver_type)
      end if
      if (coupled) then
        write(*,'(A)') " - Drag    --> "//trim(obj_time_scheme%drag)
        write(*,'(A)') " - Heat    --> "//trim(obj_time_scheme%heat)//" ("//heat_formula(obj_time_scheme%heatSelect)//")"
        do m = 1, nmat
          write(*,'(A,I0,A)') " - Evap    --> material ", m, " ("//trim(obj_condensed(m)%name)//"): "// &
            trim(obj_condensed(m)%evapWord)
          if (obj_condensed(m)%evapSelect /= 0) &
            write(*,'(A)') " - Interf  --> "//trim(obj_condensed(m)%intfWord)
        enddo
        if (any(obj_condensed(:)%evapSelect /= 0)) &
          write(*,'(A)') " - Blowing --> "//trim(obj_time_scheme%blowing)
      end if
      do m = 1, nmat
        if (obj_condensed(m)%solid) &
          write(*,'(A,I0,A,F0.2,A,F0.2,A,ES12.5,A,F0.2,A)') " - Solid   --> material ", m, " ("// &
            trim(obj_condensed(m)%name)//"): T-melt ", obj_condensed(m)%Tmelt, " K, T-nuc ", obj_condensed(m)%Tnuc, &
            " K, h-fus ", obj_condensed(m)%hFus, " J/kg, cp-solid ", obj_condensed(m)%cpSol, " J/(kg K)"
      enddo

    end subroutine print_simulation_info


  end subroutine ICE_setup


end module ICE_Wrap_Setup
