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
    use ICE_Load_Table,     only: Load_Table
    use ICE_IO_Probes,      only: Setup_Probes
    use ICE_Mod_BC_Fluxes
    use ICE_Mod_Phase
    use ICE_Mod_Multigrid,  only: Setup_Multigrid
    use ICE_Mod_MPI,        only: mpi_is_root, partition_blocks
    use ICE_Mod_GhostExchange, only: build_ghost_schedule
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

    ! Load optional temperature-dependent property tables (printed here, after check section)
    call Load_Table()

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
        call read_bck(sim%ODP(1))
        call copyORION(IOgas, sim%OCP)
      else
        call read_bck(sim%ODP(1), sim%OCP)
      end if
    else
      call read_bck(sim%ODP(1))
    end if

    ! Allocate data structure for fine level
    call allocate_data(sim%domain(1), coupled, sim%ODP(1))

    ! Setup metrics for fine level
    call setup_metrics(sim%domain(1), sim%ODP(1))

    ! Print simulation info onto the logfile/shell
    if (mpi_is_root) call print_simulation_info()

    ! Setup boundaries (fine level only)
    call Setup_BC(sim%domain(1))

    ! Distribute the blocks over the MPI ranks (every rank keeps the whole domain but
    ! updates only its own blocks) and build the halo exchange for the ghost fill
    call partition_blocks(sim%domain(1)%nb, &
                          [(product(sim%domain(1)%blk(b)%dim), b = 1, sim%domain(1)%nb)])
    call build_ghost_schedule(sim%domain(1))

    ! Setup gaseous phase (fine level only)
    if (coupled) call setup_gas(sim%domain(1), sim%OCP)

    ! Setup condensed phase (fine level)
    call setup_cond(sim%domain(1), sim%ODP(1))

    ! With multigrid, allocate coarse grid levels and compute their metrics
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


  contains


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
      if (index(obj_space_scheme%space_reconstruction, 'MUSCL') == 0) then
        write(*,'(A)') " - Space   --> I order"
      else if (obj_space_scheme%SD) then
        write(*,'(A)') " - Space   --> MUSCL-SD with "//trim(obj_space_scheme%flux_limiter)//" flux limiter"
      else
        write(*,'(A)') " - Space   --> MUSCL with "//trim(obj_space_scheme%flux_limiter)//" flux limiter"
      end if
      if (trim(obj_time_scheme%solver_type) == '1') then
        write(*,'(A)') " - Time    --> Explicit Euler"
      else
        write(*,'(A)') " - Time    --> Explicit Runge-Kutta "//trim(obj_time_scheme%solver_type)
      end if
      if (coupled) then
        write(*,'(A)') " - Drag    --> "//trim(obj_time_scheme%drag)
        write(*,'(A)') " - Heat    --> "//trim(obj_time_scheme%heat)
      end if

    end subroutine print_simulation_info


  end subroutine ICE_setup


end module ICE_Wrap_Setup
