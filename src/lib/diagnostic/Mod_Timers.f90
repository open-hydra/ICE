!>@brief Wall-clock and core-cycle instrumentation of the ICE iteration loop.
!>
!> ICE carried no timers at all until the scaling campaign needed them. The
!> design is MOSE's `MOSE_Mod_Timers` — every rank accumulates the time spent
!> in one iteration, the time blocked on the halo exchange and the time blocked
!> on collectives, and the accumulators are reduced over the ranks — with two
!> additions:
!>
!>   * the core-cycle counter from hydra-MF, so that scaling can be reported
!>     with the processor's frequency scaling divided out;
!>   * a split of the compute into `source`, `flux` and `halo`, which is what
!>     says *where* an ICE iteration goes.
!>
!> Below those three phases sit the regions (`timer_region_begin/end`): the
!> steps of `Explicit_Step` one by one, and the pack, wait and unpack of the
!> halo exchange. They are what tells a per-block barrier from a serial copy,
!> or a message from its packing, and they are printed as the `ICE Detail`
!> line beside the phases. The phases keep their old extent, so a report from
!> before the regions reads the same way.
!>
!> Times are reported as the **maximum** over the ranks, which is the critical
!> path and therefore what sets the time to solution; the spread between the
!> slowest and the mean rank is the load imbalance. Reporting the mean would
!> hide exactly the effect a decomposition study is looking for.
!>
!> Everything is off unless `timers = true` is set in `[ICE-IO]`: with the
!> switch off every entry point returns at once and no clock is read.
module ICE_Mod_Timers
  use iso_fortran_env, only: R8 => real64, I4 => int32, I8 => int64
  use iso_c_binding,   only: c_int, c_int64_t
#ifdef USE_MPI
  use mpi
#endif

  implicit none
  private

  !> Core-cycles, so that scaling can be reported with the clock divided out.
  !> A single active core turbos well above the all-core clock, so an efficiency
  !> built from seconds mixes the chip's power management into the solver's
  !> parallel behaviour. Perfect parallelism keeps TOTAL core-cycles per
  !> iteration constant however many cores are used, so cycles are the
  !> clock-free measure. This replaces MOSE's ballast method, which hydra-MF
  !> found produced parallel efficiencies above 100 %: the spinner is more
  !> power-hungry than a memory-bound solver and pushed the 1-core anchor below
  !> the full-node clock instead of pinning it there. Implementation and its
  !> pitfalls: `ice_cycles.c`.
#ifdef ICE_HAVE_CYCLES
  interface
    subroutine ice_cyc_open_thread(slot) bind(C, name='ice_cyc_open_thread')
      import :: c_int
      integer(c_int), intent(in) :: slot
    end subroutine ice_cyc_open_thread
    subroutine ice_cyc_count(n) bind(C, name='ice_cyc_count')
      import :: c_int
      integer(c_int), intent(out) :: n
    end subroutine ice_cyc_count
    subroutine ice_cyc_read_all(total) bind(C, name='ice_cyc_read_all')
      import :: c_int64_t
      integer(c_int64_t), intent(out) :: total
    end subroutine ice_cyc_read_all
  end interface
#endif

  logical     :: on        = .false.  !< master switch, from [ICE-IO] timers
  logical     :: cyc_ok    = .false.  !< counter available on this build/kernel
  integer(I8) :: c_win_beg = 0_I8     !< cycles at the start of the window

  !> The regions of an iteration. `dt` to `diag` are the steps of Explicit_Step
  !> in the order they run (each once per family and stage, `dt`, `copy` and
  !> `diag` once per step); `pack`, `wait`, `unpack` are the parts of the halo
  !> exchange, inside `ghost1`; `sync` is the three collectives, the same time
  !> the `collective wait` figure reports.
  integer, parameter, public :: TR_DT = 1, TR_COPY = 2, TR_ZERO = 3, TR_SOURCE = 4, &
                                TR_GHOST1 = 5, TR_GHOST2 = 6, TR_BOUND = 7, TR_FLUX = 8, &
                                TR_RESID = 9, TR_IRS = 10, TR_UPDATE = 11, TR_DIAG = 12, &
                                TR_PACK = 13, TR_WAIT = 14, TR_UNPACK = 15, TR_SYNC = 16
  integer, parameter, public :: NREG = 16
  character(len=6), parameter :: reg_name(NREG) = [character(len=6) :: &
    'dt', 'copy', 'zero', 'source', 'ghost1', 'ghost2', 'bound', 'flux', &
    'resid', 'irs', 'update', 'diag', 'pack', 'wait', 'unpack', 'sync']

  !> Reduced timings for one window of iterations, per iteration.
  type :: stats_type
    real(R8) :: tmax = 0.0_R8       !< iteration time on the slowest rank
    real(R8) :: tmin = 0.0_R8       !< iteration time on the fastest rank
    real(R8) :: tavg = 0.0_R8       !< iteration time, mean over ranks
    real(R8) :: imbalance = 0.0_R8  !< (tmax-tavg)/tavg, %
    real(R8) :: commfrac = 0.0_R8   !< time blocked on the halo exchange, %
    real(R8) :: syncfrac = 0.0_R8   !< time blocked on collectives, %
    real(R8) :: wmax = 0.0_R8       !< compute time (iteration minus waits), max
    real(R8) :: wmin = 0.0_R8       !< compute time, min
    real(R8) :: wavg = 0.0_R8       !< compute time, mean
    real(R8) :: wspread = 0.0_R8    !< (wmax-wavg)/wavg, %
    integer  :: rmax = 0            !< rank holding wmax
    integer  :: rmin = 0            !< rank holding wmin
    real(R8) :: source = 0.0_R8     !< source terms, max over ranks
    real(R8) :: flux = 0.0_R8       !< transport operator, max over ranks
    real(R8) :: halo = 0.0_R8       !< ghost and boundary fill, max over ranks
    real(R8) :: reg(NREG) = 0.0_R8  !< the regions, max over ranks
  end type stats_type

  real(R8) :: t_iter_beg = 0.0_R8, t_comm_beg = 0.0_R8, t_sync_beg = 0.0_R8
  real(R8) :: t_src_beg  = 0.0_R8, t_flx_beg  = 0.0_R8, t_halo_beg = 0.0_R8
  real(R8) :: t_reg_beg(NREG) = 0.0_R8

  real(R8) :: t_iter_acc = 0.0_R8, t_comm_acc = 0.0_R8, t_sync_acc = 0.0_R8
  real(R8) :: t_src_acc  = 0.0_R8, t_flx_acc  = 0.0_R8, t_halo_acc = 0.0_R8
  real(R8) :: t_reg_acc(NREG) = 0.0_R8
  integer  :: n_iter_acc = 0

  real(R8) :: t_run_beg = 0.0_R8
  real(R8) :: t_run_iter = 0.0_R8, t_run_comm = 0.0_R8, t_run_sync = 0.0_R8
  real(R8) :: t_run_src  = 0.0_R8, t_run_flx  = 0.0_R8, t_run_halo = 0.0_R8
  real(R8) :: t_run_reg(NREG) = 0.0_R8
  integer  :: n_run = 0

  !> Cells in the whole domain, for the throughput figure in the summary.
  integer(I8), public :: n_cells_total = 0_I8

  public :: timers_enabled, timer_wtime
  public :: timer_run_begin, timer_report, timer_summary
  public :: timer_iter_begin, timer_iter_end
  public :: timer_comm_begin, timer_comm_end
  public :: timer_sync_begin, timer_sync_end
  public :: timer_source_begin, timer_source_end
  public :: timer_flux_begin, timer_flux_end
  public :: timer_halo_begin, timer_halo_end
  public :: timer_region_begin, timer_region_end

contains

  pure logical function timers_enabled()
    timers_enabled = on
  end function timers_enabled


  !> Wall-clock reading, in seconds. MPI_WTIME when available (monotonic and
  !> consistent across the ranks of a job), SYSTEM_CLOCK otherwise.
  real(R8) function timer_wtime()
#ifdef USE_MPI
    timer_wtime = MPI_WTIME()
#else
    integer(I8) :: count, rate
    call system_clock(count, rate)
    timer_wtime = real(count, R8) / real(rate, R8)
#endif
  end function timer_wtime


  !> Mark the start of the iteration loop, i.e. the end of set-up. `enable` is
  !> the `[ICE-IO] timers` key; everything below is inert when it is false.
  subroutine timer_run_begin(enable, ncells_total)
    ! The `!$` sentinel, not #ifdef USE_OPENMP: USE_OPENMP is a CMake option and
    ! is never passed to the compiler as a definition, so an #ifdef on it is
    ! always false -- which would give every thread slot 0 and count one.
    !$ use omp_lib, only: omp_get_thread_num
    logical,     intent(in) :: enable
    integer(I8), intent(in) :: ncells_total
    integer(c_int)     :: slot, nc
    integer(c_int64_t) :: c0

    on = enable
    if (.not. on) return

    n_cells_total = ncells_total
    t_run_beg = timer_wtime()

#ifdef ICE_HAVE_CYCLES
    ! Attach a counter to every thread of the pool. Done here, after set-up, so
    ! the pool already has its full complement; the counters then persist for
    ! the run and each window is read as a difference.
    !$omp parallel private(slot)
    slot = 0_c_int
    !$ slot = int(omp_get_thread_num(), c_int)
    call ice_cyc_open_thread(slot)
    !$omp end parallel

    call ice_cyc_count(nc)
    cyc_ok = (nc > 0)
    if (cyc_ok) then
      call ice_cyc_read_all(c0)
      c_win_beg = int(c0, I8)
    end if
#else
    slot = 0_c_int; nc = 0_c_int; c0 = 0_c_int64_t
    cyc_ok = .false.
#endif

    t_run_iter = 0.0_R8; t_run_comm = 0.0_R8; t_run_sync = 0.0_R8
    t_run_src  = 0.0_R8; t_run_flx  = 0.0_R8; t_run_halo = 0.0_R8
    t_run_reg  = 0.0_R8
    n_run = 0
    t_iter_acc = 0.0_R8; t_comm_acc = 0.0_R8; t_sync_acc = 0.0_R8
    t_src_acc  = 0.0_R8; t_flx_acc  = 0.0_R8; t_halo_acc = 0.0_R8
    t_reg_acc  = 0.0_R8
    n_iter_acc = 0
  end subroutine timer_run_begin


  subroutine timer_iter_begin()
    if (on) t_iter_beg = timer_wtime()
  end subroutine timer_iter_begin

  subroutine timer_iter_end()
    real(R8) :: dt
    if (.not. on) return
    dt = timer_wtime() - t_iter_beg
    t_iter_acc = t_iter_acc + dt; t_run_iter = t_run_iter + dt
    n_iter_acc = n_iter_acc + 1;  n_run = n_run + 1
  end subroutine timer_iter_end

  !> Blocked on the neighbour exchange. A subset of the iteration, and of `halo`.
  subroutine timer_comm_begin()
    if (on) t_comm_beg = timer_wtime()
  end subroutine timer_comm_begin

  subroutine timer_comm_end()
    real(R8) :: dt
    if (.not. on) return
    dt = timer_wtime() - t_comm_beg
    t_comm_acc = t_comm_acc + dt; t_run_comm = t_run_comm + dt
  end subroutine timer_comm_end

  !> Blocked on a collective: the dt minimum, the residual sum, the stop flag.
  subroutine timer_sync_begin()
    if (on) t_sync_beg = timer_wtime()
  end subroutine timer_sync_begin

  subroutine timer_sync_end()
    real(R8) :: dt
    if (.not. on) return
    dt = timer_wtime() - t_sync_beg
    t_sync_acc = t_sync_acc + dt; t_run_sync = t_run_sync + dt
    t_reg_acc(TR_SYNC) = t_reg_acc(TR_SYNC) + dt; t_run_reg(TR_SYNC) = t_run_reg(TR_SYNC) + dt
  end subroutine timer_sync_end

  subroutine timer_source_begin()
    if (on) t_src_beg = timer_wtime()
  end subroutine timer_source_begin

  subroutine timer_source_end()
    real(R8) :: dt
    if (.not. on) return
    dt = timer_wtime() - t_src_beg
    t_src_acc = t_src_acc + dt; t_run_src = t_run_src + dt
  end subroutine timer_source_end

  subroutine timer_flux_begin()
    if (on) t_flx_beg = timer_wtime()
  end subroutine timer_flux_begin

  subroutine timer_flux_end()
    real(R8) :: dt
    if (.not. on) return
    dt = timer_wtime() - t_flx_beg
    t_flx_acc = t_flx_acc + dt; t_run_flx = t_run_flx + dt
  end subroutine timer_flux_end

  !> Whole of the ghost and boundary fill, of which `comm` is the MPI part.
  subroutine timer_halo_begin()
    if (on) t_halo_beg = timer_wtime()
  end subroutine timer_halo_begin

  subroutine timer_halo_end()
    real(R8) :: dt
    if (.not. on) return
    dt = timer_wtime() - t_halo_beg
    t_halo_acc = t_halo_acc + dt; t_run_halo = t_run_halo + dt
  end subroutine timer_halo_end

  !> One region of the iteration; regions nest inside the phases and do not
  !> overlap one another, so their sum is the time of the steps they name.
  subroutine timer_region_begin(id)
    integer, intent(in) :: id
    if (on) t_reg_beg(id) = timer_wtime()
  end subroutine timer_region_begin

  subroutine timer_region_end(id)
    integer, intent(in) :: id
    real(R8) :: dt
    if (.not. on) return
    dt = timer_wtime() - t_reg_beg(id)
    t_reg_acc(id) = t_reg_acc(id) + dt; t_run_reg(id) = t_run_reg(id) + dt
  end subroutine timer_region_end


  !> Reduce the accumulators over the ranks. Collective: every rank must call
  !> it on the same iteration.
  function reduce_times(t_iter, t_comm, t_sync, t_src, t_flx, t_hal, t_reg, n) result(s)
    use ICE_Mod_MPI, only: mpi_size_, mpi_allreduce_sum_r8, &
                           mpi_allreduce_min_r8, mpi_allreduce_max_r8, &
                           mpi_allreduce_max_r8_array, mpi_gather_r8
    real(R8), intent(in) :: t_iter, t_comm, t_sync, t_src, t_flx, t_hal
    real(R8), intent(in) :: t_reg(NREG)
    integer,  intent(in) :: n
    type(stats_type)     :: s
    real(R8) :: tsum, csum, ssum, t_work
    real(R8), allocatable :: work(:)
    integer  :: r, ni

    ni = max(n, 1)

    ! Time this rank spent on its own cells: everything not blocked waiting for
    ! a neighbour's halo or for a collective.
    t_work = t_iter - t_comm - t_sync

    call mpi_allreduce_max_r8(t_iter, s%tmax)
    call mpi_allreduce_min_r8(t_iter, s%tmin)
    call mpi_allreduce_sum_r8(t_iter, tsum)
    call mpi_allreduce_sum_r8(t_comm, csum)
    call mpi_allreduce_sum_r8(t_sync, ssum)
    call mpi_allreduce_max_r8(t_src,  s%source)
    call mpi_allreduce_max_r8(t_flx,  s%flux)
    call mpi_allreduce_max_r8(t_hal,  s%halo)
    s%reg = t_reg
    call mpi_allreduce_max_r8_array(s%reg, NREG)

    allocate(work(max(mpi_size_, 1)))
    work = 0.0_R8
    call mpi_gather_r8(t_work, work)

    s%tavg = tsum / real(mpi_size_, R8)
    if (s%tavg > 0.0_R8) s%imbalance = (s%tmax - s%tavg) / s%tavg * 100.0_R8
    if (tsum   > 0.0_R8) s%commfrac  = csum / tsum * 100.0_R8
    if (tsum   > 0.0_R8) s%syncfrac  = ssum / tsum * 100.0_R8

    s%wmin = work(1); s%wmax = work(1)
    do r = 1, mpi_size_
      if (work(r) > s%wmax) then; s%wmax = work(r); s%rmax = r - 1; end if
      if (work(r) < s%wmin) then; s%wmin = work(r); s%rmin = r - 1; end if
      s%wavg = s%wavg + work(r)
    end do
    s%wavg = s%wavg / real(mpi_size_, R8)
    if (s%wavg > 0.0_R8) s%wspread = (s%wmax - s%wavg) / s%wavg * 100.0_R8

    deallocate(work)

    ! Per-iteration figures
    s%tmax = s%tmax / real(ni, R8); s%tmin = s%tmin / real(ni, R8)
    s%tavg = s%tavg / real(ni, R8)
    s%wmax = s%wmax / real(ni, R8); s%wmin = s%wmin / real(ni, R8)
    s%wavg = s%wavg / real(ni, R8)
    s%source = s%source / real(ni, R8)
    s%flux   = s%flux   / real(ni, R8)
    s%halo   = s%halo   / real(ni, R8)
    s%reg    = s%reg    / real(ni, R8)
  end function reduce_times


  !> Sum a per-rank cycle count over the communicator, as a real so that the
  !> 64-bit totals at high core counts stay exact enough to print.
  real(R8) function sum_cycles_over_ranks(c) result(total)
    integer(I8), intent(in) :: c
    real(R8) :: local
#ifdef USE_MPI
    real(R8) :: g
    integer  :: ierr
#endif
    local = real(c, R8)
#ifdef USE_MPI
    call MPI_ALLREDUCE(local, g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    total = g
#else
    total = local
#endif
  end function sum_cycles_over_ranks


  !> Percentage of `whole`, guarded against a zero denominator.
  pure real(R8) function pct(part, whole)
    real(R8), intent(in) :: part, whole
    if (whole > 0.0_R8) then
      pct = part / whole * 100.0_R8
    else
      pct = 0.0_R8
    end if
  end function pct


  !> Report the current window and reset it. Collective: every rank must call.
  subroutine timer_report(iter)
    use ICE_Mod_MPI, only: mpi_is_root, mpi_size_
    integer, intent(in) :: iter
    type(stats_type) :: s
    integer(I8) :: c_now, c_win
    real(R8)    :: cyc_iter, ghz
    integer     :: r

    ! n_iter_acc == 0 means this window has already been reported; say nothing
    ! rather than print a window of zero iterations.
    if (.not. on .or. n_iter_acc == 0) return
    s = reduce_times(t_iter_acc, t_comm_acc, t_sync_acc, &
                     t_src_acc, t_flx_acc, t_halo_acc, t_reg_acc, n_iter_acc)

    ! Total core-cycles spent in this window, summed over every rank and the
    ! threads it owns. Divided by the iterations in the window this is the
    ! clock-free cost of an iteration: compare it across core counts and the
    ! ratio is parallel efficiency with the frequency scaling removed.
    cyc_iter = 0.0_R8
    ghz      = 0.0_R8
    if (cyc_ok) then
#ifdef ICE_HAVE_CYCLES
      block
        integer(c_int64_t) :: craw
        call ice_cyc_read_all(craw)
        c_now = int(craw, I8)
      end block
#else
      c_now = c_win_beg
#endif
      c_win = c_now - c_win_beg
      c_win_beg = c_now
      cyc_iter = sum_cycles_over_ranks(c_win) / real(max(n_iter_acc, 1), R8)
      if (s%tmax > 0.0_R8) ghz = cyc_iter / s%tmax / 1.0e9_R8
    end if

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A,ES11.4,A,ES11.4,A,ES11.4,A,F7.1,A,F7.1,A,F7.1,A)') &
        ' ICE Timing | Iter ', iter, ' | ', max(n_iter_acc, 1), &
        ' iters | wall/iter ', s%tmax, &
        ' s | rank min ', s%tmin, ' avg ', s%tavg, &
        ' | imbalance ', s%imbalance, &
        ' % | exchange wait ', s%commfrac, &
        ' % | collective wait ', s%syncfrac, ' %'

      write(*,'(A,ES11.4,A,F5.1,A,ES11.4,A,F5.1,A,ES11.4,A,F5.1,A)') &
        ' ICE Phases | source ', s%source, ' s (', pct(s%source, s%tmax), &
        ' %) | flux ', s%flux, ' s (', pct(s%flux, s%tmax), &
        ' %) | halo ', s%halo, ' s (', pct(s%halo, s%tmax), ' %)'

      ! Seconds per iteration of every region, max over ranks; the harness
      ! splits this line on blanks into name/value pairs.
      write(*,'(A,16(1X,A,1X,ES10.3))') &
        ' ICE Detail |', (trim(reg_name(r)), s%reg(r), r = 1, NREG)

      if (mpi_size_ > 1) &
        write(*,'(A,ES11.4,A,I0,A,ES11.4,A,I0,A,ES11.4,A,F7.1,A)') &
          ' ICE Ranks  | compute/iter max ', s%wmax, &
          ' s (rank ', s%rmax, ') | min ', s%wmin, &
          ' s (rank ', s%rmin, ') | mean ', s%wavg, &
          ' s | spread ', s%wspread, ' %'

      if (cyc_ok) &
        write(*,'(A,ES11.4,A,F9.2,A)') &
          ' ICE Cycles | core-cycles/iter ', cyc_iter, &
          ' | aggregate ', ghz, ' GHz'
    end if

    t_iter_acc = 0.0_R8; t_comm_acc = 0.0_R8; t_sync_acc = 0.0_R8
    t_src_acc  = 0.0_R8; t_flx_acc  = 0.0_R8; t_halo_acc = 0.0_R8
    t_reg_acc  = 0.0_R8
    n_iter_acc = 0
  end subroutine timer_report


  !> End-of-run summary. Collective, so every rank must call it.
  subroutine timer_summary()
    use ICE_Mod_MPI, only: mpi_is_root, mpi_size_
    type(stats_type) :: s
    real(R8) :: elapsed, cells_per_cs
    integer  :: r

    if (.not. on) return

    s = reduce_times(t_run_iter, t_run_comm, t_run_sync, &
                     t_run_src, t_run_flx, t_run_halo, t_run_reg, n_run)
    elapsed = timer_wtime() - t_run_beg

    if (mpi_is_root) then
      cells_per_cs = 0.0_R8
      if (s%tmax > 0.0_R8) &
        cells_per_cs = real(n_cells_total, R8) / s%tmax / real(max(mpi_size_, 1), R8)
      write(*,'(A)') ''
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' ICE Timing'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A,T35,I0)')       '   Iterations', n_run
      write(*,'(A,T35,ES12.5,A)') '   Solver', s%tmax * real(max(n_run,1), R8), ' s'
      write(*,'(A,T35,ES12.5,A)') '   Solver per iteration', s%tmax, ' s'
      write(*,'(A,T35,ES12.5,A)') '   Loop elapsed, with I/O', elapsed, ' s'
      write(*,'(A,T35,ES12.5,A)') '   Source per iteration', s%source, ' s'
      write(*,'(A,T35,ES12.5,A)') '   Flux per iteration', s%flux, ' s'
      write(*,'(A,T35,ES12.5,A)') '   Halo per iteration', s%halo, ' s'
      if (mpi_size_ > 1) then
        write(*,'(A,T35,F7.1,A)') '   Load imbalance', s%imbalance, ' %'
        write(*,'(A,T35,F7.1,A)') '   Exchange wait', s%commfrac, ' %'
        write(*,'(A,T35,F7.1,A)') '   Collective wait', s%syncfrac, ' %'
        write(*,'(A,T35,F7.1,A)') '   Compute spread over ranks', s%wspread, ' %'
      end if
      write(*,'(A,T35,ES12.5)')   '   Cells per rank-second', cells_per_cs
      write(*,'(A)') '   Regions per iteration, max over ranks'
      do r = 1, NREG
        write(*,'(A,A,T35,ES12.5,A)') '     ', reg_name(r), s%reg(r), ' s'
      end do
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ''
    end if
  end subroutine timer_summary

end module ICE_Mod_Timers
