!>@brief MPI infrastructure module for ICE.
!> Provides MPI environment management, block partitioning, and
!> global reduction wrappers. When compiled without USE_MPI, all
!> routines are no-ops or identity operations (serial fallback).
module ICE_Mod_MPI
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I4 => int32

  implicit none
  private

  ! --- Public state ---
  integer, public :: mpi_rank_ = 0        !< This process rank (0 in serial)
  integer, public :: mpi_size_ = 1        !< Total number of MPI processes
  logical, public :: mpi_is_root = .true.  !< True on rank 0
#ifdef USE_MPI
  integer, public :: ice_comm = MPI_COMM_WORLD  !< The communicator ICE runs on: the world, or the one a host hands to mpi_init_env
#else
  integer, public :: ice_comm = 0
#endif
  logical, private :: mpi_init_here = .false.  !< ICE called MPI_INIT itself, so it is ICE that finalizes

  ! --- Block-to-rank mapping ---
  integer, allocatable, public :: block_owner(:)     !< block_owner(b) = rank owning block b
  integer, allocatable, public :: local_block_ids(:) !< Global block IDs owned by this rank
  integer, public :: n_local_blocks = 0

  ! --- Public procedures ---
  public :: mpi_init_env, mpi_finalize_env
  public :: is_local_block
  public :: mpi_gather_r8
  public :: partition_blocks, lpt_assign
  public :: mpi_allreduce_sum_r8, mpi_allreduce_min_r8, mpi_allreduce_max_r8
  public :: mpi_allreduce_sum_r8_array, mpi_allreduce_max_r8_array
  public :: mpi_reduce_sum_r8, mpi_reduce_sum_r8_array
  public :: mpi_allreduce_norm2
  public :: mpi_bcast_logical, mpi_bcast_integer
  public :: mpi_abort_all
  public :: check_mpi_error

contains


  !> Initialize MPI environment with thread support. Call early in main program.
  !> Uses MPI_THREAD_FUNNELED: only master thread makes MPI calls. A host program
  !> that has initialized MPI itself may pass the communicator ICE is to run on
  !> (every collective and message of ICE then stays inside it); the default is
  !> MPI_COMM_WORLD. Errors abort the world either way.
  subroutine mpi_init_env(comm)
    integer, intent(in), optional :: comm
#ifdef USE_MPI
    integer :: ierr, provided
    logical :: already

    call MPI_Initialized(already, ierr)
    if (.not. already) then
      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
      call check_mpi_error(ierr)
      if (provided < MPI_THREAD_FUNNELED) then
        write(*,'(A)') ' ERROR: MPI does not support the required threading level (MPI_THREAD_FUNNELED)'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if
      mpi_init_here = .true.
    end if
    if (present(comm)) ice_comm = comm
    call MPI_COMM_RANK(ice_comm, mpi_rank_, ierr)
    call check_mpi_error(ierr)
    call MPI_COMM_SIZE(ice_comm, mpi_size_, ierr)
    call check_mpi_error(ierr)
    mpi_is_root = (mpi_rank_ == 0)
#else
    mpi_rank_ = 0
    mpi_size_ = 1
    mpi_is_root = .true.
#endif
  end subroutine mpi_init_env


  !> Finalize MPI environment. Call at end of main program. MPI is finalized only
  !> when mpi_init_env initialized it (a host program that did finalizes itself).
  subroutine mpi_finalize_env()
#ifdef USE_MPI
    integer :: ierr
    call MPI_BARRIER(ice_comm, ierr)
    if (mpi_init_here) call MPI_FINALIZE(ierr)
#endif
  end subroutine mpi_finalize_env


  !> Check if block b is owned by this rank.
  logical function is_local_block(b)
    integer, intent(in) :: b

    if (.not. allocated(block_owner)) then
      is_local_block = .true.  ! Serial mode or not yet partitioned
      return
    end if
    is_local_block = (block_owner(b) == mpi_rank_)
  end function is_local_block


  !> Partition nb blocks across MPI ranks using greedy load balancing.
  !> blk_ncells(b) = total number of cells in block b.
  !> Blocks are walked largest-first (LPT), so the balance depends only on the cell
  !> counts and not on the order the blocks appear in the mesh file.
  subroutine partition_blocks(nb, blk_ncells, quiet)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    logical, intent(in), optional :: quiet   !< no balance report (an early call; the report comes with the later one)
    ! Local
    integer :: b, nloc
    integer, allocatable :: rank_load(:)

    if (allocated(block_owner)) deallocate(block_owner)
    allocate(block_owner(nb))
    allocate(rank_load(0:mpi_size_-1))

    call lpt_assign(nb, blk_ncells, mpi_size_, block_owner, rank_load)

    if (.not. present(quiet)) then
      call report_partition_balance(nb, blk_ncells, rank_load)
    else if (.not. quiet) then
      call report_partition_balance(nb, blk_ncells, rank_load)
    end if

    ! Build local block list
    nloc = count(block_owner == mpi_rank_)
    n_local_blocks = nloc
    if (allocated(local_block_ids)) deallocate(local_block_ids)
    allocate(local_block_ids(nloc))
    nloc = 0
    do b = 1, nb
      if (block_owner(b) == mpi_rank_) then
        nloc = nloc + 1
        local_block_ids(nloc) = b
      end if
    end do

    deallocate(rank_load)
  end subroutine partition_blocks


  !> Longest-processing-time assignment of n items of weight w to nbin bins:
  !> items walked by descending weight (insertion sort, ties keep their order, so
  !> the result is reproducible), each given to the least-loaded bin (the lowest
  !> index among equals). owner(i) in 0..nbin-1; load(0:nbin-1) the bins' sums.
  !> Used for blocks over ranks and, inside a rank, for blocks over thread groups.
  subroutine lpt_assign(n, w, nbin, owner, load)
    integer, intent(in)  :: n, nbin
    integer, intent(in)  :: w(n)
    integer, intent(out) :: owner(n)
    integer, intent(out) :: load(0:nbin-1)
    integer :: b, i, j, r, tmp
    integer, allocatable :: order(:)

    allocate(order(n))
    load = 0
    do b = 1, n
      order(b) = b
    end do
    do i = 2, n
      tmp = order(i)
      j = i - 1
      do while (j >= 1)
        if (w(order(j)) >= w(tmp)) exit
        order(j+1) = order(j)
        j = j - 1
      end do
      order(j+1) = tmp
    end do

    do i = 1, n
      b = order(i)
      r = minloc(load, dim=1) - 1   ! bin with minimum load (0-indexed)
      owner(b) = r
      load(r) = load(r) + w(b)
    end do
    deallocate(order)
  end subroutine lpt_assign


  subroutine report_partition_balance(nb, blk_ncells, rank_load)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    integer, intent(in) :: rank_load(0:)
    ! Local
    real(R8) :: ideal, eff

    if (mpi_rank_ /= 0) return
    if (mpi_size_ <= 1) return

    ideal = real(sum(blk_ncells), R8) / real(mpi_size_, R8)
    eff   = ideal / real(maxval(rank_load), R8) * 100.0_R8

    write(*,'(A,I0,A,I0,A,F5.1,A)') '  MPI partition: ', nb, ' blocks over ', &
      mpi_size_, ' ranks, balance ', eff, '% of ideal'

    if (mpi_size_ > nb) then
      write(*,'(A,I0,A,I0,A)') '  WARNING: blocks are indivisible - only ', nb, &
        ' of ', mpi_size_, ' ranks have work, the rest idle. Reduce ranks or split the mesh.'
    else if (real(maxval(blk_ncells), R8) > ideal) then
      write(*,'(A,I0,A)') '  NOTE: balance is capped by the largest block (', &
        maxval(blk_ncells), ' cells); only a finer mesh split can improve it.'
    end if
  end subroutine report_partition_balance


  !> MPI_ALLREDUCE with MPI_SUM for a scalar real(R8).
  subroutine mpi_allreduce_sum_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ice_comm, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_sum_r8


  !> MPI_ALLREDUCE with MPI_MIN for a scalar real(R8).
  subroutine mpi_allreduce_min_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ice_comm, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_min_r8


  !> MPI_ALLREDUCE with MPI_MAX for a scalar real(R8).
  !> MPI_GATHER of one real(R8) per rank onto root. Used by the timers to show
  !> the spread of compute time across the ranks, which is the load imbalance a
  !> max-over-ranks figure hides.
  subroutine mpi_gather_r8(local_val, arr)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: arr(:)
#ifdef USE_MPI
    integer :: ierr
    call MPI_GATHER(local_val, 1, MPI_DOUBLE_PRECISION, &
                    arr, 1, MPI_DOUBLE_PRECISION, 0, ice_comm, ierr)
    call check_mpi_error(ierr)
#else
    arr(1) = local_val
#endif
  end subroutine mpi_gather_r8


  subroutine mpi_allreduce_max_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ice_comm, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_max_r8


  !> MPI_ALLREDUCE with MPI_SUM for an array of real(R8).
  !> Works in-place: local_arr is overwritten with the global result.
  subroutine mpi_allreduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(MPI_IN_PLACE, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, ice_comm, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_allreduce_sum_r8_array


  !> MPI_ALLREDUCE with MPI_MAX for an array of real(R8), in place.
  subroutine mpi_allreduce_max_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(MPI_IN_PLACE, arr, n, MPI_DOUBLE_PRECISION, MPI_MAX, ice_comm, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_allreduce_max_r8_array


  !> Distributed L2 norm: local sum of squares + MPI_ALLREDUCE SUM + sqrt.
  function mpi_allreduce_norm2(x, n) result(global_norm)
    integer, intent(in) :: n
    real(R8), intent(in) :: x(n)
    real(R8) :: global_norm
    real(R8) :: local_ss

    local_ss = dot_product(x, x)
#ifdef USE_MPI
    call mpi_allreduce_sum_r8(local_ss, global_norm)
    global_norm = sqrt(global_norm)
#else
    global_norm = sqrt(local_ss)
#endif
  end function mpi_allreduce_norm2


  !> MPI_REDUCE with MPI_SUM for a scalar real(R8) to root rank 0.
  subroutine mpi_reduce_sum_r8(val)
    real(R8), intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    if (mpi_rank_ == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, ice_comm, ierr)
    else
      call MPI_REDUCE(val, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, ice_comm, ierr)
    end if
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8


  !> MPI_REDUCE with MPI_SUM for an array of real(R8) to root rank 0.
  !> Result is only valid on root. More efficient than ALLREDUCE when
  !> only root needs the result (e.g. for residual output).
  subroutine mpi_reduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    if (mpi_rank_ == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, ice_comm, ierr)
    else
      call MPI_REDUCE(arr, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, ice_comm, ierr)
    end if
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8_array


  !> Broadcast a logical value from root (rank 0) to all ranks.
  subroutine mpi_bcast_logical(val)
    logical, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_LOGICAL, 0, ice_comm, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_logical


  !> Broadcast an integer value from root (rank 0) to all ranks.
  subroutine mpi_bcast_integer(val)
    integer, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_INTEGER, 0, ice_comm, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_integer


  !> Check MPI return code. Abort all ranks on error.
  subroutine check_mpi_error(ierr)
    integer, intent(in) :: ierr
#ifdef USE_MPI
    integer :: abort_ierr
    if (ierr /= 0) then
      write(*,'(A,I4,A,I6)') ' [RANK ', mpi_rank_, '] MPI error code: ', ierr
      call MPI_ABORT(MPI_COMM_WORLD, 1, abort_ierr)
    end if
#endif
  end subroutine check_mpi_error


  !> MPI-aware abort: prints message, then aborts all ranks.
  subroutine mpi_abort_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] ABORT: ', trim(message)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#else
    write(*,'(A,A)') ' ABORT: ', trim(message)
    error stop
#endif
  end subroutine mpi_abort_all


end module ICE_Mod_MPI
