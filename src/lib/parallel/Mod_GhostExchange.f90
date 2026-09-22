!>@brief Ghost cell MPI communication for inter-block boundaries in ICE.
!> Builds a communication schedule from the BC connectivity array (type=101/201)
!> and provides persistent non-blocking exchange routines for the condensed-phase
!> primitive-variable field (prim). Messages are aggregated per remote rank.
!> In serial mode (USE_MPI not defined), all routines are no-ops.
module ICE_Mod_GhostExchange
  use iso_fortran_env, only: R8 => real64, I4 => int32
  use ICE_Mod_MPI
  use ICE_Global_m, only: ngroups, ncond, guide

  implicit none
  private

  !> Single ghost cell exchange entry (one BC face-cell, one particle group).
  type :: ghost_entry_type
    integer :: bc_idx         !< Index into grid%bc(:)
    integer :: local_block    !< Block on this rank
    integer :: remote_block   !< Block on the remote rank
    integer :: remote_rank    !< MPI rank of the remote block
    integer :: src_i, src_j, src_k    !< Source cell in remote block
    integer :: dst_i, dst_j, dst_k, dst_f  !< Dest cell and face in local block
    integer :: grp            !< Particle group index (bc%p)
    integer :: nvar           !< ncond(grp) — number of variables for this entry
  end type ghost_entry_type

  !> Per-rank aggregation: contiguous run of entries sharing the same remote rank.
  type :: rank_group_type
    integer :: rank   = -1  !< Remote MPI rank
    integer :: offset =  0  !< First entry index in send/recv list (1-based)
    integer :: count  =  0  !< Number of entries for this rank
  end type rank_group_type

  !> Communication schedule.
  type :: ghost_schedule_type
    integer :: n_send = 0, n_recv = 0
    type(ghost_entry_type), allocatable :: send_list(:)
    type(ghost_entry_type), allocatable :: recv_list(:)
    logical :: built = .false.

    ! Per-rank aggregation
    integer :: n_send_ranks = 0, n_recv_ranks = 0
    type(rank_group_type), allocatable :: send_groups(:)
    type(rank_group_type), allocatable :: recv_groups(:)

    ! Pre-allocated MPI buffers (one contiguous flat buffer, ncond_max values per entry)
    integer  :: entry_size = 0  !< maxval(ncond(1:ngroups))
    real(R8), allocatable :: send_buf(:)  !< (entry_size * n_send)
    real(R8), allocatable :: recv_buf(:)  !< (entry_size * n_recv)
    integer,  allocatable :: send_req(:)  !< (n_send_ranks)
    integer,  allocatable :: recv_req(:)  !< (n_recv_ranks)
    integer,  allocatable :: mpi_stat(:,:)!< (MPI_STATUS_SIZE, max(n_send_ranks,n_recv_ranks))

    ! Persistent MPI requests for prim exchange (initialized once, started each step)
    integer, allocatable :: prim_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable :: prim_recv_req_pers(:)  !< (n_recv_ranks)
    logical :: persistent_init = .false.
  end type ghost_schedule_type

  type(ghost_schedule_type), public :: ghost_sched

  public :: build_ghost_schedule
  public :: cleanup_ghost_schedule
  public :: exchange_ghost_prim
  public :: gather_prim_to_root
  public :: gather_diagnostic_to_root
  public :: scatter_prim_from_root
  public :: mpi_io_barrier

contains


  !> Build the communication schedule by scanning all BC entries of type 101/201 (connection).
  !> Entries are sorted by remote rank for aggregated MPI messaging.
  !> Must be called after partition_blocks and grid allocation.
  subroutine build_ghost_schedule(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer :: i, ns, nr, bm, bs, pm

    if (mpi_size_ <= 1) then
      ghost_sched%built = .true.
      return
    end if

#ifdef USE_MPI
    ghost_sched%entry_size = maxval(ncond(1:ngroups))

    ! Count send and recv entries
    ns = 0; nr = 0
    do i = 1, size(grid%bc)
      if (grid%bc(i)%type /= 101 .and. grid%bc(i)%type /= 201) cycle
      bm = grid%bc(i)%b
      bs = grid%bc(i)%bs
      if (.not. allocated(block_owner)) cycle
      if (is_local_block(bm) .and. (.not. is_local_block(bs))) nr = nr + 1
      if (is_local_block(bs) .and. (.not. is_local_block(bm))) ns = ns + 1
    end do

    ghost_sched%n_send = ns
    ghost_sched%n_recv = nr
    if (allocated(ghost_sched%send_list)) deallocate(ghost_sched%send_list)
    if (allocated(ghost_sched%recv_list)) deallocate(ghost_sched%recv_list)
    allocate(ghost_sched%send_list(ns))
    allocate(ghost_sched%recv_list(nr))

    ! Fill send and recv lists
    ns = 0; nr = 0
    do i = 1, size(grid%bc)
      if (grid%bc(i)%type /= 101 .and. grid%bc(i)%type /= 201) cycle
      bm = grid%bc(i)%b
      bs = grid%bc(i)%bs
      pm = grid%bc(i)%p

      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
        ghost_sched%recv_list(nr)%bc_idx       = i
        ghost_sched%recv_list(nr)%local_block  = bm
        ghost_sched%recv_list(nr)%remote_block = bs
        ghost_sched%recv_list(nr)%remote_rank  = block_owner(bs)
        ghost_sched%recv_list(nr)%src_i        = grid%bc(i)%is
        ghost_sched%recv_list(nr)%src_j        = grid%bc(i)%js
        ghost_sched%recv_list(nr)%src_k        = grid%bc(i)%ks
        ghost_sched%recv_list(nr)%dst_i        = grid%bc(i)%i
        ghost_sched%recv_list(nr)%dst_j        = grid%bc(i)%j
        ghost_sched%recv_list(nr)%dst_k        = grid%bc(i)%k
        ghost_sched%recv_list(nr)%dst_f        = grid%bc(i)%f
        ghost_sched%recv_list(nr)%grp          = pm
        ghost_sched%recv_list(nr)%nvar         = ncond(pm)
      end if

      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
        ghost_sched%send_list(ns)%bc_idx       = i
        ghost_sched%send_list(ns)%local_block  = bs
        ghost_sched%send_list(ns)%remote_block = bm
        ghost_sched%send_list(ns)%remote_rank  = block_owner(bm)
        ghost_sched%send_list(ns)%src_i        = grid%bc(i)%is
        ghost_sched%send_list(ns)%src_j        = grid%bc(i)%js
        ghost_sched%send_list(ns)%src_k        = grid%bc(i)%ks
        ghost_sched%send_list(ns)%dst_i        = grid%bc(i)%i
        ghost_sched%send_list(ns)%dst_j        = grid%bc(i)%j
        ghost_sched%send_list(ns)%dst_k        = grid%bc(i)%k
        ghost_sched%send_list(ns)%dst_f        = grid%bc(i)%f
        ghost_sched%send_list(ns)%grp          = pm
        ghost_sched%send_list(ns)%nvar         = ncond(pm)
      end if
    end do

    ! Sort by remote rank and build aggregation groups
    call sort_entries_by_rank(ghost_sched%send_list, ghost_sched%n_send)
    call sort_entries_by_rank(ghost_sched%recv_list, ghost_sched%n_recv)
    call build_rank_groups(ghost_sched%send_list, ghost_sched%n_send, &
                           ghost_sched%send_groups, ghost_sched%n_send_ranks)
    call build_rank_groups(ghost_sched%recv_list, ghost_sched%n_recv, &
                           ghost_sched%recv_groups, ghost_sched%n_recv_ranks)

    ! Pre-allocate flat MPI buffers (reused every RK stage)
    if (ghost_sched%n_send > 0 .or. ghost_sched%n_recv > 0) then
      block
        use mpi, only: MPI_STATUS_SIZE
        integer :: es, max_ranks
        es        = ghost_sched%entry_size
        max_ranks = max(ghost_sched%n_send_ranks, ghost_sched%n_recv_ranks, 1)
        allocate(ghost_sched%send_buf (es * max(ghost_sched%n_send, 1)))
        allocate(ghost_sched%recv_buf (es * max(ghost_sched%n_recv, 1)))
        allocate(ghost_sched%send_req (max(ghost_sched%n_send_ranks, 1)))
        allocate(ghost_sched%recv_req (max(ghost_sched%n_recv_ranks, 1)))
        allocate(ghost_sched%mpi_stat (MPI_STATUS_SIZE, max_ranks))
      end block
    end if

    call init_persistent_requests()

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A)') ' ICE ghost schedule: ', ghost_sched%n_send_ranks, &
        ' send ranks, ', ghost_sched%n_recv_ranks, ' recv ranks'
    end if
#endif

    ghost_sched%built = .true.
  end subroutine build_ghost_schedule


  !> Free persistent MPI requests and deallocate schedule arrays.
  subroutine cleanup_ghost_schedule()
    implicit none
#ifdef USE_MPI
    call cleanup_persistent_requests()
#endif
  end subroutine cleanup_ghost_schedule


  !> Blocking exchange of the condensed-phase prim field across MPI ranks.
  !> Covers all particle groups. Call once per ghost-cell update (each RK stage).
  subroutine exchange_ghost_prim(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_prim_begin(grid)
    call exchange_prim_end(grid)
#endif
  end subroutine exchange_ghost_prim


  !> Gather all blocks' prim interior data to rank 0 for I/O.
  !> Uses non-blocking ISEND/IRECV + WAITALL. All ranks must call this (collective).
  subroutine gather_prim_to_root(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer :: b, p, ni, nj, nk

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, offset, total_buf
      integer :: blk_size
      integer, allocatable :: reqs(:), blk_offset(:,:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Per-block size: sum of ncond(p)*ni*nj*nk over all groups
      allocate(blk_ncells(grid%nb), blk_offset(grid%nb, ngroups))

      total_buf = 0
      do b = 1, grid%nb
        ni = grid%blk(b)%dim(1)
        nj = grid%blk(b)%dim(2)
        nk = grid%blk(b)%dim(3)
        blk_ncells(b) = 0
        do p = 1, ngroups
          blk_offset(b, p) = total_buf
          blk_size = ncond(p) * ni * nj * nk
          if (mpi_is_root) then
            if (.not. is_local_block(b)) total_buf = total_buf + blk_size
          else
            if (is_local_block(b)) total_buf = total_buf + blk_size
          end if
          blk_ncells(b) = blk_ncells(b) + blk_size
        end do
      end do

      allocate(buf(max(total_buf, 1)))
      allocate(reqs(grid%nb * ngroups))
      allocate(stats(MPI_STATUS_SIZE, grid%nb * ngroups))

      nreq = 0
      if (mpi_is_root) then
        ! Root: post receives for remote blocks
        do b = 1, grid%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b,1)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack and send local blocks
        do b = 1, grid%nb
          if (is_local_block(b)) then
            ni = grid%blk(b)%dim(1)
            nj = grid%blk(b)%dim(2)
            nk = grid%blk(b)%dim(3)
            offset = blk_offset(b, 1)
            do p = 1, ngroups
              blk_size = ncond(p) * ni * nj * nk
              buf(offset+1 : offset+blk_size) = &
                reshape(grid%blk(b)%cond_phase(p)%prim(1:ncond(p), 1:ni, 1:nj, 1:nk), [blk_size])
              offset = offset + blk_size
            end do
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b,1)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data
      if (mpi_is_root) then
        do b = 1, grid%nb
          if (.not. is_local_block(b)) then
            ni = grid%blk(b)%dim(1)
            nj = grid%blk(b)%dim(2)
            nk = grid%blk(b)%dim(3)
            offset = blk_offset(b, 1)
            do p = 1, ngroups
              blk_size = ncond(p) * ni * nj * nk
              grid%blk(b)%cond_phase(p)%prim(1:ncond(p), 1:ni, 1:nj, 1:nk) = &
                reshape(buf(offset+1 : offset+blk_size), [ncond(p), ni, nj, nk])
              offset = offset + blk_size
            end do
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_ncells, blk_offset)
    end block
#endif
  end subroutine gather_prim_to_root


  !> Gather diagnostic fields (residual, dt) from all ranks to root for writing.
  !> All ranks must call this (collective).
  subroutine gather_diagnostic_to_root(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer :: b, p, ni, nj, nk

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, offset, total_buf
      integer :: blk_size_r, blk_size_dt
      integer, allocatable :: reqs(:), blk_offset(:,:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      allocate(blk_ncells(grid%nb), blk_offset(grid%nb, ngroups+1))

      total_buf = 0
      do b = 1, grid%nb
        ni = grid%blk(b)%dim(1)
        nj = grid%blk(b)%dim(2)
        nk = grid%blk(b)%dim(3)
        blk_ncells(b) = 0
        do p = 1, ngroups
          blk_offset(b, p) = total_buf
          blk_size_r = ncond(p) * ni * nj * nk
          if (mpi_is_root) then
            if (.not. is_local_block(b)) total_buf = total_buf + blk_size_r
          else
            if (is_local_block(b)) total_buf = total_buf + blk_size_r
          end if
          blk_ncells(b) = blk_ncells(b) + blk_size_r
        end do
        blk_offset(b, ngroups+1) = total_buf  ! dt offset
        blk_size_dt = ngroups * ni * nj * nk
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf = total_buf + blk_size_dt
        else
          if (is_local_block(b)) total_buf = total_buf + blk_size_dt
        end if
        blk_ncells(b) = blk_ncells(b) + blk_size_dt
      end do

      allocate(buf(max(total_buf, 1)))
      allocate(reqs(grid%nb))
      allocate(stats(MPI_STATUS_SIZE, grid%nb))

      nreq = 0
      if (mpi_is_root) then
        do b = 1, grid%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b,1)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        do b = 1, grid%nb
          if (is_local_block(b)) then
            ni = grid%blk(b)%dim(1)
            nj = grid%blk(b)%dim(2)
            nk = grid%blk(b)%dim(3)
            offset = blk_offset(b, 1)
            do p = 1, ngroups
              blk_size_r = ncond(p) * ni * nj * nk
              buf(offset+1 : offset+blk_size_r) = &
                reshape(grid%blk(b)%cond_phase(p)%residual(1:ncond(p), 1:ni, 1:nj, 1:nk), [blk_size_r])
              offset = offset + blk_size_r
            end do
            do p = 1, ngroups
              blk_size_dt = ni * nj * nk
              buf(offset+1 : offset+blk_size_dt) = &
                reshape(grid%blk(b)%cond_phase(p)%dt(1:ni, 1:nj, 1:nk), [blk_size_dt])
              offset = offset + blk_size_dt
            end do
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b,1)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      if (mpi_is_root) then
        do b = 1, grid%nb
          if (.not. is_local_block(b)) then
            ni = grid%blk(b)%dim(1)
            nj = grid%blk(b)%dim(2)
            nk = grid%blk(b)%dim(3)
            offset = blk_offset(b, 1)
            do p = 1, ngroups
              blk_size_r = ncond(p) * ni * nj * nk
              grid%blk(b)%cond_phase(p)%residual(1:ncond(p), 1:ni, 1:nj, 1:nk) = &
                reshape(buf(offset+1 : offset+blk_size_r), [ncond(p), ni, nj, nk])
              offset = offset + blk_size_r
            end do
            do p = 1, ngroups
              blk_size_dt = ni * nj * nk
              grid%blk(b)%cond_phase(p)%dt(1:ni, 1:nj, 1:nk) = &
                reshape(buf(offset+1 : offset+blk_size_dt), [ni, nj, nk])
              offset = offset + blk_size_dt
            end do
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_ncells, blk_offset)
    end block
#endif
  end subroutine gather_diagnostic_to_root


  !> Scatter all blocks' prim data from rank 0 to owning ranks after reading.
  !> Uses non-blocking ISEND/IRECV + WAITALL. All ranks must call this (collective).
  subroutine scatter_prim_from_root(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer :: b, p, ni, nj, nk

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, offset, total_buf
      integer :: blk_size
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      allocate(blk_ncells(grid%nb), blk_offset(grid%nb))

      total_buf = 0
      do b = 1, grid%nb
        ni = grid%blk(b)%dim(1)
        nj = grid%blk(b)%dim(2)
        nk = grid%blk(b)%dim(3)
        blk_ncells(b) = 0
        do p = 1, ngroups
          blk_ncells(b) = blk_ncells(b) + ncond(p) * ni * nj * nk
        end do
        blk_offset(b) = total_buf
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf = total_buf + blk_ncells(b)
        else
          if (is_local_block(b)) total_buf = total_buf + blk_ncells(b)
        end if
      end do

      allocate(buf(max(total_buf, 1)))
      allocate(reqs(grid%nb))
      allocate(stats(MPI_STATUS_SIZE, grid%nb))

      nreq = 0
      if (mpi_is_root) then
        ! Root: pack and send to non-local blocks' owners
        do b = 1, grid%nb
          if (.not. is_local_block(b)) then
            ni = grid%blk(b)%dim(1)
            nj = grid%blk(b)%dim(2)
            nk = grid%blk(b)%dim(3)
            offset = blk_offset(b)
            do p = 1, ngroups
              blk_size = ncond(p) * ni * nj * nk
              buf(offset+1 : offset+blk_size) = &
                reshape(grid%blk(b)%cond_phase(p)%prim(1:ncond(p), 1:ni, 1:nj, 1:nk), [blk_size])
              offset = offset + blk_size
            end do
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: post receives for local blocks
        do b = 1, grid%nb
          if (is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Non-root: unpack received data into prim
      if (.not. mpi_is_root) then
        do b = 1, grid%nb
          if (is_local_block(b)) then
            ni = grid%blk(b)%dim(1)
            nj = grid%blk(b)%dim(2)
            nk = grid%blk(b)%dim(3)
            offset = blk_offset(b)
            do p = 1, ngroups
              blk_size = ncond(p) * ni * nj * nk
              grid%blk(b)%cond_phase(p)%prim(1:ncond(p), 1:ni, 1:nj, 1:nk) = &
                reshape(buf(offset+1 : offset+blk_size), [ncond(p), ni, nj, nk])
              offset = offset + blk_size
            end do
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_ncells, blk_offset)
    end block
#endif
  end subroutine scatter_prim_from_root


  !> MPI barrier for synchronizing ranks after I/O operations.
  subroutine mpi_io_barrier()
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end block
#endif
  end subroutine mpi_io_barrier


#ifdef USE_MPI
  !> Initialize persistent MPI requests for prim exchange.
  subroutine init_persistent_requests()
    use mpi
    implicit none
    integer :: r, ierr, tag, buf_pos, es

    es = ghost_sched%entry_size
    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%prim_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%prim_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Persistent receives
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * es + 1
      tag = ghost_sched%recv_groups(r)%rank
      call MPI_RECV_INIT(ghost_sched%recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * es, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%prim_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Persistent sends
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * es + 1
      tag = mpi_rank_
      call MPI_SEND_INIT(ghost_sched%send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * es, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%prim_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_init = .true.
  end subroutine init_persistent_requests


  !> Free all persistent MPI requests.
  subroutine cleanup_persistent_requests()
    use mpi
    implicit none
    integer :: ierr, r

    if (.not. ghost_sched%persistent_init) return
    do r = 1, ghost_sched%n_send_ranks
      call MPI_REQUEST_FREE(ghost_sched%prim_send_req_pers(r), ierr)
    end do
    do r = 1, ghost_sched%n_recv_ranks
      call MPI_REQUEST_FREE(ghost_sched%prim_recv_req_pers(r), ierr)
    end do
    ghost_sched%persistent_init = .false.
  end subroutine cleanup_persistent_requests


  !> Post receives, pack and send prim field.
  subroutine exchange_prim_begin(grid)
    use ICE_Advanced_Types_m
    use mpi
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer :: i, r, ierr, buf_pos, es
    integer :: bs, Is, Js, Ks, pm, nv

    es = ghost_sched%entry_size

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    ! Start persistent receives
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%prim_recv_req_pers, ierr)
      call check_mpi_error(ierr)
    end if

    ! Pack send buffer: for each send entry copy prim(:, is, js, ks) from local source block
    do i = 1, ghost_sched%n_send
      bs  = ghost_sched%send_list(i)%local_block
      Is  = ghost_sched%send_list(i)%src_i
      Js  = ghost_sched%send_list(i)%src_j
      Ks  = ghost_sched%send_list(i)%src_k
      pm  = ghost_sched%send_list(i)%grp
      nv  = ghost_sched%send_list(i)%nvar
      buf_pos = (i - 1) * es
      ghost_sched%send_buf(buf_pos+1 : buf_pos+nv) = &
        grid%blk(bs)%cond_phase(pm)%prim(1:nv, Is, Js, Ks)
    end do

    ! Start persistent sends
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%prim_send_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_prim_begin


  !> Wait for receives and unpack prim into ghost cells; then wait for sends.
  subroutine exchange_prim_end(grid)
    use ICE_Advanced_Types_m
    use mpi
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer :: i, ierr, buf_pos, es
    integer :: bm, Im, Jm, Km, Fm, Ig, Jg, Kg, pm, nv

    es = ghost_sched%entry_size

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    ! Wait for all receives
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%prim_recv_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack recv buffer into ghost cells
    do i = 1, ghost_sched%n_recv
      bm  = ghost_sched%recv_list(i)%local_block
      Im  = ghost_sched%recv_list(i)%dst_i
      Jm  = ghost_sched%recv_list(i)%dst_j
      Km  = ghost_sched%recv_list(i)%dst_k
      Fm  = ghost_sched%recv_list(i)%dst_f
      pm  = ghost_sched%recv_list(i)%grp
      nv  = ghost_sched%recv_list(i)%nvar
      ! Ghost cell is one step inward from the boundary face
      Ig  = Im - guide(Fm, 1)
      Jg  = Jm - guide(Fm, 2)
      Kg  = Km - guide(Fm, 3)
      buf_pos = (i - 1) * es
      grid%blk(bm)%cond_phase(pm)%prim(1:nv, Ig, Jg, Kg) = &
        ghost_sched%recv_buf(buf_pos+1 : buf_pos+nv)
    end do

    ! Wait for sends to complete
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%prim_send_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_prim_end


  !> Sort ghost entries by remote_rank using insertion sort (lists are short).
  subroutine sort_entries_by_rank(list, n)
    implicit none
    type(ghost_entry_type), intent(inout) :: list(:)
    integer, intent(in) :: n
    integer :: i, j
    type(ghost_entry_type) :: tmp

    do i = 2, n
      tmp = list(i)
      j = i - 1
      do while (j >= 1 .and. list(j)%remote_rank > tmp%remote_rank)
        list(j+1) = list(j)
        j = j - 1
      end do
      list(j+1) = tmp
    end do
  end subroutine sort_entries_by_rank


  !> Build rank groups from a sorted entry list.
  subroutine build_rank_groups(list, n, groups, n_groups)
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(rank_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ng = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%rank   = list(1)%remote_rank
    groups(1)%offset = 1
    groups(1)%count  = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) then
        ng = ng + 1
        groups(ng)%rank   = list(i)%remote_rank
        groups(ng)%offset = i
        groups(ng)%count  = 1
      else
        groups(ng)%count = groups(ng)%count + 1
      end if
    end do
  end subroutine build_rank_groups
#endif


end module ICE_Mod_GhostExchange
