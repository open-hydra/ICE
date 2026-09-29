!>@brief Inter-rank halo exchange for the condensed-phase primitive field of ICE.
!> Every rank holds the whole domain but updates only the blocks it owns. The ghost
!> cells of an owned block are filled from interior cells of other blocks: the two
!> source cells of a connection (101/201) and the donor cells of a chimera face (102).
!> When such a cell lies in a block owned by another rank, that rank sends it before
!> each ghost fill and it is written into the local copy of the remote block; the
!> ghost routines then read it exactly as in a serial run. This is how MOSE moves
!> chimera donor cells, applied here to connections as well.
!> In serial mode (USE_MPI not defined, or one rank), all routines are no-ops.
module ICE_Mod_GhostExchange
  use iso_fortran_env, only: R8 => real64, I4 => int32
  use ICE_Mod_MPI
  use ICE_Global_m, only: ngroups, ncond, guide

  implicit none
  private

  !> One interior cell that crosses a rank boundary (one particle group).
  type :: halo_cell_type
    integer :: b, p, i, j, k  !< Block, particle group and cell
    integer :: off = 0        !< Offset of its first value in the send/recv buffer
  end type halo_cell_type

  !> Contiguous run of cells exchanged with one remote rank (one message).
  type :: rank_group_type
    integer :: rank  = -1  !< Remote MPI rank
    integer :: first =  0  !< First cell in the send/recv list
    integer :: count =  0  !< Number of cells
    integer :: off   =  0  !< Buffer offset of the message
    integer :: len   =  0  !< Message length (reals)
  end type rank_group_type

  !> Communication schedule for the fine-grid BC list.
  type :: ghost_schedule_type
    integer :: n_send = 0, n_recv = 0
    type(halo_cell_type), allocatable :: send_list(:), recv_list(:)
    integer :: n_send_ranks = 0, n_recv_ranks = 0
    type(rank_group_type), allocatable :: send_groups(:), recv_groups(:)
    real(R8), allocatable :: send_buf(:), recv_buf(:)
    integer,  allocatable :: send_req(:), recv_req(:)  !< Persistent requests
    logical :: persistent_init = .false.
    logical :: built = .false.
  end type ghost_schedule_type

  integer, parameter :: halo_tag = 7101

  type(ghost_schedule_type), allocatable, target, private :: mg_sched(:)
  type(ghost_schedule_type), pointer,     public          :: ghost_sched => null()

  public :: build_ghost_schedule
  public :: build_local_bc_index
  public :: select_ghost_level
  public :: cleanup_ghost_schedule
  public :: exchange_ghost_prim
  public :: gather_prim_to_root
  public :: gather_diagnostic_to_root
  public :: scatter_prim_from_root
  public :: mpi_io_barrier

contains


  subroutine build_local_bc_index(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer :: i, n

    n = 0
    do i = 1, size(grid%bc)
      if (is_local_block(grid%bc(i)%b)) n = n + 1
    end do

    grid%n_local_bc = n
    if (allocated(grid%local_bc_idx)) deallocate(grid%local_bc_idx)
    allocate(grid%local_bc_idx(n))

    n = 0
    do i = 1, size(grid%bc)
      if (is_local_block(grid%bc(i)%b)) then
        n = n + 1
        grid%local_bc_idx(n) = i
      end if
    end do

  end subroutine build_local_bc_index


  !> Build the halo schedule from the BC list. Every rank scans the same list in the
  !> same order, so the cells a rank sends to a neighbour are listed in the order the
  !> neighbour expects them. Must be called after partition_blocks and Setup_BC.
  subroutine build_ghost_schedule(grid, level)
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_multigrid
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer, intent(in)               :: level  !< Multigrid level this schedule serves
#ifdef USE_MPI
    integer, allocatable :: send_per_rank(:), recv_per_rank(:), send_pos(:), recv_pos(:)
#endif

    if (.not. allocated(mg_sched)) allocate(mg_sched(max(1, obj_multigrid%MGL)))
    if (level < 1 .or. level > size(mg_sched)) &
      call mpi_abort_all('build_ghost_schedule: grid level outside 1..MG-levels')
    ghost_sched => mg_sched(level)

    if (mpi_size_ <= 1) then
      ghost_sched%built = .true.
      return
    end if

#ifdef USE_MPI
    block
      use mpi
      integer :: pass, n, c, fs, r, ierr
      integer :: ns, nbad
      integer, allocatable :: expect(:)

      allocate(send_per_rank(0:mpi_size_-1), recv_per_rank(0:mpi_size_-1), expect(0:mpi_size_-1))
      allocate(send_pos(0:mpi_size_-1), recv_pos(0:mpi_size_-1))

      ! Pass 1 counts cells per remote rank, pass 2 places them grouped by rank,
      ! keeping BC order inside each group.
      do pass = 1, 2
        send_per_rank = 0 ; recv_per_rank = 0
        do n = 1, size(grid%bc)
          select case (grid%bc(n)%type)
          case (101, 201)
            fs = grid%bc(n)%fs
            call add_cell(pass, grid%bc(n), grid%bc(n)%bs, grid%bc(n)%is, grid%bc(n)%js, grid%bc(n)%ks)
            call add_cell(pass, grid%bc(n), grid%bc(n)%bs, grid%bc(n)%is + guide(fs,1), &
                          grid%bc(n)%js + guide(fs,2), grid%bc(n)%ks + guide(fs,3))
          case (102)
            do c = 1, sum(grid%bc(n)%ni)
              call add_cell(pass, grid%bc(n), grid%bc(n)%donorID(c,1), grid%bc(n)%donorID(c,2), &
                            grid%bc(n)%donorID(c,3), grid%bc(n)%donorID(c,4))
            end do
          end select
        end do

        if (pass == 1) then
          ghost_sched%n_send = sum(send_per_rank)
          ghost_sched%n_recv = sum(recv_per_rank)
          allocate(ghost_sched%send_list(ghost_sched%n_send), ghost_sched%recv_list(ghost_sched%n_recv))
          send_pos(0) = 0 ; recv_pos(0) = 0
          do r = 1, mpi_size_-1
            send_pos(r) = send_pos(r-1) + send_per_rank(r-1)
            recv_pos(r) = recv_pos(r-1) + recv_per_rank(r-1)
          end do
        end if
      end do

      ! Each rank must receive exactly what its neighbours send it
      call MPI_ALLTOALL(send_per_rank, 1, MPI_INTEGER, expect, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
      nbad = count(expect /= recv_per_rank)
      if (nbad > 0) call mpi_abort_all('build_ghost_schedule: send and receive lists disagree')

      call build_rank_groups(ghost_sched%send_list, ghost_sched%n_send, send_per_rank, &
                             ghost_sched%send_groups, ghost_sched%n_send_ranks)
      call build_rank_groups(ghost_sched%recv_list, ghost_sched%n_recv, recv_per_rank, &
                             ghost_sched%recv_groups, ghost_sched%n_recv_ranks)

      allocate(ghost_sched%send_buf(max(1, sum(ghost_sched%send_groups(:)%len))))
      allocate(ghost_sched%recv_buf(max(1, sum(ghost_sched%recv_groups(:)%len))))

      call init_persistent_requests()

      ns = ghost_sched%n_send
      call MPI_ALLREDUCE(MPI_IN_PLACE, ns, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (mpi_is_root) write(*,'(A,I0,A)') '  MPI halo: ', ns, ' cells exchanged per ghost fill'

      deallocate(send_per_rank, recv_per_rank, expect, send_pos, recv_pos)
    end block
#endif

    ghost_sched%built = .true.

#ifdef USE_MPI
  contains

    !> Cell (b,i,j,k) is read by BC record bc. If its block and the ghost's block have
    !> different owners, the owner of b sends it to the owner of bc%b.
    subroutine add_cell(pass, bc, b, i, j, k)
      integer,           intent(in) :: pass, b, i, j, k
      type(ICE_bc_type), intent(in) :: bc
      integer :: owner_g, owner_s, pos

      owner_g = block_owner(bc%b)
      owner_s = block_owner(b)
      if (owner_g == owner_s) return

      if (owner_g == mpi_rank_) then
        recv_per_rank(owner_s) = recv_per_rank(owner_s) + 1
        if (pass == 2) then
          pos = recv_pos(owner_s) + recv_per_rank(owner_s)
          ghost_sched%recv_list(pos) = halo_cell_type(b, bc%p, i, j, k)
        end if
      else if (owner_s == mpi_rank_) then
        send_per_rank(owner_g) = send_per_rank(owner_g) + 1
        if (pass == 2) then
          pos = send_pos(owner_g) + send_per_rank(owner_g)
          ghost_sched%send_list(pos) = halo_cell_type(b, bc%p, i, j, k)
        end if
      end if
    end subroutine add_cell
#endif

  end subroutine build_ghost_schedule


  !> Point the halo routines at one level's schedule. Call whenever the solver moves
  !> between multigrid levels; build_ghost_schedule has to have run for that level.
  subroutine select_ghost_level(level)
    implicit none
    integer, intent(in) :: level

    if (.not. allocated(mg_sched)) return
    if (level < 1 .or. level > size(mg_sched)) &
      call mpi_abort_all('select_ghost_level: grid level outside 1..MG-levels')
    ghost_sched => mg_sched(level)

  end subroutine select_ghost_level


  !> Free persistent MPI requests, on every level that built any.
  subroutine cleanup_ghost_schedule()
    implicit none
#ifdef USE_MPI
    integer :: m

    if (.not. allocated(mg_sched)) return
    do m = 1, size(mg_sched)
      ghost_sched => mg_sched(m)
      call cleanup_persistent_requests()
    end do
    ghost_sched => mg_sched(1)
#endif
  end subroutine cleanup_ghost_schedule


  !> Bring the remote interior cells read by the ghost fill up to date. Call on every
  !> rank before compute_ghost, with ghost_sched pointing at this grid's level.
  subroutine exchange_ghost_prim(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid

    if (mpi_size_ <= 1) return
    if (size(grid%bc) == 0) return
#ifdef USE_MPI
    block
      use mpi
      integer :: n, nv, off, ierr
      integer, allocatable :: stats(:,:)

      if (ghost_sched%n_recv_ranks > 0) then
        call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, ierr)
        call check_mpi_error(ierr)
      end if

      do n = 1, ghost_sched%n_send
        associate (c => ghost_sched%send_list(n))
          nv = ncond(c%p) ; off = c%off
          ghost_sched%send_buf(off+1:off+nv) = grid%blk(c%b)%cond_phase(c%p)%prim(1:nv, c%i, c%j, c%k)
        end associate
      end do

      if (ghost_sched%n_send_ranks > 0) then
        call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%send_req, ierr)
        call check_mpi_error(ierr)
      end if

      allocate(stats(MPI_STATUS_SIZE, max(1, ghost_sched%n_send_ranks, ghost_sched%n_recv_ranks)))

      if (ghost_sched%n_recv_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, stats, ierr)
        call check_mpi_error(ierr)
      end if

      do n = 1, ghost_sched%n_recv
        associate (c => ghost_sched%recv_list(n))
          nv = ncond(c%p) ; off = c%off
          grid%blk(c%b)%cond_phase(c%p)%prim(1:nv, c%i, c%j, c%k) = ghost_sched%recv_buf(off+1:off+nv)
        end associate
      end do

      if (ghost_sched%n_send_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%send_req, stats, ierr)
        call check_mpi_error(ierr)
      end if
    end block
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
  !> Initialize one persistent receive and one persistent send per neighbour rank.
  subroutine init_persistent_requests()
    use mpi
    implicit none
    integer :: r, ierr

    allocate(ghost_sched%recv_req(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%send_req(max(ghost_sched%n_send_ranks, 1)))

    do r = 1, ghost_sched%n_recv_ranks
      associate (g => ghost_sched%recv_groups(r))
        call MPI_RECV_INIT(ghost_sched%recv_buf(g%off+1), g%len, MPI_DOUBLE_PRECISION, &
                           g%rank, halo_tag, MPI_COMM_WORLD, ghost_sched%recv_req(r), ierr)
        call check_mpi_error(ierr)
      end associate
    end do

    do r = 1, ghost_sched%n_send_ranks
      associate (g => ghost_sched%send_groups(r))
        call MPI_SEND_INIT(ghost_sched%send_buf(g%off+1), g%len, MPI_DOUBLE_PRECISION, &
                           g%rank, halo_tag, MPI_COMM_WORLD, ghost_sched%send_req(r), ierr)
        call check_mpi_error(ierr)
      end associate
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
      call MPI_REQUEST_FREE(ghost_sched%send_req(r), ierr)
    end do
    do r = 1, ghost_sched%n_recv_ranks
      call MPI_REQUEST_FREE(ghost_sched%recv_req(r), ierr)
    end do
    ghost_sched%persistent_init = .false.
  end subroutine cleanup_persistent_requests


  !> Split a rank-grouped cell list into one message per remote rank and assign each
  !> cell its buffer offset (cells carry ncond(p) values, which differs between groups).
  subroutine build_rank_groups(list, n, per_rank, groups, n_groups)
    implicit none
    type(halo_cell_type), intent(inout) :: list(:)
    integer, intent(in) :: n
    integer, intent(in) :: per_rank(0:)
    type(rank_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: r, g, c, first, off

    n_groups = count(per_rank > 0)
    allocate(groups(n_groups))

    g = 0 ; first = 1 ; off = 0
    do r = 0, size(per_rank)-1
      if (per_rank(r) == 0) cycle
      g = g + 1
      groups(g)%rank  = r
      groups(g)%first = first
      groups(g)%count = per_rank(r)
      groups(g)%off   = off
      do c = first, first + per_rank(r) - 1
        list(c)%off = off
        off = off + ncond(list(c)%p)
      end do
      groups(g)%len = off - groups(g)%off
      first = first + per_rank(r)
    end do
  end subroutine build_rank_groups
#endif


end module ICE_Mod_GhostExchange
