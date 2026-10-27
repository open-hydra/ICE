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
  use ICE_Mod_Timers, only: timer_region_begin, timer_region_end, TR_PACK, TR_WAIT, TR_UNPACK
  use ICE_Mod_ThreadGroups, only: n_tgroups

  implicit none
  private

  !> One interior cell that crosses a rank boundary (one particle group).
  type :: halo_cell_type
    integer :: b, p, i, j, k  !< Block, particle group and cell
    integer :: off = 0        !< Offset of its first value in the send/recv buffer
  end type halo_cell_type

  !> The cells of one block a rank run has already listed (unique_cells).
  type :: mask_type
    logical, allocatable :: seen(:,:,:)
  end type mask_type

  !> Contiguous run of cells exchanged with one remote rank (one message).
  type :: rank_group_type
    integer :: rank  = -1  !< Remote MPI rank
    integer :: first =  0  !< First cell in the send/recv list
    integer :: count =  0  !< Number of cells
    integer :: off   =  0  !< Buffer offset of the message
    integer :: len   =  0  !< Message length (reals)
  end type rank_group_type

  !> Communication schedule of one family on one grid level.
  type :: ghost_schedule_type
    integer :: group = 1                  !< the family, and the message tag offset
    integer :: n_send = 0, n_recv = 0
    type(halo_cell_type), allocatable :: send_list(:), recv_list(:)
    integer :: n_send_ranks = 0, n_recv_ranks = 0
    type(rank_group_type), allocatable :: send_groups(:), recv_groups(:)
    real(R8), allocatable :: send_buf(:), recv_buf(:)
    integer,  allocatable :: send_req(:), recv_req(:)  !< Persistent requests
    integer,  allocatable :: stats(:,:)                !< Statuses of their waits (no MPI_STATUSES_IGNORE)
    logical :: persistent_init = .false.
    !> The connection faces whose two records lie in blocks of different ranks: the rank of the owner
    !> record sends the face's flux (ncond values) to the rank of the partner (exchange_conn_flux).
    !> Face numbers (bc%cf), grouped by remote rank, table order inside a group.
    integer :: cf_n_send = 0, cf_n_recv = 0
    integer, allocatable :: cf_send(:), cf_recv(:)
    integer :: cf_n_send_ranks = 0, cf_n_recv_ranks = 0
    type(rank_group_type), allocatable :: cf_send_groups(:), cf_recv_groups(:)
    real(R8), allocatable :: cf_send_buf(:), cf_recv_buf(:)
    integer,  allocatable :: cf_send_req(:), cf_recv_req(:), cf_stats(:,:)
    logical :: cf_persistent_init = .false.
    !> The same faces the other way, before the owner pass (exchange_conn_beta): the partner's rank
    !> sends its cell's detector weight (one value per face, cf_recv's faces) to the owner's rank,
    !> which receives it for cf_send's faces.
    real(R8), allocatable :: cb_send_buf(:), cb_recv_buf(:)
    integer,  allocatable :: cb_send_req(:), cb_recv_req(:)
    logical :: built = .false.
  end type ghost_schedule_type

  !> Message tags, within 32767, the least MPI_TAG_UB the MPI standard allows. The
  !> halo exchange of family p uses halo_tag + p, its connection-face fluxes
  !> conn_tag + p and the partner cells' weights cbeta_tag + p. The gathers and the scatter of the
  !> output use io_tag for every block: each receive names its source rank, and both
  !> sides post the messages of a pair of ranks in increasing block order, which MPI
  !> matches in the order they were started (MPI 4.1, Sec. 3.7.4), so the k-th
  !> receive from a rank takes the k-th block that rank sends.
  integer, parameter :: halo_tag = 7101
  integer, parameter :: io_tag   = 7100
  integer, parameter :: conn_tag = 7600
  integer, parameter :: cbeta_tag = 7700

  !> One schedule per (level, family): a stage of family p fills only p's ghosts,
  !> so only p's cells cross. ghost_sched points at the one being used.
  type(ghost_schedule_type), allocatable, target, private :: mg_sched(:,:)
  type(ghost_schedule_type), pointer,     public          :: ghost_sched => null()
  integer, private :: level_cur = 1

  public :: build_ghost_schedule
  public :: build_local_bc_index
  public :: report_bc_mix
  public :: select_ghost_level
  public :: cleanup_ghost_schedule
  public :: exchange_ghost_prim
  public :: exchange_conn_flux
  public :: exchange_conn_beta
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

    ! The same records grouped by family, in table order inside each group: a
    ! stage of family p fills p's ghosts only, and reads its range of this list.
    if (allocated(grid%local_bc_grp)) deallocate(grid%local_bc_grp)
    if (allocated(grid%grp_first))    deallocate(grid%grp_first)
    allocate(grid%local_bc_grp(max(n, 1)), grid%grp_first(ngroups + 1))
    block
      integer, allocatable :: pos(:)
      integer :: p
      allocate(pos(ngroups))
      pos = 0
      do i = 1, n
        p = grid%bc(grid%local_bc_idx(i))%p
        pos(p) = pos(p) + 1
      end do
      grid%grp_first(1) = 1
      do p = 1, ngroups
        grid%grp_first(p+1) = grid%grp_first(p) + pos(p)
      end do
      pos = grid%grp_first(1:ngroups)
      do i = 1, n
        p = grid%bc(grid%local_bc_idx(i))%p
        grid%local_bc_grp(pos(p)) = grid%local_bc_idx(i)
        pos(p) = pos(p) + 1
      end do
      deallocate(pos)
    end block

    ! One list of boundary cells per family, each cell with its records in
    ! face order. compute_bound then accumulates a cell's boundary
    ! fluxes from one thread in one fixed order: an edge or corner cell carries
    ! two or three records, and (a+b)+c is not (a+c)+b.
    if (allocated(grid%bcells)) deallocate(grid%bcells)
    allocate(grid%bcells(ngroups))
    do i = 1, ngroups
      call collect_boundary_cells(grid, i, grid%bcells(i))
    end do

    ! The connection faces whose flux this rank computes, one per face (compute_bound)
    if (.not. allocated(grid%conn)) allocate(grid%conn(ngroups))
    do i = 1, ngroups
      call collect_owner_records(grid, i, grid%conn(i))
    end do

    ! With thread groups, the lists again by the group owning each record's block
    if (n_tgroups > 1) call build_group_order(grid)

  end subroutine build_local_bc_index


  !> The owner records of family p's connection faces in this rank's blocks, in table order, and
  !> the flux and the partner cell's weight of every face of the family: a face keeps its number on
  !> every rank, and a rank holds the faces it computes and those it receives (exchange_conn_flux,
  !> exchange_conn_beta).
  subroutine collect_owner_records(grid, p, conn)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type),     intent(in)    :: grid
    integer,                   intent(in)    :: p
    type(ICE_conn_faces_type), intent(inout) :: conn
    integer :: ii, n, k

    if (allocated(conn%F))   deallocate(conn%F)
    if (allocated(conn%bp))  deallocate(conn%bp)
    if (allocated(conn%own)) deallocate(conn%own)
    allocate(conn%F(ncond(p), max(conn%n, 1)), conn%bp(max(conn%n, 1)))
    conn%F = 0._R8
    conn%bp = 1._R8
    k = 0
    do ii = grid%grp_first(p), grid%grp_first(p+1) - 1
      n = grid%local_bc_grp(ii)
      if ((grid%bc(n)%type == 101 .or. grid%bc(n)%type == 201) .and. grid%bc(n)%owner) k = k + 1
    end do
    conn%nown = k
    allocate(conn%own(max(k, 1)))
    k = 0
    do ii = grid%grp_first(p), grid%grp_first(p+1) - 1
      n = grid%local_bc_grp(ii)
      if ((grid%bc(n)%type == 101 .or. grid%bc(n)%type == 201) .and. grid%bc(n)%owner) then
        k = k + 1
        conn%own(k) = n
      end if
    end do
  end subroutine collect_owner_records


  !> The records of local_bc_grp, the boundary cells and the connection owner records of
  !> every family, re-ordered by the thread group owning their block (stable: table order
  !> inside a group), so a group's threads sweep only their own blocks' ghosts and
  !> boundary fluxes.
  subroutine build_group_order(grid)
    use ICE_Global_m,         only: ngroups
    use ICE_Advanced_Types_m
    use ICE_Mod_ThreadGroups, only: block_group
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer :: p, g, ii, c, pos
    integer, allocatable :: own(:)

    if (allocated(grid%tg_bc))       deallocate(grid%tg_bc)
    if (allocated(grid%tg_bc_first)) deallocate(grid%tg_bc_first)
    allocate(grid%tg_bc(max(grid%n_local_bc, 1)), grid%tg_bc_first(0:n_tgroups, ngroups))
    do p = 1, ngroups
      pos = grid%grp_first(p)
      do g = 0, n_tgroups - 1
        grid%tg_bc_first(g, p) = pos
        do ii = grid%grp_first(p), grid%grp_first(p+1) - 1
          if (block_group(grid%bc(grid%local_bc_grp(ii))%b) /= g) cycle
          grid%tg_bc(pos) = ii
          pos = pos + 1
        end do
      end do
      grid%tg_bc_first(n_tgroups, p) = pos
    end do

    do p = 1, ngroups
      associate (cells => grid%bcells(p))
      if (allocated(cells%tg_cell))  deallocate(cells%tg_cell)
      if (allocated(cells%tg_first)) deallocate(cells%tg_first)
      allocate(cells%tg_cell(max(cells%n, 1)), cells%tg_first(0:n_tgroups))
      pos = 1
      do g = 0, n_tgroups - 1
        cells%tg_first(g) = pos
        do c = 1, cells%n
          if (block_group(grid%bc(cells%rec(cells%first(c)))%b) /= g) cycle
          cells%tg_cell(pos) = c
          pos = pos + 1
        end do
      end do
      cells%tg_first(n_tgroups) = pos
      end associate
    end do

    do p = 1, ngroups
      associate (cn => grid%conn(p))
      if (allocated(cn%tg_first)) deallocate(cn%tg_first)
      allocate(cn%tg_first(0:n_tgroups), own(max(cn%nown, 1)))
      pos = 1
      do g = 0, n_tgroups - 1
        cn%tg_first(g) = pos
        do ii = 1, cn%nown
          if (block_group(grid%bc(cn%own(ii))%b) /= g) cycle
          own(pos) = cn%own(ii)
          pos = pos + 1
        end do
      end do
      cn%tg_first(n_tgroups) = pos
      cn%own(1:cn%nown) = own(1:cn%nown)
      deallocate(own)
      end associate
    end do
  end subroutine build_group_order


  !> The boundary cells of family p on this rank, block by block in cell order,
  !> each with its records in face order. Type-0 records (no condition) are
  !> left out. Linear in the records and the cells of the owned blocks.
  subroutine collect_boundary_cells(grid, p, cells)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type),   intent(in)    :: grid
    integer,                 intent(in)    :: p
    type(ICE_bc_cells_type), intent(inout) :: cells
    integer, allocatable :: head(:,:,:), tail(:,:,:), next(:)
    integer :: b, nn, n, i, j, k, c, r, ncell, nrec, m, l

    ! Pass 1 counts, pass 2 fills; the chains are rebuilt per block each pass.
    nrec = 0 ; ncell = 0
    allocate(next(max(grid%n_local_bc, 1)))
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      allocate(head(grid%blk(b)%dim(1), grid%blk(b)%dim(2), grid%blk(b)%dim(3)))
      head = 0
      do nn = 1, grid%n_local_bc
        n = grid%local_bc_idx(nn)
        if (grid%bc(n)%b /= b .or. grid%bc(n)%p /= p .or. grid%bc(n)%type == 0) cycle
        nrec = nrec + 1
        if (head(grid%bc(n)%i, grid%bc(n)%j, grid%bc(n)%k) == 0) ncell = ncell + 1
        head(grid%bc(n)%i, grid%bc(n)%j, grid%bc(n)%k) = 1
      end do
      deallocate(head)
    end do

    cells%n = ncell
    if (allocated(cells%first)) deallocate(cells%first)
    if (allocated(cells%rec))   deallocate(cells%rec)
    allocate(cells%first(ncell + 1), cells%rec(max(nrec, 1)))

    c = 0 ; r = 0
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      allocate(head(grid%blk(b)%dim(1), grid%blk(b)%dim(2), grid%blk(b)%dim(3)))
      allocate(tail(grid%blk(b)%dim(1), grid%blk(b)%dim(2), grid%blk(b)%dim(3)))
      head = 0 ; tail = 0
      do nn = 1, grid%n_local_bc
        n = grid%local_bc_idx(nn)
        if (grid%bc(n)%b /= b .or. grid%bc(n)%p /= p .or. grid%bc(n)%type == 0) cycle
        i = grid%bc(n)%i ; j = grid%bc(n)%j ; k = grid%bc(n)%k
        next(nn) = 0
        if (head(i,j,k) == 0) then
          head(i,j,k) = nn
        else
          next(tail(i,j,k)) = nn
        end if
        tail(i,j,k) = nn
      end do
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)
        if (head(i,j,k) == 0) cycle
        c = c + 1
        cells%first(c) = r + 1
        nn = head(i,j,k)
        do while (nn /= 0)
          r = r + 1
          cells%rec(r) = grid%local_bc_idx(nn)
          nn = next(nn)
        end do
        ! The cell's records by face, stably (at most six): a file that lists them in
        ! another order gives the sum of a file that lists them by face
        do m = cells%first(c) + 1, r
          n = cells%rec(m)
          l = m - 1
          do while (l >= cells%first(c))
            if (grid%bc(cells%rec(l))%f <= grid%bc(n)%f) exit
            cells%rec(l+1) = cells%rec(l)
            l = l - 1
          end do
          cells%rec(l+1) = n
        end do
      end do ; end do ; end do
      deallocate(head, tail)
    end do
    cells%first(ncell + 1) = r + 1
    deallocate(next)
  end subroutine collect_boundary_cells


  !> Build the halo schedules of one grid level, one per family, from the BC list.
  !> Every rank scans the same list in the same order, so the cells a rank sends to
  !> a neighbour are listed in the order the neighbour expects them. Must be called
  !> after partition_blocks and Setup_BC.
  subroutine build_ghost_schedule(grid, level)
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_multigrid
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer, intent(in)               :: level  !< Multigrid level this schedule serves
    integer :: p
#ifdef USE_MPI
    integer, allocatable :: send_per_rank(:), recv_per_rank(:), send_pos(:), recv_pos(:)
    integer :: ns_total, ndup, nf_total
#endif

    if (.not. allocated(mg_sched)) allocate(mg_sched(max(1, obj_multigrid%MGL), ngroups))
    if (level < 1 .or. level > size(mg_sched, 1)) &
      call mpi_abort_all('build_ghost_schedule: grid level outside 1..MG-levels')
    level_cur = level
    ghost_sched => mg_sched(level, 1)

    if (mpi_size_ <= 1) then
      do p = 1, ngroups
        mg_sched(level, p)%group = p
        mg_sched(level, p)%built = .true.
      end do
      return
    end if

#ifdef USE_MPI
    ns_total = 0 ; ndup = 0 ; nf_total = 0
    do p = 1, ngroups
      ghost_sched => mg_sched(level, p)
      ghost_sched%group = p
      block
        use mpi
        integer :: pass, n, c, fs, r, ierr
        integer :: ns, nbad
        integer, allocatable :: expect(:)

        allocate(send_per_rank(0:mpi_size_-1), recv_per_rank(0:mpi_size_-1), expect(0:mpi_size_-1))
        allocate(send_pos(0:mpi_size_-1), recv_pos(0:mpi_size_-1))

        ! Pass 1 counts cells per remote rank, pass 2 places them grouped by rank,
        ! keeping BC order inside each group. Only this family's records.
        do pass = 1, 2
          send_per_rank = 0 ; recv_per_rank = 0
          do n = 1, size(grid%bc)
            if (grid%bc(n)%p /= p) cycle
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

        ! A cell read by several records (a block edge seen from two faces, a donor
        ! of several receivers) crosses once: keep the first occurrence in each
        ! rank's run. Both ends of a message hold the same run in the same order,
        ! so the same rule on both keeps them matched.
        ndup = ndup + ghost_sched%n_send
        call unique_cells(ghost_sched%send_list, send_per_rank, ghost_sched%n_send)
        ndup = ndup - ghost_sched%n_send
        call unique_cells(ghost_sched%recv_list, recv_per_rank, ghost_sched%n_recv)

        ! Each rank must receive exactly what its neighbours send it
        call MPI_ALLTOALL(send_per_rank, 1, MPI_INTEGER, expect, 1, MPI_INTEGER, ice_comm, ierr)
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

        call MPI_ALLREDUCE(ghost_sched%n_send, ns, 1, MPI_INTEGER, MPI_SUM, ice_comm, ierr)
        call check_mpi_error(ierr)
        ns_total = ns_total + ns

        call build_conn_lists(p)
        call MPI_ALLREDUCE(ghost_sched%cf_n_send, ns, 1, MPI_INTEGER, MPI_SUM, ice_comm, ierr)
        call check_mpi_error(ierr)
        nf_total = nf_total + ns

        deallocate(send_per_rank, recv_per_rank, expect, send_pos, recv_pos)
      end block
      ghost_sched%built = .true.
    end do
    ! The sum over the families: what a set-up fill exchanges; a stage of
    ! family p exchanges p's share of it. Sent from a copy: no MPI_IN_PLACE
    ! (see mpi_allreduce_sum_r8_array).
    block
      use mpi
      integer :: ierr, ndup_rank
      ndup_rank = ndup
      call MPI_ALLREDUCE(ndup_rank, ndup, 1, MPI_INTEGER, MPI_SUM, ice_comm, ierr)
      call check_mpi_error(ierr)
    end block
    if (mpi_is_root) then
      write(*,'(A,I0,A)') '  MPI halo: ', ns_total, ' cells exchanged per ghost fill'
      if (ndup > 0) write(*,'(A,I0,A)') '  MPI halo: ', ndup, ' records that read a cell already sent, merged'
      if (nf_total > 0) write(*,'(A,I0,A)') '  MPI halo: ', nf_total, ' connection faces whose flux crosses ranks'
    end if
    ghost_sched => mg_sched(level, 1)
#endif

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


    !> Keep the first occurrence of every cell in each rank's run of list, in
    !> place; per_rank and n become the unique counts. The mask covers a block's
    !> interior and both ghost layers (a connection of a one-cell block reads its
    !> own ghost) and is cleared run by run, so the cost is linear in the list.
    subroutine unique_cells(list, per_rank, n)
      type(halo_cell_type), intent(inout) :: list(:)
      integer, intent(inout) :: per_rank(0:)
      integer, intent(inout) :: n
      type(mask_type), allocatable :: m(:)
      integer :: r, a, c, w, kept, b

      allocate(m(grid%nb))
      w = 0 ; a = 1
      do r = 0, size(per_rank) - 1
        kept = 0
        do c = a, a + per_rank(r) - 1
          b = list(c)%b
          if (.not. allocated(m(b)%seen)) then
            allocate(m(b)%seen(-1:grid%blk(b)%dim(1)+2, -1:grid%blk(b)%dim(2)+2, -1:grid%blk(b)%dim(3)+2))
            m(b)%seen = .false.
          end if
          if (any([list(c)%i, list(c)%j, list(c)%k] < -1) .or. &
              any([list(c)%i, list(c)%j, list(c)%k] > grid%blk(b)%dim(1:3) + 2)) &
            call mpi_abort_all('build_ghost_schedule: a halo cell outside its block and its ghost layers')
          if (m(b)%seen(list(c)%i, list(c)%j, list(c)%k)) cycle
          m(b)%seen(list(c)%i, list(c)%j, list(c)%k) = .true.
          w = w + 1 ; kept = kept + 1
          list(w) = list(c)
        end do
        do c = w - kept + 1, w
          m(list(c)%b)%seen(list(c)%i, list(c)%j, list(c)%k) = .false.
        end do
        a = a + per_rank(r)
        per_rank(r) = kept
      end do
      n = w
    end subroutine unique_cells


    !> The connection faces of family p whose owner record and partner record lie in blocks of
    !> different ranks: the owner's rank sends the face's flux to the partner's. Pass 1 counts per
    !> rank, pass 2 places the face numbers grouped by rank, table order inside a group; every rank
    !> scans the same table in the same order, so both ends of a message list the same faces.
    subroutine build_conn_lists(p)
      use mpi
      integer, intent(in) :: p
      integer, allocatable :: ns(:), nr(:), spos(:), rpos(:), expect(:)
      integer :: pass, n, q, ro, rp, r, nv, ierr

      allocate(ns(0:mpi_size_-1), nr(0:mpi_size_-1), spos(0:mpi_size_-1), rpos(0:mpi_size_-1), &
               expect(0:mpi_size_-1))
      do pass = 1, 2
        ns = 0 ; nr = 0
        do n = 1, size(grid%bc)
          if (grid%bc(n)%p /= p) cycle
          if (grid%bc(n)%type /= 101 .and. grid%bc(n)%type /= 201) cycle
          if (.not. grid%bc(n)%owner) cycle
          q  = grid%bc(n)%partner
          ro = block_owner(grid%bc(n)%b)
          rp = block_owner(grid%bc(q)%b)
          if (ro == rp) cycle
          if (ro == mpi_rank_) then
            ns(rp) = ns(rp) + 1
            if (pass == 2) ghost_sched%cf_send(spos(rp) + ns(rp)) = grid%bc(n)%cf
          else if (rp == mpi_rank_) then
            nr(ro) = nr(ro) + 1
            if (pass == 2) ghost_sched%cf_recv(rpos(ro) + nr(ro)) = grid%bc(n)%cf
          end if
        end do
        if (pass == 1) then
          ghost_sched%cf_n_send = sum(ns) ; ghost_sched%cf_n_recv = sum(nr)
          allocate(ghost_sched%cf_send(max(1, ghost_sched%cf_n_send)), ghost_sched%cf_recv(max(1, ghost_sched%cf_n_recv)))
          spos(0) = 0 ; rpos(0) = 0
          do r = 1, mpi_size_-1
            spos(r) = spos(r-1) + ns(r-1)
            rpos(r) = rpos(r-1) + nr(r-1)
          end do
        end if
      end do

      ! Each rank must receive exactly the faces its neighbours send it
      call MPI_ALLTOALL(ns, 1, MPI_INTEGER, expect, 1, MPI_INTEGER, ice_comm, ierr)
      call check_mpi_error(ierr)
      if (count(expect /= nr) > 0) call mpi_abort_all('build_ghost_schedule: connection-face send and receive lists disagree')

      nv = ncond(p)
      call face_groups(ns, nv, ghost_sched%cf_send_groups, ghost_sched%cf_n_send_ranks)
      call face_groups(nr, nv, ghost_sched%cf_recv_groups, ghost_sched%cf_n_recv_ranks)
      allocate(ghost_sched%cf_send_buf(max(1, nv*ghost_sched%cf_n_send)))
      allocate(ghost_sched%cf_recv_buf(max(1, nv*ghost_sched%cf_n_recv)))
      allocate(ghost_sched%cf_send_req(max(1, ghost_sched%cf_n_send_ranks)))
      allocate(ghost_sched%cf_recv_req(max(1, ghost_sched%cf_n_recv_ranks)))
      allocate(ghost_sched%cf_stats(MPI_STATUS_SIZE, max(1, ghost_sched%cf_n_send_ranks, ghost_sched%cf_n_recv_ranks)))
      do r = 1, ghost_sched%cf_n_recv_ranks
        associate (g => ghost_sched%cf_recv_groups(r))
          call MPI_RECV_INIT(ghost_sched%cf_recv_buf(g%off+1), g%len, MPI_DOUBLE_PRECISION, &
                             g%rank, conn_tag + p, ice_comm, ghost_sched%cf_recv_req(r), ierr)
          call check_mpi_error(ierr)
        end associate
      end do
      do r = 1, ghost_sched%cf_n_send_ranks
        associate (g => ghost_sched%cf_send_groups(r))
          call MPI_SEND_INIT(ghost_sched%cf_send_buf(g%off+1), g%len, MPI_DOUBLE_PRECISION, &
                             g%rank, conn_tag + p, ice_comm, ghost_sched%cf_send_req(r), ierr)
          call check_mpi_error(ierr)
        end associate
      end do
      ! The partner cells' weights take the same messages the other way, one value per face: the owner's
      ! rank receives on its send groups, the partner's rank sends on its receive groups
      allocate(ghost_sched%cb_send_buf(max(1, ghost_sched%cf_n_recv)), ghost_sched%cb_recv_buf(max(1, ghost_sched%cf_n_send)))
      allocate(ghost_sched%cb_send_req(max(1, ghost_sched%cf_n_recv_ranks)))
      allocate(ghost_sched%cb_recv_req(max(1, ghost_sched%cf_n_send_ranks)))
      do r = 1, ghost_sched%cf_n_send_ranks
        associate (g => ghost_sched%cf_send_groups(r))
          call MPI_RECV_INIT(ghost_sched%cb_recv_buf(g%first), g%count, MPI_DOUBLE_PRECISION, &
                             g%rank, cbeta_tag + p, ice_comm, ghost_sched%cb_recv_req(r), ierr)
          call check_mpi_error(ierr)
        end associate
      end do
      do r = 1, ghost_sched%cf_n_recv_ranks
        associate (g => ghost_sched%cf_recv_groups(r))
          call MPI_SEND_INIT(ghost_sched%cb_send_buf(g%first), g%count, MPI_DOUBLE_PRECISION, &
                             g%rank, cbeta_tag + p, ice_comm, ghost_sched%cb_send_req(r), ierr)
          call check_mpi_error(ierr)
        end associate
      end do
      ghost_sched%cf_persistent_init = .true.
      deallocate(ns, nr, spos, rpos, expect)
    end subroutine build_conn_lists


    !> One message per remote rank with faces to exchange: per_rank(r) faces of nv values each,
    !> in rank order, as build_conn_lists placed them.
    subroutine face_groups(per_rank, nv, groups, n_groups)
      integer, intent(in) :: per_rank(0:), nv
      type(rank_group_type), allocatable, intent(out) :: groups(:)
      integer, intent(out) :: n_groups
      integer :: r, g, first

      n_groups = count(per_rank > 0)
      allocate(groups(n_groups))
      g = 0 ; first = 1
      do r = 0, size(per_rank) - 1
        if (per_rank(r) == 0) cycle
        g = g + 1
        groups(g) = rank_group_type(r, first, per_rank(r), (first - 1) * nv, per_rank(r) * nv)
        first = first + per_rank(r)
      end do
    end subroutine face_groups
#endif

  end subroutine build_ghost_schedule


  !> Point the halo routines at one level's schedules. Call whenever the solver moves
  !> between multigrid levels; build_ghost_schedule has to have run for that level.
  subroutine select_ghost_level(level)
    implicit none
    integer, intent(in) :: level

    if (.not. allocated(mg_sched)) return
    if (level < 1 .or. level > size(mg_sched, 1)) &
      call mpi_abort_all('select_ghost_level: grid level outside 1..MG-levels')
    level_cur = level
    ghost_sched => mg_sched(level, 1)

  end subroutine select_ghost_level


  !> Free persistent MPI requests, on every level and family that built any.
  subroutine cleanup_ghost_schedule()
    implicit none
#ifdef USE_MPI
    integer :: m, p

    if (.not. allocated(mg_sched)) return
    do m = 1, size(mg_sched, 1)
      do p = 1, size(mg_sched, 2)
        ghost_sched => mg_sched(m, p)
        call cleanup_persistent_requests()
      end do
    end do
    ghost_sched => mg_sched(1, 1)
#endif
  end subroutine cleanup_ghost_schedule


  !> Bring the remote interior cells read by the ghost fill up to date: family p's
  !> when p is given (the stage of one family), every family's otherwise (set-up).
  !> Call on every rank before compute_ghost, on the level select_ghost_level chose.
  subroutine exchange_ghost_prim(grid, p)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer, intent(in), optional        :: p
    integer :: q

    if (mpi_size_ <= 1) return
    if (size(grid%bc) == 0) return
    if (present(p)) then
      call exchange_family(grid, p)
    else
      do q = 1, ngroups
        call exchange_family(grid, q)
      end do
    end if
    ghost_sched => mg_sched(level_cur, 1)
  end subroutine exchange_ghost_prim


  !> The exchange of one family: start the receives, pack, start the sends, wait,
  !> unpack. The schedule of (level_cur, p) is used.
  subroutine exchange_family(grid, p)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer, intent(in)                  :: p

    ghost_sched => mg_sched(level_cur, p)
#ifdef USE_MPI
    block
      use mpi
      integer :: n, nv, off, ierr

      if (ghost_sched%n_recv_ranks > 0) then
        call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, ierr)
        call check_mpi_error(ierr)
      end if

      ! Every cell owns a slice of the buffer, so the threads write disjoint ranges;
      ! the MPI calls stay on the thread that calls this routine (FUNNELED).
      call timer_region_begin(TR_PACK)
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ghost_sched, grid, ncond) PRIVATE(n, nv, off) SCHEDULE(STATIC)
      do n = 1, ghost_sched%n_send
        associate (c => ghost_sched%send_list(n))
          nv = ncond(c%p) ; off = c%off
          ghost_sched%send_buf(off+1:off+nv) = grid%blk(c%b)%cond_phase(c%p)%prim(1:nv, c%i, c%j, c%k)
        end associate
      end do
      !$OMP END PARALLEL DO
      call timer_region_end(TR_PACK)

      if (ghost_sched%n_send_ranks > 0) then
        call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%send_req, ierr)
        call check_mpi_error(ierr)
      end if

      call timer_region_begin(TR_WAIT)
      if (ghost_sched%n_recv_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, ghost_sched%stats, ierr)
        call check_mpi_error(ierr)
      end if
      call timer_region_end(TR_WAIT)

      ! A received cell is listed once (unique_cells) and belongs to one source
      ! rank, so the threads write disjoint cells.
      call timer_region_begin(TR_UNPACK)
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ghost_sched, grid, ncond) PRIVATE(n, nv, off) SCHEDULE(STATIC)
      do n = 1, ghost_sched%n_recv
        associate (c => ghost_sched%recv_list(n))
          nv = ncond(c%p) ; off = c%off
          grid%blk(c%b)%cond_phase(c%p)%prim(1:nv, c%i, c%j, c%k) = ghost_sched%recv_buf(off+1:off+nv)
        end associate
      end do
      !$OMP END PARALLEL DO
      call timer_region_end(TR_UNPACK)

      call timer_region_begin(TR_WAIT)
      if (ghost_sched%n_send_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%send_req, ghost_sched%stats, ierr)
        call check_mpi_error(ierr)
      end if
      call timer_region_end(TR_WAIT)
    end block
#endif
  end subroutine exchange_family


  !> The connection-face fluxes of family p whose owner record lies on another rank than its
  !> partner: the owner's rank sends the flux its owner pass computed (compute_bound), the
  !> partner's rank writes it into its copy of the face, where the partner record reads it.
  !> Call on every rank, outside a parallel region, between that pass and the boundary cells,
  !> on the level select_ghost_level chose.
  subroutine exchange_conn_flux(grid, p)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer, intent(in)                  :: p

    if (mpi_size_ <= 1) return
    ghost_sched => mg_sched(level_cur, p)
#ifdef USE_MPI
    block
      use mpi
      integer :: n, nv, ierr

      nv = ncond(p)
      if (ghost_sched%cf_n_recv_ranks > 0) then
        call MPI_STARTALL(ghost_sched%cf_n_recv_ranks, ghost_sched%cf_recv_req, ierr)
        call check_mpi_error(ierr)
      end if
      do n = 1, ghost_sched%cf_n_send
        ghost_sched%cf_send_buf((n-1)*nv+1:n*nv) = grid%conn(p)%F(1:nv, ghost_sched%cf_send(n))
      end do
      if (ghost_sched%cf_n_send_ranks > 0) then
        call MPI_STARTALL(ghost_sched%cf_n_send_ranks, ghost_sched%cf_send_req, ierr)
        call check_mpi_error(ierr)
      end if
      if (ghost_sched%cf_n_recv_ranks > 0) then
        call MPI_WAITALL(ghost_sched%cf_n_recv_ranks, ghost_sched%cf_recv_req, ghost_sched%cf_stats, ierr)
        call check_mpi_error(ierr)
      end if
      do n = 1, ghost_sched%cf_n_recv
        grid%conn(p)%F(1:nv, ghost_sched%cf_recv(n)) = ghost_sched%cf_recv_buf((n-1)*nv+1:n*nv)
      end do
      if (ghost_sched%cf_n_send_ranks > 0) then
        call MPI_WAITALL(ghost_sched%cf_n_send_ranks, ghost_sched%cf_send_req, ghost_sched%cf_stats, ierr)
        call check_mpi_error(ierr)
      end if
    end block
#endif
    ghost_sched => mg_sched(level_cur, 1)
  end subroutine exchange_conn_flux


  !> The detector weights of family p's partner cells that the owner pass needs from another rank:
  !> the rank of a face's partner record sends the weight of its cell, computed by compute_beta in
  !> this stage, to the rank of the owner record, which keeps it in conn(p)%bp. Call on every rank,
  !> outside a parallel region, after compute_beta and before the owner pass (compute_bound).
  subroutine exchange_conn_beta(grid, p)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer, intent(in)                  :: p

    if (mpi_size_ <= 1) return
    ghost_sched => mg_sched(level_cur, p)
#ifdef USE_MPI
    block
      use mpi
      integer :: n, q, ierr

      if (ghost_sched%cf_n_send_ranks > 0) then
        call MPI_STARTALL(ghost_sched%cf_n_send_ranks, ghost_sched%cb_recv_req, ierr)
        call check_mpi_error(ierr)
      end if
      do n = 1, ghost_sched%cf_n_recv
        q = grid%conn(p)%prt(ghost_sched%cf_recv(n))
        ghost_sched%cb_send_buf(n) = grid%blk(grid%bc(q)%b)%cond_phase(p)%beta(grid%bc(q)%i, grid%bc(q)%j, grid%bc(q)%k)
      end do
      if (ghost_sched%cf_n_recv_ranks > 0) then
        call MPI_STARTALL(ghost_sched%cf_n_recv_ranks, ghost_sched%cb_send_req, ierr)
        call check_mpi_error(ierr)
      end if
      if (ghost_sched%cf_n_send_ranks > 0) then
        call MPI_WAITALL(ghost_sched%cf_n_send_ranks, ghost_sched%cb_recv_req, ghost_sched%cf_stats, ierr)
        call check_mpi_error(ierr)
      end if
      do n = 1, ghost_sched%cf_n_send
        grid%conn(p)%bp(ghost_sched%cf_send(n)) = ghost_sched%cb_recv_buf(n)
      end do
      if (ghost_sched%cf_n_recv_ranks > 0) then
        call MPI_WAITALL(ghost_sched%cf_n_recv_ranks, ghost_sched%cb_send_req, ghost_sched%cf_stats, ierr)
        call check_mpi_error(ierr)
      end if
    end block
#endif
    ghost_sched => mg_sched(level_cur, 1)
  end subroutine exchange_conn_beta


  !> How many of this rank's boundary records carry each type. The face mix is
  !> what compute_bound's cost follows -- a symmetry face is a reflection, a
  !> connection a copy -- and a decomposition that balances cells does not
  !> balance it; the decomposer's cost model is calibrated on these counts.
  !> Printed once at set-up, when the timers are on, as max and sum over ranks.
  subroutine report_bc_mix(grid)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    ! Type 0 is a face without a condition (the k faces of a 2-D case).
    integer, parameter :: codes(11) = [0, 101, 102, 200, 201, 300, 301, 400, 401, 402, 403]
    real(R8) :: cmax(12), csum(12)
    integer  :: i, n, c, t
    character(len=8) :: label(12)

    cmax = 0.0_R8
    do n = 1, grid%n_local_bc
      i = grid%local_bc_idx(n)
      t = 12
      do c = 1, size(codes)
        if (grid%bc(i)%type == codes(c)) then
          t = c
          exit
        end if
      end do
      cmax(t) = cmax(t) + 1.0_R8
    end do
    csum = cmax
    call mpi_allreduce_max_r8_array(cmax, 12)
    call mpi_allreduce_sum_r8_array(csum, 12)
    if (.not. mpi_is_root) return
    do c = 1, size(codes)
      write(label(c), '(I0)') codes(c)
    end do
    label(12) = 'other'
    write(*,'(A,12(1X,A,1X,I0,A,I0))') ' ICE BCmix  |', &
      (trim(label(c)), nint(cmax(c)), '/', nint(csum(c)), c = 1, 12)
  end subroutine report_bc_mix


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
                           block_owner(b), io_tag, ice_comm, reqs(nreq), ierr)
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
                           0, io_tag, ice_comm, reqs(nreq), ierr)
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
                           block_owner(b), io_tag, ice_comm, reqs(nreq), ierr)
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
                           0, io_tag, ice_comm, reqs(nreq), ierr)
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
                           block_owner(b), io_tag, ice_comm, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: post receives for local blocks
        do b = 1, grid%nb
          if (is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, io_tag, ice_comm, reqs(nreq), ierr)
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
      call MPI_BARRIER(ice_comm, ierr)
      call check_mpi_error(ierr)
    end block
#endif
  end subroutine mpi_io_barrier


#ifdef USE_MPI
  !> Initialize one persistent receive and one persistent send per neighbour rank,
  !> and the statuses their waits fill. MPI_STATUSES_IGNORE, like MPI_IN_PLACE,
  !> lives in a COMMON block of the mpi module, which -fno-underscoring renames
  !> (see mpi_allreduce_sum_r8_array): the library would write the statuses
  !> through an ordinary address.
  subroutine init_persistent_requests()
    use mpi
    implicit none
    integer :: r, ierr

    allocate(ghost_sched%recv_req(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%send_req(max(ghost_sched%n_send_ranks, 1)))
    allocate(ghost_sched%stats(MPI_STATUS_SIZE, max(ghost_sched%n_recv_ranks, ghost_sched%n_send_ranks, 1)))

    do r = 1, ghost_sched%n_recv_ranks
      associate (g => ghost_sched%recv_groups(r))
        call MPI_RECV_INIT(ghost_sched%recv_buf(g%off+1), g%len, MPI_DOUBLE_PRECISION, &
                           g%rank, halo_tag + ghost_sched%group, ice_comm, ghost_sched%recv_req(r), ierr)
        call check_mpi_error(ierr)
      end associate
    end do

    do r = 1, ghost_sched%n_send_ranks
      associate (g => ghost_sched%send_groups(r))
        call MPI_SEND_INIT(ghost_sched%send_buf(g%off+1), g%len, MPI_DOUBLE_PRECISION, &
                           g%rank, halo_tag + ghost_sched%group, ice_comm, ghost_sched%send_req(r), ierr)
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

    if (ghost_sched%cf_persistent_init) then
      do r = 1, ghost_sched%cf_n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%cf_send_req(r), ierr)
        call MPI_REQUEST_FREE(ghost_sched%cb_recv_req(r), ierr)
      end do
      do r = 1, ghost_sched%cf_n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%cf_recv_req(r), ierr)
        call MPI_REQUEST_FREE(ghost_sched%cb_send_req(r), ierr)
      end do
      ghost_sched%cf_persistent_init = .false.
    end if
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
