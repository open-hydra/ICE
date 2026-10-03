!>@brief Thread groups: the blocks of a rank owned by groups of its threads.
!> A rank whose team spans several sockets (one rank of 80 threads on four sockets)
!> would otherwise spread every block over every thread: each block's pages are then
!> shared by all the sockets, and on many small blocks the kernel's NUMA balancer
!> keeps migrating them. With ICE_THREAD_GROUPS = G the team is cut into G groups of
!> contiguous thread ids -- with OMP_PROC_BIND=close a group is one socket -- and each
!> local block is owned by one group (longest-processing-time on cell counts, as the
!> blocks over the ranks): its fields are first touched, and its cells swept, only by
!> the threads of that group. Which thread computes a cell is not arithmetic, so no
!> result depends on G. Unset or `auto` (the default): with the threads pinned
!> (OMP_PROC_BIND), one group per socket the team spans, read from where each thread
!> runs; unpinned, one group. G = 1 is the whole team on every block.
module ICE_Mod_ThreadGroups
  use ICE_Mod_MPI, only: lpt_assign, local_block_ids, n_local_blocks, mpi_is_root
  use iso_c_binding, only: c_int
  !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_proc_bind, omp_proc_bind_false
  implicit none
  private

  integer, public :: n_tgroups = 1                    !< Thread groups in use (1: every block over the whole team)
  integer, public :: n_tthreads = 1                   !< Threads of the team the groups cut
  integer, allocatable, public :: tgroup_of_block(:)  !< Group (0-based) owning block b; -1 for a block of another rank
  integer, allocatable, public :: tg_first(:)         !< Blocks of group g: tg_blocks(tg_first(g):tg_first(g+1)-1)
  integer, allocatable, public :: tg_blocks(:)        !< ... in mesh order
  integer, allocatable, public :: tg_thr_lo(:)        !< First thread id of group g
  integer, allocatable, public :: tg_thr_hi(:)        !< Last thread id of group g
  integer, public :: tg_rounds = 0                    !< Most blocks owned by one group

  interface
    function sched_getcpu() bind(C, name='sched_getcpu')
      import :: c_int
      integer(c_int) :: sched_getcpu
    end function sched_getcpu
  end interface

  public :: setup_thread_groups, thread_group, split_range, block_group, cell_slice

contains

  !> Read ICE_THREAD_GROUPS and assign the rank's blocks to the groups. Called once,
  !> after partition_blocks; nb and blk_ncells as partition_blocks takes them. Fewer
  !> local blocks than groups, or fewer threads than groups, means one group.
  subroutine setup_thread_groups(nb, blk_ncells)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    character(len=32) :: env
    integer :: status, ios, asked, g, i, nloc
    integer, allocatable :: owner(:), load(:), w(:)

    n_tthreads = 1
    !$ n_tthreads = omp_get_max_threads()

    asked = 0                                   ! 0: auto
    call get_environment_variable('ICE_THREAD_GROUPS', env, status=status)
    if (status == 0 .and. trim(adjustl(env)) /= 'auto') then
      read(env, *, iostat=ios) asked
      if (ios /= 0 .or. asked < 1) asked = 1
    end if

    nloc = n_local_blocks
    if (allocated(tgroup_of_block)) deallocate(tgroup_of_block, tg_first, tg_blocks, tg_thr_lo, tg_thr_hi)
    if (asked == 0) then
      call socket_runs(n_tgroups)               ! sets tg_thr_lo/hi when it finds more than one socket
    else
      n_tgroups = asked
    end if
    if (n_tgroups > n_tthreads .or. n_tgroups > nloc) n_tgroups = 1
    if (asked /= 0 .or. n_tgroups == 1) then
      if (allocated(tg_thr_lo)) deallocate(tg_thr_lo, tg_thr_hi)
      allocate(tg_thr_lo(0:n_tgroups-1), tg_thr_hi(0:n_tgroups-1))
      do g = 0, n_tgroups - 1
        tg_thr_lo(g) = (g * n_tthreads) / n_tgroups
        tg_thr_hi(g) = ((g + 1) * n_tthreads) / n_tgroups - 1
      end do
    end if

    allocate(tgroup_of_block(nb), tg_first(0:n_tgroups), tg_blocks(max(nloc, 1)))
    tgroup_of_block = -1

    allocate(owner(max(nloc, 1)), load(0:n_tgroups-1), w(max(nloc, 1)))
    do i = 1, nloc
      w(i) = blk_ncells(local_block_ids(i))
    end do
    call lpt_assign(nloc, w(1:nloc), n_tgroups, owner(1:nloc), load)
    do i = 1, nloc
      tgroup_of_block(local_block_ids(i)) = owner(i)
    end do

    ! Blocks of each group in mesh order (local_block_ids is in mesh order)
    tg_first(0) = 1
    do g = 0, n_tgroups - 1
      tg_first(g+1) = tg_first(g)
      do i = 1, nloc
        if (owner(i) /= g) cycle
        tg_blocks(tg_first(g+1)) = local_block_ids(i)
        tg_first(g+1) = tg_first(g+1) + 1
      end do
    end do
    tg_rounds = 0
    do g = 0, n_tgroups - 1
      tg_rounds = max(tg_rounds, tg_first(g+1) - tg_first(g))
    end do

    if (mpi_is_root .and. (status == 0 .or. n_tgroups > 1)) then
      if (asked == 0) then
        write(*,'(A,I0,A,I0,A)') '  OpenMP: thread groups ', n_tgroups, ' (auto: sockets spanned) over ', &
          n_tthreads, ' threads'
      else
        write(*,'(A,I0,A,I0,A,I0,A)') '  OpenMP: thread groups ', n_tgroups, ' (asked ', asked, ') over ', &
          n_tthreads, ' threads'
      end if
      if (n_tgroups > 1) then
        do g = 0, n_tgroups - 1
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '    group ', g, ': threads ', tg_thr_lo(g), '-', tg_thr_hi(g), &
            ', blocks ', tg_first(g+1) - tg_first(g), ', cells ', load(g)
        end do
      end if
    end if
    deallocate(owner, load, w)
  end subroutine setup_thread_groups


  !> Auto mode: the sockets the team spans, from where each thread runs. Only with the
  !> threads pinned (an unpinned thread's CPU at set-up says nothing about later), and
  !> only when the threads of each socket are consecutive ids (OMP_PROC_BIND=close or
  !> spread over a socket-major core list); otherwise one group. On success the
  !> groups are the runs of consecutive threads on one socket: tg_thr_lo/hi are set.
  subroutine socket_runs(ng)
    integer, intent(out) :: ng
    integer, allocatable :: cpu(:), pkg(:)
    integer :: t, u, ios, nrun
    character(len=96) :: path
    logical :: pinned

    ng = 1
    pinned = .false.
    !$ pinned = omp_get_proc_bind() /= omp_proc_bind_false
    if (.not. pinned .or. n_tthreads < 2) return

    allocate(cpu(0:n_tthreads-1), pkg(0:n_tthreads-1))
    cpu = -1
    !$OMP PARALLEL
    !$ cpu(omp_get_thread_num()) = int(sched_getcpu())
    !$OMP END PARALLEL
    do t = 0, n_tthreads - 1
      if (cpu(t) < 0) return
      write(path,'(A,I0,A)') '/sys/devices/system/cpu/cpu', cpu(t), '/topology/physical_package_id'
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read(u, *, iostat=ios) pkg(t)
      close(u)
      if (ios /= 0) return
    end do

    ! Runs of one socket; a socket that comes back later means the threads are not
    ! grouped by socket: one group
    nrun = 1
    do t = 1, n_tthreads - 1
      if (pkg(t) /= pkg(t-1)) then
        if (any(pkg(0:t-1) == pkg(t))) return
        nrun = nrun + 1
      end if
    end do
    if (nrun == 1) return

    if (allocated(tg_thr_lo)) deallocate(tg_thr_lo, tg_thr_hi)
    allocate(tg_thr_lo(0:nrun-1), tg_thr_hi(0:nrun-1))
    nrun = 0
    tg_thr_lo(0) = 0
    do t = 1, n_tthreads - 1
      if (pkg(t) /= pkg(t-1)) then
        tg_thr_hi(nrun) = t - 1
        nrun = nrun + 1
        tg_thr_lo(nrun) = t
      end if
    end do
    tg_thr_hi(nrun) = n_tthreads - 1
    ng = nrun + 1
  end subroutine socket_runs


  !> The group owning block b, or -1 when the whole team works on it: one group,
  !> no groups set up yet, or a block of another rank.
  pure integer function block_group(b)
    integer, intent(in) :: b
    block_group = -1
    if (n_tgroups <= 1 .or. .not. allocated(tgroup_of_block)) return
    if (b < 1 .or. b > size(tgroup_of_block)) return
    block_group = tgroup_of_block(b)
  end function block_group


  !> The group of thread tid (0-based), and the thread's index m inside it and the
  !> group's size s.
  pure subroutine thread_group(tid, g, m, s)
    integer, intent(in)  :: tid
    integer, intent(out) :: g, m, s
    g = n_tgroups - 1
    do while (g > 0)
      if (tg_thr_lo(g) <= tid) exit
      g = g - 1
    end do
    m = tid - tg_thr_lo(g)
    s = tg_thr_hi(g) - tg_thr_lo(g) + 1
  end subroutine thread_group


  !> The m-th (0-based) of s contiguous, near-equal pieces of a..b: i1..i2 (empty
  !> when i1 > i2). Static, so the same thread gets the same planes in every phase.
  pure subroutine split_range(a, b, m, s, i1, i2)
    integer, intent(in)  :: a, b, m, s
    integer, intent(out) :: i1, i2
    integer :: n
    n  = b - a + 1
    i1 = a + (m * n) / s
    i2 = a + ((m + 1) * n) / s - 1
  end subroutine split_range


  !> The cells of a block (dims) that thread m of a group of s sweeps: the m-th of s
  !> contiguous, near-equal pieces of its cells in storage order (i fastest) -- the
  !> split a static collapse(3) schedule makes. From (i1,j1,k1) to (i2,j2,k2) in that
  !> order; empty when k1 > k2. Walk it as
  !>   do k = k1, k2;  j from (j1 if k == k1 else 1) to (j2 if k == k2 else nj)
  !>     i from (i1 if k == k1 and j == j1 else 1) to (i2 if k == k2 and j == j2 else ni)
  pure subroutine cell_slice(dims, m, s, i1, j1, k1, i2, j2, k2)
    use iso_fortran_env, only: int64
    integer, intent(in)  :: dims(3), m, s
    integer, intent(out) :: i1, j1, k1, i2, j2, k2
    integer(int64) :: n, nij, c1, c2
    nij = int(dims(1), int64) * dims(2)
    n   = nij * dims(3)
    c1  = (int(m, int64) * n) / s
    c2  = (int(m + 1, int64) * n) / s - 1
    if (c1 > c2) then
      i1 = 1 ; j1 = 1 ; k1 = 1 ; i2 = 0 ; j2 = 0 ; k2 = 0
      return
    end if
    k1 = int(c1 / nij) + 1 ; j1 = int(mod(c1, nij) / dims(1)) + 1 ; i1 = int(mod(c1, int(dims(1), int64))) + 1
    k2 = int(c2 / nij) + 1 ; j2 = int(mod(c2, nij) / dims(1)) + 1 ; i2 = int(mod(c2, int(dims(1), int64))) + 1
  end subroutine cell_slice

end module ICE_Mod_ThreadGroups
