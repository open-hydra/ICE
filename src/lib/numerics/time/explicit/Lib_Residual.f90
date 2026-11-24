module ICE_Lib_Residual
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Global_m
  use ICE_Advanced_Types_m

  implicit none
  private
  public :: zero_residual, compute_residual, add_births

contains

  subroutine zero_residual(grid, p)
    use ICE_Mod_ThreadGroups, only: n_tgroups
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k

    if (n_tgroups > 1) then
      call zero_residual_grouped(grid, p)
      return
    end if

    !$OMP PARALLEL
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$OMP DO COLLAPSE (3) schedule(runtime)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) = 0._R8

      enddo ; enddo ; enddo
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL

  end subroutine zero_residual


  !> The residual of the stage. The slots frozen at set-up (setup_freeze) are not integrated: their
  !> residual is zero, and residual smoothing keeps it zero.
  subroutine compute_residual(grid, p)
    use ICE_Mod_ThreadGroups, only: n_tgroups
    use ICE_Mod_Freeze,       only: nfrozen, frozen
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, nf

    if (n_tgroups > 1) then
      call compute_residual_grouped(grid, p)
      return
    end if

    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)

    !$OMP PARALLEL
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$omp do collapse (3) schedule(runtime)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        call compute_residual_(grid%blk(b)%cond_phase(p)%source(1:ncond(p),i,j,k),  &
                               grid%blk(b)%cond_phase(p)%dt(i,j,k),                 &
                               grid%blk(b)%vol(i,j,k),                              &
                               grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) )
        if (nf > 0) grid%blk(b)%cond_phase(p)%residual(frozen(1:nf,p),i,j,k) = 0._R8

      enddo ; enddo ; enddo
      !$omp end do nowait
    enddo
    !$OMP END PARALLEL

  end subroutine compute_residual


  !> zero_residual with thread groups: each thread sweeps its slice of its group's blocks.
  subroutine zero_residual_grouped(grid, p)
    use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, thread_group, cell_slice
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, n, tid, g, m, s, i1, j1, k1, i2, j2, k2

    !$OMP PARALLEL PRIVATE(b, i, j, k, n, tid, g, m, s, i1, j1, k1, i2, j2, k2)
    tid = 0
    !$ tid = omp_get_thread_num()
    call thread_group(tid, g, m, s)
    do n = tg_first(g), tg_first(g+1) - 1
      b = tg_blocks(n)
      call cell_slice(grid%blk(b)%dim, m, s, i1, j1, k1, i2, j2, k2)
      do k = k1, k2
      do j = merge(j1, 1, k == k1), merge(j2, grid%blk(b)%dim(2), k == k2)
      do i = merge(i1, 1, k == k1 .and. j == j1), merge(i2, grid%blk(b)%dim(1), k == k2 .and. j == j2)
          grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) = 0._R8
      enddo ; enddo ; enddo
    enddo
    !$OMP END PARALLEL

  end subroutine zero_residual_grouped


  !> compute_residual with thread groups: each thread sweeps its slice of its group's blocks.
  subroutine compute_residual_grouped(grid, p)
    use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, thread_group, cell_slice
    use ICE_Mod_Freeze,       only: nfrozen, frozen
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, n, tid, g, m, s, i1, j1, k1, i2, j2, k2, nf

    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)

    !$OMP PARALLEL PRIVATE(b, i, j, k, n, tid, g, m, s, i1, j1, k1, i2, j2, k2)
    tid = 0
    !$ tid = omp_get_thread_num()
    call thread_group(tid, g, m, s)
    do n = tg_first(g), tg_first(g+1) - 1
      b = tg_blocks(n)
      call cell_slice(grid%blk(b)%dim, m, s, i1, j1, k1, i2, j2, k2)
      do k = k1, k2
      do j = merge(j1, 1, k == k1), merge(j2, grid%blk(b)%dim(2), k == k2)
      do i = merge(i1, 1, k == k1 .and. j == j1), merge(i2, grid%blk(b)%dim(1), k == k2 .and. j == j2)
          call compute_residual_(grid%blk(b)%cond_phase(p)%source(1:ncond(p),i,j,k),  &
                                 grid%blk(b)%cond_phase(p)%dt(i,j,k),                 &
                                 grid%blk(b)%vol(i,j,k),                              &
                                 grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) )
          if (nf > 0) grid%blk(b)%cond_phase(p)%residual(frozen(1:nf,p),i,j,k) = 0._R8
      enddo ; enddo ; enddo
    enddo
    !$OMP END PARALLEL

  end subroutine compute_residual_grouped


  subroutine compute_residual_(source, dt, volume, residual)
    implicit none
    real(kind=R8), dimension(:), intent(in)    :: source
    real(kind=R8),               intent(in)    :: dt
    real(kind=R8),               intent(in)    :: volume
    real(kind=R8), dimension(:), intent(inout) :: residual

    residual = (residual/volume + source) * dt

  end subroutine compute_residual_


  !> The births in force on family p, added to the residual of the stage over the family's own dt, as its source
  !> is: those a coupler set (births_ext, add_birth_rate) and, in a steady run, those of the mass ICE's own burning
  !> families generated (births_gen, gen_birth). Each is a birth packet per unit volume and time, [mass, momentum (3),
  !> energy, number], which goes to slots 1:4, nbase-1 and nbase. The frozen slots stay zero. Neither is the family's
  !> source, so neither reaches the gas.
  subroutine add_births(grid, p)
    use ICE_Mod_Freeze, only: nfrozen, frozen
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, nf

    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      if (allocated(grid%blk(b)%cond_phase(p)%births_ext)) &
        call add_packets(grid%blk(b), grid%blk(b)%cond_phase(p)%births_ext)
      if (allocated(grid%blk(b)%cond_phase(p)%births_gen)) &
        call add_packets(grid%blk(b), grid%blk(b)%cond_phase(p)%births_gen)
    enddo

  contains

    !> Adds the packets pk (6, cells) of block blk times the cells' dt to family p's residual
    subroutine add_packets(blk, pk)
      type(ICE_block_type), intent(inout) :: blk
      real(kind=R8),        intent(in)    :: pk(:,:,:,:)
      integer(kind=I4) :: i, j, k, ien, inum
      real(kind=R8)    :: dt

      ! the slots of the energy and of the number
      ien  = nbase(p) - 1
      inum = nbase(p)
      !$OMP PARALLEL DO COLLAPSE (3) PRIVATE(dt) schedule(runtime)
      do k = 1, blk%dim(3)
      do j = 1, blk%dim(2)
      do i = 1, blk%dim(1)
        dt = blk%cond_phase(p)%dt(i,j,k)
        blk%cond_phase(p)%residual(1:4,i,j,k)  = blk%cond_phase(p)%residual(1:4,i,j,k)  + pk(1:4,i,j,k)*dt
        blk%cond_phase(p)%residual(ien,i,j,k)  = blk%cond_phase(p)%residual(ien,i,j,k)  + pk(5,i,j,k)*dt
        blk%cond_phase(p)%residual(inum,i,j,k) = blk%cond_phase(p)%residual(inum,i,j,k) + pk(6,i,j,k)*dt
        if (nf > 0) blk%cond_phase(p)%residual(frozen(1:nf,p),i,j,k) = 0._R8
      enddo ; enddo ; enddo
      !$OMP END PARALLEL DO
    end subroutine add_packets

  end subroutine add_births


end module ICE_Lib_Residual
