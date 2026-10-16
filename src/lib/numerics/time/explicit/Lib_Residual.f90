module ICE_Lib_Residual
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Global_m
  use ICE_Advanced_Types_m

  implicit none
  private
  public :: zero_residual, compute_residual

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


  subroutine compute_residual(grid, p)
    use ICE_Mod_ThreadGroups, only: n_tgroups
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k

    if (n_tgroups > 1) then
      call compute_residual_grouped(grid, p)
      return
    end if

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
          call compute_residual_(grid%blk(b)%cond_phase(p)%source(1:ncond(p),i,j,k),  &
                                 grid%blk(b)%cond_phase(p)%dt(i,j,k),                 &
                                 grid%blk(b)%vol(i,j,k),                              &
                                 grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) )
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


end module ICE_Lib_Residual
