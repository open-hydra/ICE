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
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k

    !$OMP PARALLEL
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$OMP DO COLLAPSE (3)
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
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k

    !$OMP PARALLEL
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$omp do collapse (3)
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


  subroutine compute_residual_(source, dt, volume, residual)
    implicit none
    real(kind=R8), dimension(:), intent(in)    :: source
    real(kind=R8),               intent(in)    :: dt
    real(kind=R8),               intent(in)    :: volume
    real(kind=R8), dimension(:), intent(inout) :: residual

    residual = (residual/volume + source) * dt

  end subroutine compute_residual_


end module ICE_Lib_Residual
