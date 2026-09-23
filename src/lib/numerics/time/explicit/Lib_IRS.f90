module ICE_Lib_IRS
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Global_m
  use ICE_Advanced_Types_m

  implicit none
  private
  public :: residual_smoothing

  integer(kind=I4), parameter, private :: njb = 2

contains

  subroutine residual_smoothing(grid, p)
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: d, b

    do d = 1, ndir
      do b = 1, grid%nb
        if (.not. is_local_block(b)) cycle
        call residual_smoothing_(grid%blk(b)%cond_phase(p)%residual, &
                                 grid%blk(b)%dim(1),                 &
                                 grid%blk(b)%dim(2),                 &
                                 grid%blk(b)%dim(3),                 &
                                 d                                   )
      end do
    enddo
  end subroutine residual_smoothing


  subroutine residual_smoothing_(residual, Nx, Ny, Nz, d)
    use ICE_Config_Types_m, only: obj_irs
    implicit none
    real(kind=R8), dimension(:,:,:,:), intent(inout) :: residual
    integer(kind=I4), intent(in) :: Nx, Ny, Nz, d
    real(kind=R8), allocatable   :: residual_star(:,:,:,:), residual_new(:,:,:,:)
    integer(kind=I4)             :: i, j, k, sjb
    integer(kind=I4)             :: i1, j1, k1, i2, j2, k2
    real(kind=R8)                :: beta_irs

    beta_irs = obj_irs%beta

    allocate(residual_star(size(residual,1), 0:Nx+1, 0:Ny+1, 0:Nz+1))
    allocate(residual_new (size(residual,1),   1:Nx,   1:Ny,   1:Nz))

    ! Copy interior, then fill ghost cells with zero-gradient (Neumann) extrapolation
    ! Ghost cells are fixed for all Jacobi iterations (matches MOSE's Ghost_Residual_Extrapolation)
    residual_star = 0._R8
    residual_star(:,1:Nx,1:Ny,1:Nz) = residual
    select case (d)
      case (1)
        residual_star(:,0,   1:Ny,1:Nz) = residual_star(:,1, 1:Ny,1:Nz)
        residual_star(:,Nx+1,1:Ny,1:Nz) = residual_star(:,Nx,1:Ny,1:Nz)
      case (2)
        residual_star(:,1:Nx,0,   1:Nz) = residual_star(:,1:Nx,1, 1:Nz)
        residual_star(:,1:Nx,Ny+1,1:Nz) = residual_star(:,1:Nx,Ny,1:Nz)
      case (3)
        residual_star(:,1:Nx,1:Ny,0   ) = residual_star(:,1:Nx,1:Ny,1 )
        residual_star(:,1:Nx,1:Ny,Nz+1) = residual_star(:,1:Nx,1:Ny,Nz)
    end select

    do sjb = 1, njb

      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i1,j1,k1,i2,j2,k2)
      do k = 1, Nz
      do j = 1, Ny
      do i = 1, Nx
        call jacobi_stencil(i, j, k, d, i1, j1, k1, i2, j2, k2)
        residual_new(:,i,j,k) = (residual(:,i,j,k) + beta_irs*(residual_star(:,i1,j1,k1) + &
                                 residual_star(:,i2,j2,k2))) / (1._R8 + 2._R8*beta_irs)
      end do ; end do ; end do
      !$OMP END PARALLEL DO

      residual_star(:,1:Nx,1:Ny,1:Nz) = residual_new
      ! Ghost cells are NOT updated between Jacobi iterations (fixed Neumann from initial Res)

    end do

    residual = residual_new

    deallocate(residual_star, residual_new)

  end subroutine residual_smoothing_


  subroutine jacobi_stencil(i, j, k, d, i1, j1, k1, i2, j2, k2)
    implicit none
    integer(kind=I4), intent(in)  :: i, j, k, d
    integer(kind=I4), intent(out) :: i1, j1, k1, i2, j2, k2

    select case (d)
      case (1)
        i1 = i-1 ; j1 = j   ; k1 = k
        i2 = i+1 ; j2 = j   ; k2 = k
      case (2)
        i1 = i   ; j1 = j-1 ; k1 = k
        i2 = i   ; j2 = j+1 ; k2 = k
      case (3)
        i1 = i   ; j1 = j   ; k1 = k-1
        i2 = i   ; j2 = j   ; k2 = k+1
    end select

  end subroutine jacobi_stencil


end module ICE_Lib_IRS
