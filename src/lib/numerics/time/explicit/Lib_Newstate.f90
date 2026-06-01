module ICE_Lib_Newstate
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Global_m
  use ICE_Advanced_Types_m

  implicit none
  private
  public :: state_copy, state_update

contains

  subroutine state_copy(grid)
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, p, i, j, k

    do b = 1, grid%nb
      do p = 1, ngroups
        !$omp do collapse(3)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)
          grid%blk(b)%cond_phase(p)%prim_old(:,i,j,k) = grid%blk(b)%cond_phase(p)%prim(:,i,j,k)
        end do ; enddo ; enddo
        !$omp end do
      enddo
    enddo

  end subroutine state_copy


  subroutine state_update(grid, p, srk)
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4), intent(in)  :: srk
    integer(kind=I4) :: b, i, j, k
    logical :: prim_status

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(b,i,j,k,prim_status) SHARED(grid,p,srk)
    do b = 1, grid%nb

      !$OMP DO COLLAPSE (3)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        call state_update_(grid%blk(b)%cond_phase(p)%prim(:,i,j,k),     &
                           grid%blk(b)%cond_phase(p)%prim_old(:,i,j,k), &
                           grid%blk(b)%cond_phase(p)%residual(:,i,j,k), &
                           srk, prim_status)

        if (.not. prim_status) then
          write(*,'(A,I4,A,3I4,A,I4)') " Error in block ", b, " cell ", i, j, k, " phase ", p
          write(*,*) grid%blk(b)%cond_phase(p)%prim_old(:,i,j,k)
          write(*,*) grid%blk(b)%cond_phase(p)%prim(:,i,j,k)
          write(*,*) grid%blk(b)%cond_phase(p)%residual(:,i,j,k)
          write(*,*) grid%blk(b)%cond_phase(p)%source(:,i,j,k)
          stop
        endif

      end do ; end do ; end do
      !$OMP END DO
    enddo
    !$OMP END PARALLEL

  end subroutine state_update


  subroutine state_update_(prim, prim_old, residual, srk, prim_status)
    use ICE_Lib_RK,    only: RK_stage
    use ICE_Lib_Model
    implicit none
    real(kind=R8), dimension(:), intent(inout) :: prim
    real(kind=R8), dimension(:), intent(in)    :: prim_old, residual
    integer(kind=I4),            intent(in)    :: srk
    logical,                     intent(out)   :: prim_status

    real(kind=R8), dimension(size(prim)) :: cons, cons_old

    cons_old    = prim_2_cons(prim_old)
    cons        = prim_2_cons(prim)
    cons        = RK_stage(srk, nrk, cons, cons_old, residual)
    prim        = cons_2_prim(cons)
    prim_status = check_prim(prim)

  end subroutine state_update_


end module ICE_Lib_Newstate
