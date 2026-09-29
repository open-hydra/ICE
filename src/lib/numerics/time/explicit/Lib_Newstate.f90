module ICE_Lib_Newstate
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
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
      if (.not. is_local_block(b)) cycle
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
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4), intent(in)  :: srk
    integer(kind=I4) :: b, i, j, k
    logical :: prim_status, bad_state

    bad_state = .false.
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(b,i,j,k,prim_status) SHARED(grid,p,srk,bad_state,obj_condensed,mat_of)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$OMP DO COLLAPSE (3)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        call state_update_(grid%blk(b)%cond_phase(p)%prim(:,i,j,k),     &
                           grid%blk(b)%cond_phase(p)%prim_old(:,i,j,k), &
                           grid%blk(b)%cond_phase(p)%residual(:,i,j,k), &
                           srk, prim_status, obj_condensed(mat_of(p)))

        if (.not. prim_status) then
          !$OMP CRITICAL (ICE_state_report)
          if (.not. bad_state) then                     !> dump the first failing cell only
            write(*,'(A,I4,A,3I4,A,I4)') " Error in block ", b, " cell ", i, j, k, " phase ", p
            write(*,*) grid%blk(b)%cond_phase(p)%prim_old(:,i,j,k)
            write(*,*) grid%blk(b)%cond_phase(p)%prim(:,i,j,k)
            write(*,*) grid%blk(b)%cond_phase(p)%residual(:,i,j,k)
            write(*,*) grid%blk(b)%cond_phase(p)%source(:,i,j,k)
          endif
          bad_state = .true.
          !$OMP END CRITICAL (ICE_state_report)
        endif

      end do ; end do ; end do
      !$OMP END DO
    enddo
    !$OMP END PARALLEL
    if (bad_state) then
      write(*,'(A)') ' [ERROR] [ICE::state_update] invalid state after the update (negative density or NaN)'
      error stop 1
    endif

  end subroutine state_update


  subroutine state_update_(prim, prim_old, residual, srk, prim_status, mat)
    use ICE_Lib_RK,    only: RK_stage
    use ICE_Lib_Model
    use ICE_Config_Types_m, only: condensed_phase_t
    implicit none
    real(kind=R8), dimension(:), intent(inout) :: prim
    real(kind=R8), dimension(:), intent(in)    :: prim_old, residual
    integer(kind=I4),            intent(in)    :: srk
    logical,                     intent(out)   :: prim_status
    type(condensed_phase_t),     intent(in)    :: mat

    real(kind=R8), dimension(size(prim)) :: cons, cons_old

    cons_old    = prim_2_cons(prim_old, mat)
    cons        = prim_2_cons(prim, mat)
    cons        = RK_stage(srk, nrk, cons, cons_old, residual)
    prim        = cons_2_prim(cons, mat)
    prim_status = check_prim(prim)

  end subroutine state_update_


end module ICE_Lib_Newstate
