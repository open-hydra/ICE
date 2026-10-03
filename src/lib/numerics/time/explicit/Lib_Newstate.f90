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
    use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, thread_group, cell_slice
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, p, i, j, k
    integer(kind=I4) :: n, tid, g, m, s, i1, j1, k1, i2, j2, k2

    !> Thread groups: each thread copies its slice of its group's blocks (no worksharing:
    !> this routine is orphaned, every thread of the caller's region runs it)
    if (n_tgroups > 1) then
      tid = 0
      !$ tid = omp_get_thread_num()
      call thread_group(tid, g, m, s)
      do n = tg_first(g), tg_first(g+1) - 1
        b = tg_blocks(n)
        call cell_slice(grid%blk(b)%dim, m, s, i1, j1, k1, i2, j2, k2)
        do p = 1, ngroups
        do k = k1, k2
        do j = merge(j1, 1, k == k1), merge(j2, grid%blk(b)%dim(2), k == k2)
        do i = merge(i1, 1, k == k1 .and. j == j1), merge(i2, grid%blk(b)%dim(1), k == k2 .and. j == j2)
            grid%blk(b)%cond_phase(p)%prim_old(:,i,j,k) = grid%blk(b)%cond_phase(p)%prim(:,i,j,k)
        enddo ; enddo ; enddo
        enddo
      enddo
      return
    end if

    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      do p = 1, ngroups
        !$omp do collapse(3) schedule(runtime)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)
          grid%blk(b)%cond_phase(p)%prim_old(:,i,j,k) = grid%blk(b)%cond_phase(p)%prim(:,i,j,k)
        end do ; enddo ; enddo
        !$omp end do nowait
      enddo
    enddo

  end subroutine state_copy


  subroutine state_update(grid, p, srk)
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, thread_group, cell_slice
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4), intent(in)  :: srk
    integer(kind=I4) :: b, i, j, k
    integer(kind=I4) :: n, tid, g, m, s, i1, j1, k1, i2, j2, k2
    logical :: prim_status, bad_state

    bad_state = .false.
    !> Thread groups: each thread updates its slice of its group's blocks
    if (n_tgroups > 1) then
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(b,i,j,k,prim_status,n, tid, g, m, s, i1, j1, k1, i2, j2, k2) &
    !$OMP SHARED(grid,p,srk,bad_state,obj_condensed,mat_of,tg_first,tg_blocks)
      tid = 0
      !$ tid = omp_get_thread_num()
      call thread_group(tid, g, m, s)
      do n = tg_first(g), tg_first(g+1) - 1
        b = tg_blocks(n)
        call cell_slice(grid%blk(b)%dim, m, s, i1, j1, k1, i2, j2, k2)
        do k = k1, k2
        do j = merge(j1, 1, k == k1), merge(j2, grid%blk(b)%dim(2), k == k2)
        do i = merge(i1, 1, k == k1 .and. j == j1), merge(i2, grid%blk(b)%dim(1), k == k2 .and. j == j2)
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
        enddo ; enddo ; enddo
      enddo
    !$OMP END PARALLEL
    else
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(b,i,j,k,prim_status) SHARED(grid,p,srk,bad_state,obj_condensed,mat_of)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$OMP DO COLLAPSE (3) schedule(runtime)
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
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL
    end if
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
