!> Unit test of the source a coupler sets on a family (cond_phase(p)%ext, ICE_Lib_Residual's add_external_source):
!> one block of 3 x 2 x 1 cells, family 1 of 6 slots (MK) and family 2 of 8 (a solidifying MK), a residual and a
!> per-cell dt in every slot.
!>   X1 family 1 with a source: each slot's residual gains ext dt (to 4 ulp: a fused multiply-add may round once)
!>   X2 family 2 with no source: its residual is untouched, bit for bit
!>   X3 a frozen slot (slot 2 of family 1) stays 0, the others gain ext dt
!>   X4 family 2 with a source on its 8 slots: the two slots past the closure's six gain it too
!> Prints every assertion and exits non-zero if any failed.
program test_external_source
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use ICE_Global_m,         only: ngroups, ncond, nbase, ncond_max
  use ICE_Advanced_Types_m, only: ICE_domain_type
  use ICE_Mod_Freeze,       only: nfrozen, frozen
  use ICE_Lib_Residual,     only: add_external_source
  implicit none
  type(ICE_domain_type) :: grid
  real(R8), allocatable :: r1(:,:,:,:), r2(:,:,:,:)
  integer :: nfail = 0, p, i, j, s

  ngroups = 2
  allocate(ncond(2), nbase(2))
  ncond = [6, 8]
  nbase = [6, 6]
  grid%nb = 1
  allocate(grid%blk(1))
  grid%blk(1)%dim = [3, 2, 1]
  allocate(grid%blk(1)%cond_phase(2))
  do p = 1, 2
    allocate(grid%blk(1)%cond_phase(p)%residual(ncond(p), 3, 2, 1), grid%blk(1)%cond_phase(p)%dt(3, 2, 1))
    do j = 1, 2
      do i = 1, 3
        grid%blk(1)%cond_phase(p)%dt(i,j,1) = 1.e-6_R8*real(i + 3*j, R8)/7._R8
        grid%blk(1)%cond_phase(p)%residual(:,i,j,1) = [(0.1_R8*real(s + 10*i + 100*j, R8)/3._R8, s = 1, ncond(p))]
      enddo
    enddo
  enddo
  r1 = grid%blk(1)%cond_phase(1)%residual
  r2 = grid%blk(1)%cond_phase(2)%residual

  ! X1, X2: family 1 has a source, family 2 none
  allocate(grid%blk(1)%cond_phase(1)%ext(6, 3, 2, 1))
  do j = 1, 2
    do i = 1, 3
      grid%blk(1)%cond_phase(1)%ext(:,i,j,1) = [(real(s, R8)*1.e3_R8/real(i + j, R8), s = 1, 6)]
    enddo
  enddo
  call add_external_source(grid, 1)
  call add_external_source(grid, 2)
  call check(near(grid%blk(1)%cond_phase(1)%residual, r1 + grid%blk(1)%cond_phase(1)%ext &
                 * spread(grid%blk(1)%cond_phase(1)%dt, 1, 6)), 'X1 family 1: residual + ext dt in every slot')
  call check(all(grid%blk(1)%cond_phase(2)%residual == r2), 'X2 family 2 without a source: residual untouched')

  ! X3: slot 2 of family 1 frozen
  allocate(nfrozen(2), frozen(ncond_max, 2))
  nfrozen = [1, 0]
  frozen = 0
  frozen(1, 1) = 2
  grid%blk(1)%cond_phase(1)%residual = r1
  call add_external_source(grid, 1)
  call check(all(grid%blk(1)%cond_phase(1)%residual(2,:,:,:) == 0._R8), 'X3 the frozen slot 2 stays 0')
  call check(near(grid%blk(1)%cond_phase(1)%residual([1,3,4,5,6],:,:,:), r1([1,3,4,5,6],:,:,:) &
                 + grid%blk(1)%cond_phase(1)%ext([1,3,4,5,6],:,:,:) &
                 * spread(grid%blk(1)%cond_phase(1)%dt, 1, 5)), 'X3 the other slots gain ext dt')

  ! X4: family 2 with a source on its 8 slots
  allocate(grid%blk(1)%cond_phase(2)%ext(8, 3, 2, 1))
  grid%blk(1)%cond_phase(2)%ext = 2.5e2_R8
  call add_external_source(grid, 2)
  call check(near(grid%blk(1)%cond_phase(2)%residual, r2 + 2.5e2_R8*spread(grid%blk(1)%cond_phase(2)%dt, 1, 8)), &
             'X4 family 2: all 8 slots gain ext dt, the two past the closure too')

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_external_source: ', nfail, ' FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_external_source: PASS'

contains

  logical function near(a, b)
    real(R8), intent(in) :: a(:,:,:,:), b(:,:,:,:)
    near = all(abs(a - b) <= 4._R8*epsilon(1._R8)*abs(b))
  end function near


  subroutine check(ok, what)
    logical,          intent(in) :: ok
    character(len=*), intent(in) :: what
    if (ok) then
      write(*,'(A)') '   OK   '//what
    else
      write(*,'(A)') '   FAIL '//what
      nfail = nfail + 1
    endif
  end subroutine check

end program test_external_source
