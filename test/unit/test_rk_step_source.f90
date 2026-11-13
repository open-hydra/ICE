!> Unit test of RK_accumulate (ICE_Lib_RK), the step's source a two-way coupled gas takes. A cell whose residual is
!> its source alone (no flux): the stages of RK_stage move it by dt times the RK_accumulate sum of the stage sources.
!>   W1 for euler, RK2 and RK3 and stage sources that differ, U_n - U_0 = dt * acc_n to round-off
!>   W2 the weights: acc_n = S_1 (euler), (S_1 + S_2)/2 (RK2), S_1/6 + S_2/6 + 2 S_3/3 (RK3)
!>   W3 with RK2 and RK3 the last stage's source alone misses the change (the coupling's old hand-off)
!> Prints every assertion and exits non-zero if any failed.
program test_rk_step_source
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use ICE_Lib_RK, only: RK_stage, RK_accumulate
  implicit none
  integer,  parameter :: NV = 3
  real(R8), parameter :: DT = 0.37_R8
  real(R8), parameter :: U0(NV) = [1.3_R8, -0.2_R8, 2.5e3_R8]
  real(R8), parameter :: S(NV,3) = reshape([0.7_R8, -1.1_R8, 4.0e2_R8,   &
                                            0.2_R8,  0.5_R8, 1.5e2_R8,   &
                                           -0.4_R8,  0.9_R8, 3.1e2_R8], [NV,3])
  real(R8) :: U(NV), acc(NV), want(NV)
  integer  :: n, k
  logical  :: ok
  character(len=8) :: tag

  ok = .true.
  do n = 1, 3
    write(tag,'(A,I0)') 'n_rk=', n
    U   = U0
    acc = 0._R8
    do k = 1, n
      U   = RK_stage(k, n, U, U0, DT*S(:,k))
      acc = RK_accumulate(k, n, acc, S(:,k))
    end do
    call check(maxval(abs(U - U0 - DT*acc)) <= 1.e-12_R8*maxval(abs(U - U0)), trim(tag)//' W1 U_n - U_0 = dt acc_n', ok)
    select case (n)
    case (1); want = S(:,1)
    case (2); want = (S(:,1) + S(:,2))/2._R8
    case (3); want = S(:,1)/6._R8 + S(:,2)/6._R8 + 2._R8*S(:,3)/3._R8
    end select
    call check(maxval(abs(acc - want)) <= 1.e-13_R8*maxval(abs(want)), trim(tag)//' W2 the weights', ok)
    if (n > 1) call check(maxval(abs(U - U0 - DT*S(:,n))) > 1.e-3_R8*maxval(abs(U - U0)), &
                          trim(tag)//' W3 the last stage alone misses the change', ok)
  end do

  if (.not. ok) then
    write(*,'(A)') 'test_rk_step_source: FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_rk_step_source: PASS'

contains

  subroutine check(cond, what, ok)
    logical,          intent(in)    :: cond
    character(len=*), intent(in)    :: what
    logical,          intent(inout) :: ok
    if (cond) then
      write(*,'(A)') '   OK   '//what
    else
      write(*,'(A)') '   FAIL '//what
      ok = .false.
    endif
  end subroutine check

end program test_rk_step_source
