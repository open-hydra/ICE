!> Unit test of the dispersion cap of near-empty IG/AG cells (cons_2_prim_IG, cons_2_prim_AG): a cell at or below
!> rho_empty keeps P/rho <= theta_cap (IG) or its tensor scaled by one factor to the trace 3 rho theta_cap (AG), a denser
!> cell is untouched, the temperature is the one of the conservative variables, a NaN dispersion stays a NaN (the cap is
!> a test, not min()), and a negative density leaves every primitive finite (the floor overwrites the tensor).
!> Prints every assertion and exits non-zero if any failed.
program test_near_empty_cap
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
  use ICE_Global_m,       only: rho_empty, theta_cap
  use ICE_Lib_IG,         only: prim_2_cons_IG, cons_2_prim_IG
  use ICE_Lib_AG,         only: prim_2_cons_AG, cons_2_prim_AG
  use ICE_Config_Types_m, only: condensed_phase_t
  implicit none
  real(R8), parameter :: cl = 1250._R8, T0 = 300._R8
  real(R8), parameter :: eps = 1e-25                 ! the conversions' floor, the same default-real literal
  type(condensed_phase_t) :: mat
  integer :: nfail = 0

  mat%cs_al = cl; mat%rho_al = 2950._R8

  call test_ig()
  call test_ag()
  call test_nan()
  call test_negative()

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_near_empty_cap: ', nfail, ' FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_near_empty_cap: PASS'

contains

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


  !> The temperature of the conservative variables, as cons_2_prim forms it, to a few ulps
  logical function same_T(T, Tref)
    real(R8), intent(in) :: T, Tref
    same_T = abs(T - Tref) <= 8._R8*spacing(Tref)
  end function same_T


  subroutine test_ig()
    real(R8) :: rho, p(7), c(7), q(7)

    ! (a) 1e-7 kg/m3 at P/rho = 10: P = rho theta_cap exactly, T that of cons
    rho = 1.e-7_R8
    p = [rho, 1._R8, 0.5_R8, 0._R8, 10._R8*rho, T0, 1.e3_R8]
    c = prim_2_cons_IG(p, mat)
    q = cons_2_prim_IG(c, mat)
    call check(q(5) == q(1)*theta_cap, 'IG rho 1e-7, P/rho 10: P = rho theta_cap exactly')
    call check(same_T(q(6), (c(6) - 0.5_R8*c(5))/(cl*c(1) + eps)) .and. same_T(q(6), T0), &
               'IG rho 1e-7: T from the conservative variables, 300 K')
    call check(q(1) == rho .and. q(2) == c(2)/(c(1) + eps), 'IG rho 1e-7: density and velocity untouched')

    ! (a') at rho_empty itself: capped (rho <= rho_empty)
    p(1) = rho_empty; p(5) = 10._R8*rho_empty
    q = cons_2_prim_IG(prim_2_cons_IG(p, mat), mat)
    call check(q(5) == q(1)*theta_cap, 'IG rho = rho_empty, P/rho 10: capped')

    ! (b) 1e-5 kg/m3 (above rho_empty) at P/rho = 10: untouched
    p(1) = 1.e-5_R8; p(5) = 10._R8*p(1)
    q = cons_2_prim_IG(prim_2_cons_IG(p, mat), mat)
    call check(abs(q(5) - p(5)) <= 1.e-12_R8*p(5), 'IG rho 1e-5 (> rho_empty), P/rho 10: P untouched')

    ! below the cap: untouched
    p(1) = rho; p(5) = 0.5_R8*rho
    q = cons_2_prim_IG(prim_2_cons_IG(p, mat), mat)
    call check(abs(q(5) - p(5)) <= 1.e-12_R8*p(5), 'IG rho 1e-7, P/rho 0.5 (below theta_cap): P untouched')
  end subroutine test_ig


  subroutine test_ag()
    real(R8) :: rho, p(12), c(12), q(12), tens(6), f
    integer  :: n
    logical  :: all_scaled

    ! (c) 1e-7 kg/m3, tensor with off-diagonal terms, trace 30 rho: every component times 3 rho theta_cap / tr = 0.1
    rho  = 1.e-7_R8
    tens = rho*[10._R8, 2._R8, 1._R8, 12._R8, 3._R8, 8._R8]          ! P11 P12 P13 P22 P23 P33
    p = [rho, 1._R8, 0.5_R8, 0.25_R8, tens, T0, 1.e3_R8]
    c = prim_2_cons_AG(p, mat)
    q = cons_2_prim_AG(c, mat)
    f = 3._R8*rho*theta_cap/(30._R8*rho)
    all_scaled = .true.
    do n = 1, 6
      all_scaled = all_scaled .and. abs(q(4+n) - f*tens(n)) <= 1.e-12_R8*abs(tens(n))
    enddo
    call check(all_scaled, 'AG rho 1e-7, trace 30 rho: all six components scaled by 0.1 (off-diagonals too)')
    call check(abs(q(5) + q(8) + q(10) - 3._R8*rho*theta_cap) <= 1.e-12_R8*rho, &
               'AG rho 1e-7: trace 3 rho theta_cap after the cap')
    call check(same_T(q(11), (c(11) - 0.5_R8*(c(5) + c(8) + c(10)))/(cl*c(1) + eps)) .and. same_T(q(11), T0), &
               'AG rho 1e-7: T from the conservative variables, 300 K')

    ! 1e-5 kg/m3: untouched
    p(1) = 1.e-5_R8; p(5:10) = 1.e-5_R8*[10._R8, 2._R8, 1._R8, 12._R8, 3._R8, 8._R8]
    q = cons_2_prim_AG(prim_2_cons_AG(p, mat), mat)
    call check(all(abs(q(5:10) - p(5:10)) <= 1.e-12_R8*1.e-4_R8), 'AG rho 1e-5 (> rho_empty): tensor untouched')
  end subroutine test_ag


  !> (d) a NaN dispersion stays NaN: the cap is a test, which a NaN fails, so check_prim sees it
  subroutine test_nan()
    real(R8) :: nan, rho, c(7), q(7), ca(12), qa(12)
    nan = ieee_value(1._R8, ieee_quiet_nan)
    rho = 1.e-7_R8
    c = prim_2_cons_IG([rho, 1._R8, 0.5_R8, 0._R8, 10._R8*rho, T0, 1.e3_R8], mat)
    c(5) = nan
    q = cons_2_prim_IG(c, mat)
    call check(ieee_is_nan(q(5)), 'IG rho 1e-7, NaN second moment: P stays NaN')
    ca = prim_2_cons_AG([rho, 1._R8, 0.5_R8, 0.25_R8, rho*[10._R8, 2._R8, 1._R8, 12._R8, 3._R8, 8._R8], T0, 1.e3_R8], &
                        mat)
    ca(5) = nan
    qa = cons_2_prim_AG(ca, mat)
    call check(ieee_is_nan(qa(5)), 'AG rho 1e-7, NaN second moment: P11 stays NaN')
  end subroutine test_nan


  !> (e) rho < 0 (the cap's test fires with a negative factor): the floor overwrites every primitive, all finite
  subroutine test_negative()
    real(R8) :: c(12), q(12)
    c = prim_2_cons_AG([1.e-7_R8, 1._R8, 0.5_R8, 0.25_R8, 1.e-7_R8*[10._R8, 2._R8, 1._R8, 12._R8, 3._R8, 8._R8], &
                        T0, 1.e3_R8], mat)
    c(1) = -1.e-10_R8
    q = cons_2_prim_AG(c, mat)
    call check(all(ieee_is_finite(q)) .and. all(q(5:10) >= 0._R8) .and. q(6) == 0._R8 .and. q(1) > 0._R8, &
               'AG rho < 0: every primitive finite, the tensor at the floor')
  end subroutine test_negative

end program test_near_empty_cap
