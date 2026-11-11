!> Unit test of the bounds of the MUSCL slopes (Lib_Reconstruction, Lib_Limiters): on a line of cells of any lengths the
!> value of a cell at each of its faces lies between the value of the cell and that of the neighbour across the face.
!>  B1 a cell 1.25 times its neighbour (values 1, 0.01, 0, 0 on lengths 1.25, 1.25, 1, 1): the face value stays in
!>     [0, 0.01] for every bounded limiter (without the bounds van Leer gives -9.88e-4);
!>  B2 random lines of five cells, lengths in ratios up to 100: every face value within its two cells, and the two faces
!>     of the middle cell with the same slope (monotone data for every limiter, data of any sign for those zero at r <= 0);
!>  B3 equal cells: the faces of the limiter's own slope, bit for bit;
!>  B4 data linear along a line of unequal cells: the exact face values, no bound acting;
!>  B5 a side whose cell value sits on the face (dll or dlr = 0, an injection ghost): the cell value, finite.
!> The variable reconstructed is a velocity component of an MK state, so check_prim never sends a face back to first
!> order. Prints every assertion and exits non-zero if any failed.
program test_muscl_bounds
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ICE_Lib_MK,             only: check_prim_MK
  use ICE_Lib_Model,          only: check_prim
  use ICE_Lib_Limiters,       only: assign_limiter, limiter
  use ICE_Lib_Reconstruction, only: state_reconstruction
  use ICE_Config_Types_m,     only: condensed_phase_t
  implicit none
  character(len=10), parameter :: zero_at_neg(8) = [character(len=10) :: 'vanleer', 'mc', 'superbee', 'koren', &
                                                    'umist', 'osher', 'sweby', 'minmod']
  character(len=10), parameter :: all_lim(10) = [character(len=10) :: zero_at_neg, 'vanalbada', 'ospre']
  integer, parameter  :: nlines = 20000
  type(condensed_phase_t) :: mat
  integer :: nfail = 0

  mat%cs_al = 1000._R8; mat%rho_al = 1000._R8
  check_prim => check_prim_MK
  call seed()

  call test_stencil()
  call test_random()
  call test_equal()
  call test_linear()
  call test_injection()

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_muscl_bounds: ', nfail, ' FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_muscl_bounds: PASS'

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

  subroutine seed()
    integer :: n
    integer, allocatable :: s(:)
    call random_seed(size=n)
    allocate(s(n)); s = 20261110
    call random_seed(put=s)
  end subroutine seed

  !> An MK state (rho, u, v, w, T, n) whose u is q
  pure function state(q) result(c)
    real(R8), intent(in) :: q
    real(R8) :: c(6)
    c = [1._R8, q, 0._R8, 0._R8, 300._R8, 1.e8_R8]
  end function state

  !> The face values of the face between cells 2 and 3 of four cells of lengths L and values q
  subroutine faces(q, L, ul, ur, beta)
    real(R8), intent(in)  :: q(4), L(4)
    real(R8), intent(out) :: ul, ur
    real(R8), intent(in), optional :: beta
    real(R8) :: pl(6), pr(6), b
    b = 1._R8; if (present(beta)) b = beta
    call state_reconstruction(state(q(1)), state(q(2)), state(q(3)), state(q(4)), 0.5_R8*(L(1)+L(2)), &
                              0.5_R8*(L(2)+L(3)), 0.5_R8*(L(3)+L(4)), 0.5_R8*L(2), 0.5_R8*L(3), pl, pr, b, mat, .false.)
    ul = pl(2); ur = pr(2)
  end subroutine faces

  !> v lies between a and b, to 8 units in the last place of the larger
  pure logical function between(v, a, b)
    real(R8), intent(in) :: v, a, b
    real(R8) :: tol
    tol = 8._R8*spacing(max(abs(a), abs(b), tiny(1._R8)))
    between = v >= min(a, b) - tol .and. v <= max(a, b) + tol
  end function between

  subroutine test_stencil()
    real(R8) :: ul, ur
    integer  :: m
    logical  :: ok
    character(len=200) :: msg
    do m = 1, size(all_lim)
      call assign_limiter(trim(all_lim(m)))
      call faces([1._R8, 0.01_R8, 0._R8, 0._R8], [1.25_R8, 1.25_R8, 1._R8, 1._R8], ul, ur)
      ok = between(ul, 0.01_R8, 0._R8) .and. between(ur, 0.01_R8, 0._R8)
      write(msg,'(A,A,A,ES11.3,A,ES11.3,A)') 'B1 ', trim(all_lim(m)), ': a cell 1.25 times its neighbour, faces ', ul, &
        ' and ', ur, ' in [0, 0.01]'
      call check(ok, trim(msg))
    enddo
  end subroutine test_stencil

  subroutine test_random()
    real(R8) :: L(5), q(5), r(10), sg(5), ul, ur, vl, vr, jump
    integer  :: m, n, bad_face, bad_slope, k
    logical  :: monotone
    character(len=200) :: msg
    do m = 1, size(all_lim)
      call assign_limiter(trim(all_lim(m)))
      bad_face = 0; bad_slope = 0
      do n = 1, nlines
        call random_number(r); call random_number(sg)
        L = 10._R8**(2._R8*r(1:5))                        ! lengths 1 to 100
        q = 10._R8**(6._R8*r(6:10) - 3._R8)               ! values over six decades
        monotone = m > size(zero_at_neg) .or. mod(n, 2) == 0
        if (monotone) then
          call sort(q)
        else
          do k = 1, 5
            if (sg(k) > 0.5_R8) q(k) = -q(k)
          enddo
        endif
        call faces(q(1:4), L(1:4), ul, ur)                ! face 2|3
        call faces(q(2:5), L(2:5), vl, vr)                ! face 3|4
        if (.not. (between(ul, q(2), q(3)) .and. between(ur, q(2), q(3)) .and. &
                   between(vl, q(3), q(4)) .and. between(vr, q(3), q(4)))) bad_face = bad_face + 1
        ! cell 3: ur is its value at face 2|3, vl at face 3|4, both L(3)/2 from its centre
        jump = max(abs(q(3) - q(2)), abs(q(4) - q(3)), tiny(1._R8))
        if (abs((vl - q(3)) - (q(3) - ur)) > max(1.e-12_R8*jump, 4._R8*spacing(abs(q(3))))) bad_slope = bad_slope + 1
      enddo
      write(msg,'(A,A,A,I0,A,I0,A,I0,A)') 'B2 ', trim(all_lim(m)), ': ', nlines, ' random lines, ', bad_face, &
        ' face values outside their cells, ', bad_slope, ' cells with two slopes'
      call check(bad_face == 0 .and. bad_slope == 0, trim(msg))
    enddo
  end subroutine test_random

  subroutine sort(q)
    real(R8), intent(inout) :: q(:)
    integer  :: i, j
    real(R8) :: t
    do i = 2, size(q)
      t = q(i); j = i - 1
      do while (j >= 1)
        if (q(j) <= t) exit
        q(j+1) = q(j); j = j - 1
      enddo
      q(j+1) = t
    enddo
  end subroutine sort

  subroutine test_equal()
    real(R8), parameter :: big = huge(1._R8)
    real(R8) :: q(4), r(4), ul, ur, el, er, h, s0, s1, s2
    integer  :: m, n, bad
    character(len=200) :: msg
    do m = 1, size(all_lim)
      call assign_limiter(trim(all_lim(m)))
      bad = 0
      do n = 1, nlines
        call random_number(r); call random_number(h)
        q = 2._R8*r - 1._R8; h = 0.01_R8 + h
        call faces(q, [h, h, h, h], ul, ur)
        s0 = (q(2) - q(1))/(0.5_R8*(h + h)); s1 = (q(3) - q(2))/(0.5_R8*(h + h)); s2 = (q(4) - q(3))/(0.5_R8*(h + h))
        el = q(2) + (limiter(s1, s0, big, big)*1._R8)*(0.5_R8*h)
        er = q(3) - (limiter(s2, s1, big, big)*1._R8)*(0.5_R8*h)
        if (ul /= el .or. ur /= er) bad = bad + 1
      enddo
      write(msg,'(A,A,A,I0,A)') 'B3 ', trim(all_lim(m)), ': equal cells, ', bad, ' faces off the unbounded slope'
      call check(bad == 0, trim(msg))
    enddo
  end subroutine test_equal

  subroutine test_linear()
    real(R8) :: L(4), x(4), r(4), ul, ur, xf, a
    integer  :: m, n, bad
    character(len=200) :: msg
    do m = 1, size(all_lim)
      call assign_limiter(trim(all_lim(m)))
      bad = 0
      do n = 1, 2000
        call random_number(r); call random_number(a)
        L = 10._R8**(2._R8*r); a = 2._R8*a - 1._R8
        x(1) = 0.5_R8*L(1)
        x(2) = x(1) + 0.5_R8*(L(1) + L(2)); x(3) = x(2) + 0.5_R8*(L(2) + L(3)); x(4) = x(3) + 0.5_R8*(L(3) + L(4))
        xf = x(2) + 0.5_R8*L(2)
        call faces(a*x, L, ul, ur)
        if (abs(ul - a*xf) > 1.e-12_R8*abs(a)*x(4) .or. abs(ur - a*xf) > 1.e-12_R8*abs(a)*x(4)) bad = bad + 1
      enddo
      write(msg,'(A,A,A,I0,A)') 'B4 ', trim(all_lim(m)), ': linear data on unequal cells, ', bad, ' inexact faces'
      call check(bad == 0, trim(msg))
    enddo
  end subroutine test_linear

  subroutine test_injection()
    real(R8) :: pl(6), pr(6)
    integer  :: m
    logical  :: ok
    ok = .true.
    do m = 1, size(all_lim)
      call assign_limiter(trim(all_lim(m)))
      ! the left cell's value on the face (dll = 0, dl1 = half the right cell), and the right one's (dlr = 0)
      call state_reconstruction(state(0._R8), state(1._R8), state(3._R8), state(4._R8), 1._R8, 0.5_R8, 1._R8, &
                                0._R8, 0.5_R8, pl, pr, 1._R8, mat, .false.)
      ok = ok .and. pl(2) == 1._R8 .and. all(ieee_is_finite(pl)) .and. all(ieee_is_finite(pr)) .and. &
           between(pr(2), 1._R8, 3._R8)
      call state_reconstruction(state(0._R8), state(1._R8), state(3._R8), state(4._R8), 1._R8, 0.5_R8, 1._R8, &
                                0.5_R8, 0._R8, pl, pr, 1._R8, mat, .false.)
      ok = ok .and. pr(2) == 3._R8 .and. all(ieee_is_finite(pl)) .and. all(ieee_is_finite(pr)) .and. &
           between(pl(2), 1._R8, 3._R8)
    enddo
    call check(ok, 'B5 a side whose value sits on the face (dll or dlr = 0) keeps it, finite, for every limiter')
  end subroutine test_injection

end program test_muscl_bounds
