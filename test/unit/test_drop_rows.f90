!> Unit test of Drop_Rows_After (ICE_IO_Probes), the rows of a history file that a restart drops, on files written
!> here with the layout of residual-history.dat (an iteration, a time, a value per row) or of a probe file (an
!> iteration or a time, then a value):
!>   D1 a run gone past its solution: iterations 1..8, the solution at 5: 1..5 stay
!>   D2 a count that started again (a restart from a solution without the counters): 1..6 then 1..4, the solution
!>      at 3 of the second run: 1..6 and 1..3 stay
!>   D3 a solution older than every row: 4..6, the solution at 2: no row stays
!>   D4 a line that starts with '#' stays, and a last row cut short goes: 1..3, then a row of 2 values of 3
!>   D5 time rows among iteration rows (a probe), the solution at t = 3e-3: the first time row after it goes with
!>      every row after it, the iteration rows before it stay
!>   D6 the iteration rows of the same file, the solution at iteration 8: the rows from iteration 12 go
!> Prints every assertion and exits non-zero if any failed.
program test_drop_rows
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use ICE_IO_Probes, only: Drop_Rows_After
  implicit none
  character(len=*), parameter :: probe(8) = [character(len=24) :: &
    '           4   1.0', '   1.2E-003   1.0', '           8   1.0', '   2.4E-003   1.0', &
    '          12   1.0', '   4.0E-003   1.0', '   4.6E-003   1.0', '          16   1.0']
  logical :: ok

  ok = .true.

  call write_rows('d1.dat', [1, 2, 3, 4, 5, 6, 7, 8])
  call Drop_Rows_After('d1.dat', 3, 1, 5._R8, times=.false.)
  call check(near(firsts('d1.dat'), [1._R8, 2._R8, 3._R8, 4._R8, 5._R8]), &
             'D1 a run gone past its solution at 5: of 1..8, 1..5 stay', ok)

  call write_rows('d2.dat', [1, 2, 3, 4, 5, 6, 1, 2, 3, 4])
  call Drop_Rows_After('d2.dat', 3, 1, 3._R8, times=.false.)
  call check(near(firsts('d2.dat'), [1._R8, 2._R8, 3._R8, 4._R8, 5._R8, 6._R8, 1._R8, 2._R8, 3._R8]), &
             'D2 a count that started again, the solution at 3: of 1..6 1..4, 1..6 1..3 stay', ok)

  call write_rows('d3.dat', [4, 5, 6])
  call Drop_Rows_After('d3.dat', 3, 1, 2._R8, times=.false.)
  call check(size(firsts('d3.dat')) == 0, 'D3 a solution at 2, older than every row of 4..6: no row stays', ok)

  call write_lines('d4.dat', [character(len=48) :: '# iteration time value', &
    '       1    0.1000000000E-02    0.5000000000E+00', '       2    0.2000000000E-02    0.5000000000E+00', &
    '       3    0.3000000000E-02    0.5000000000E+00', '       4    0.4000000000E-02'])
  call Drop_Rows_After('d4.dat', 3, 1, 10._R8, times=.false.)
  call check(near(firsts('d4.dat'), [1._R8, 2._R8, 3._R8]) .and. comments('d4.dat') == 1, &
             'D4 the # line stays and the last row, cut short, goes', ok)

  call write_lines('d5.dat', probe)
  call Drop_Rows_After('d5.dat', 2, 1, 3.e-3_R8, times=.true.)
  call check(near(firsts('d5.dat'), [4._R8, 1.2e-3_R8, 8._R8, 2.4e-3_R8, 12._R8]), &
             'D5 time rows, the solution at 3e-3: from the time row 4.0e-3 on, every row goes', ok)

  call write_lines('d6.dat', probe)
  call Drop_Rows_After('d6.dat', 2, 1, 8._R8, times=.false.)
  call check(near(firsts('d6.dat'), [4._R8, 1.2e-3_R8, 8._R8, 2.4e-3_R8]), &
             'D6 iteration rows, the solution at 8: from the iteration row 12 on, every row goes', ok)

  if (.not. ok) then
    write(*,'(A)') 'test_drop_rows: FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_drop_rows: PASS'

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


  !> Rows of a residual history: an iteration, a time, a value
  subroutine write_rows(file, its)
    character(len=*), intent(in) :: file
    integer,          intent(in) :: its(:)
    integer :: u, k
    open(newunit=u, file=file, status='replace', action='write')
    do k = 1, size(its)
      write(u, '(I8,2E20.10)') its(k), 1.e-3_R8 * k, 0.5_R8
    end do
    close(u)
  end subroutine write_rows


  subroutine write_lines(file, lines)
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: lines(:)
    integer :: u, k
    open(newunit=u, file=file, status='replace', action='write')
    do k = 1, size(lines)
      write(u, '(A)') trim(lines(k))
    end do
    close(u)
  end subroutine write_lines


  !> The first value of every line that holds values: neither blank nor a '#' line
  function firsts(file) result(out)
    character(len=*), intent(in) :: file
    real(R8), allocatable :: out(:)
    character(len=256) :: line
    integer  :: u, ios
    real(R8) :: v
    allocate(out(0))
    open(newunit=u, file=file, status='old', action='read')
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(verify(line, ' '):verify(line, ' ')) == '#') cycle
      read(line, *) v
      out = [out, v]
    end do
    close(u)
  end function firsts


  integer function comments(file)
    character(len=*), intent(in) :: file
    character(len=256) :: line
    integer :: u, ios
    comments = 0
    open(newunit=u, file=file, status='old', action='read')
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(verify(line, ' '):verify(line, ' ')) == '#') comments = comments + 1
    end do
    close(u)
  end function comments


  logical function near(a, b)
    real(R8), intent(in) :: a(:), b(:)
    near = size(a) == size(b)
    if (near) near = all(abs(a - b) <= 1.e-12_R8 * max(1._R8, abs(b)))
  end function near

end program test_drop_rows
