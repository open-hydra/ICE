!> Unit test of the OpenMP schedule ICE sets (Wrap_Setup): each omp_set_schedule, called through omp_lib as ICE
!> calls it, must read back with omp_get_schedule. Run with OMP_SCHEDULE=guided,3, so that a set the runtime
!> ignores reads back as the environment's. The kind is compared without OpenMP 5.0's monotonic bit.
!>   S1 static, chunk 0 (ICE's schedule): kind static
!>   S2 dynamic, 7 (the ICE_OMP_SCHEDULE form): kind dynamic, chunk 7
!>   S3 guided, 5: kind guided, chunk 5
!> Only a gfortran build with -fno-underscoring can fail it: the call then reaches libgomp's C entry, which takes
!> its arguments by value. ifx never had the flag.
!> Prints every assertion and exits non-zero if any failed.
program test_omp_schedule
  use omp_lib
  implicit none
  logical :: ok

  ok = .true.
  call set_and_read(omp_sched_static,  0, .false., 'S1 static, 0', ok)
  call set_and_read(omp_sched_dynamic, 7, .true.,  'S2 dynamic, 7', ok)
  call set_and_read(omp_sched_guided,  5, .true.,  'S3 guided, 5', ok)
  if (.not. ok) error stop 1

contains

  subroutine set_and_read(kind, chunk, with_chunk, what, ok)
    integer(omp_sched_kind), intent(in)    :: kind
    integer,                 intent(in)    :: chunk
    logical,                 intent(in)    :: with_chunk
    character(len=*),        intent(in)    :: what
    logical,                 intent(inout) :: ok
    integer(omp_sched_kind) :: got
    integer :: got_chunk
    character(len=64) :: txt
    call omp_set_schedule(kind, chunk)
    call omp_get_schedule(got, got_chunk)
    got = ibclr(got, 31)
    write(txt, '(A,I0,A,I0)') ': read back kind ', got, ', chunk ', got_chunk
    if (got == kind .and. (.not. with_chunk .or. got_chunk == chunk)) then
      write(*,'(A)') '   OK   '//what//trim(txt)
    else
      write(*,'(A)') '   FAIL '//what//trim(txt)
      ok = .false.
    endif
  end subroutine set_and_read

end program test_omp_schedule
