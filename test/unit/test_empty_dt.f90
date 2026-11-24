!> Unit test of the step of the empty cells of a steady family that births fill (ICE_Mod_dt's bound_empty_dt): one
!> block of 4 x 1 x 1 cells of 1 mm, family 1 (MK, 6 slots), cfl 0.4, the steps as compute_dt leaves them: cfl dx/|u|
!> in an occupied cell, the dt-max ceiling (1e-4 s) in an empty one, whose velocity bounds nothing.
!>   E1 no source: every step untouched, bit for bit
!>   E2 births of ICE's own (births_gen) at 800 m/s in empty cell 3: both empty cells step at most cfl dx/800 (to 4
!>      ulp), the occupied cells keep theirs, bit for bit
!>   E3 the same births set by a coupler (imports): the same steps
!>   E4 no occupied cell (a generation family before its first births): every cell steps at the source's cfl dx/800
!>   E5 a source in an occupied cell only: the empty cells, which the cells it fills fill, step at most the occupied
!>      cells' smallest step, bit for bit
!>   E6 the same, empty cell 4 ten times narrower (0.1 mm): it steps at most its own cfl dx/100 at the occupied cells'
!>      speed (to 4 ulp), below their smallest step; empty cell 3 keeps that step, bit for bit
!>   E7 time-accurate, the births as imports in empty cell 3, no occupied cell: the global step (dt-max, as compute_dt
!>      leaves it) lowered to the source's cfl dx/800 (to 4 ulp)
!> Prints every assertion and exits non-zero if any failed.
program test_empty_dt
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use ICE_Global_m,         only: ngroups, ncond, nbase, ndir
  use ICE_Advanced_Types_m, only: ICE_domain_type
  use ICE_Config_Types_m,   only: obj_time_scheme
  use ICE_Lib_Model,        only: assign_wavespeed_make
  use ICE_Mod_dt,           only: bound_empty_dt
  implicit none
  type(ICE_domain_type) :: grid
  real(R8), parameter :: dx = 1.e-3_R8, cfl = 0.4_R8, dtmax = 1.e-4_R8
  real(R8) :: occ, src, dt0(4,1,1), s(6)
  integer  :: nfail = 0, i

  ngroups = 1
  ndir    = 1
  allocate(ncond(1), nbase(1))
  ncond = 6
  nbase = 6
  allocate(obj_time_scheme%model(1))
  obj_time_scheme%model(1) = 'MK'
  obj_time_scheme%time_accurate = .false.
  call assign_wavespeed_make(1)
  grid%nb = 1
  allocate(grid%blk(1))
  grid%blk(1)%dim = [4, 1, 1]
  allocate(grid%blk(1)%M(4,1,1))
  do i = 1, 4
    grid%blk(1)%M(i,1,1)%c      = 0._R8
    grid%blk(1)%M(i,1,1)%c(1,1) = 1._R8/dx
  enddo
  allocate(grid%blk(1)%cond_phase(1))
  allocate(grid%blk(1)%cond_phase(1)%prim(6,4,1,1), grid%blk(1)%cond_phase(1)%dt(4,1,1))
  occ = cfl*dx/100._R8
  src = cfl*dx/800._R8
  s   = [2.e2_R8, 2.e2_R8*800._R8, 0._R8, 0._R8, 0._R8, 0._R8]

  ! cells 1 and 2 occupied at 100 m/s, 3 and 4 empty (ICE's empty state)
  call occupy(2)

  ! E1: no source
  call bound_empty_dt(1, cfl, grid)
  call check(all(grid%blk(1)%cond_phase(1)%dt == dt0), 'E1 no source: every step untouched')

  ! E2: a source at 800 m/s in empty cell 3
  allocate(grid%blk(1)%cond_phase(1)%births_gen(6,4,1,1))
  grid%blk(1)%cond_phase(1)%births_gen = 0._R8
  grid%blk(1)%cond_phase(1)%births_gen(:,3,1,1) = s
  call bound_empty_dt(1, cfl, grid)
  call check(all(grid%blk(1)%cond_phase(1)%dt(1:2,1,1) == dt0(1:2,1,1)), 'E2 the occupied cells keep their steps')
  call check(near(grid%blk(1)%cond_phase(1)%dt(3:4,1,1), src), 'E2 the empty cells step at most cfl dx/800')

  ! E3: the same births set by a coupler
  deallocate(grid%blk(1)%cond_phase(1)%births_gen)
  allocate(grid%blk(1)%cond_phase(1)%imports(6,4,1,1))
  grid%blk(1)%cond_phase(1)%imports = 0._R8
  grid%blk(1)%cond_phase(1)%imports(:,3,1,1) = s
  call occupy(2)
  call bound_empty_dt(1, cfl, grid)
  call check(all(grid%blk(1)%cond_phase(1)%dt(1:2,1,1) == dt0(1:2,1,1)) .and. &
             near(grid%blk(1)%cond_phase(1)%dt(3:4,1,1), src), 'E3 imports: the same steps')

  ! E4: no occupied cell
  call occupy(0)
  call bound_empty_dt(1, cfl, grid)
  call check(near(grid%blk(1)%cond_phase(1)%dt(:,1,1), src), 'E4 no occupied cell: every cell at the source''s step')

  ! E5: a source in occupied cell 1 only
  deallocate(grid%blk(1)%cond_phase(1)%imports)
  allocate(grid%blk(1)%cond_phase(1)%births_gen(6,4,1,1))
  grid%blk(1)%cond_phase(1)%births_gen = 0._R8
  grid%blk(1)%cond_phase(1)%births_gen(:,1,1,1) = s
  call occupy(2)
  call bound_empty_dt(1, cfl, grid)
  call check(all(grid%blk(1)%cond_phase(1)%dt(:,1,1) == occ), &
             'E5 a source in an occupied cell: the empty cells at the occupied cells'' step')

  ! E6: the same, empty cell 4 ten times narrower
  grid%blk(1)%M(4,1,1)%c(1,1) = 10._R8/dx
  call occupy(2)
  call bound_empty_dt(1, cfl, grid)
  call check(all(grid%blk(1)%cond_phase(1)%dt(1:3,1,1) == occ) .and. &
             near(grid%blk(1)%cond_phase(1)%dt(4:4,1,1), occ/10._R8), &
             'E6 a narrower empty cell steps at most its own width at the occupied cells'' speed')
  grid%blk(1)%M(4,1,1)%c(1,1) = 1._R8/dx

  ! E7: time-accurate, the births as imports in empty cell 3, no occupied cell
  deallocate(grid%blk(1)%cond_phase(1)%births_gen)
  allocate(grid%blk(1)%cond_phase(1)%imports(6,4,1,1))
  grid%blk(1)%cond_phase(1)%imports = 0._R8
  grid%blk(1)%cond_phase(1)%imports(:,3,1,1) = s
  call occupy(0)
  obj_time_scheme%time_accurate = .true.
  grid%dtglobal = dtmax
  call bound_empty_dt(1, cfl, grid)
  call check(near([grid%dtglobal], src), 'E7 time-accurate: the global step at the source''s step')
  obj_time_scheme%time_accurate = .false.

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_empty_dt: ', nfail, ' FAILED'
    write(*,'(A,4ES12.4)') '  steps: ', grid%blk(1)%cond_phase(1)%dt(:,1,1)
    error stop 1
  endif
  write(*,'(A)') 'test_empty_dt: PASS'

contains

  !> the first n cells occupied at 100 m/s, the others empty; the steps as compute_dt leaves them (dt0)
  subroutine occupy(n)
    integer, intent(in) :: n
    integer :: c
    grid%blk(1)%cond_phase(1)%prim = 1.e-25_R8
    do c = 1, 4
      dt0(c,1,1) = dtmax
      if (c <= n) then
        grid%blk(1)%cond_phase(1)%prim(:,c,1,1) = [1._R8, 100._R8, 0._R8, 0._R8, 300._R8, 1.e9_R8]
        dt0(c,1,1) = occ
      endif
    enddo
    grid%blk(1)%cond_phase(1)%dt = dt0
  end subroutine occupy


  logical function near(a, b)
    real(R8), intent(in) :: a(:), b
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

end program test_empty_dt
