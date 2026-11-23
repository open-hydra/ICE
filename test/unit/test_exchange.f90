!> Unit test of ICE's exchange with a coupler (ICE_Mod_Exchange) and of ICE's own generated mass, which takes the same
!> births route (ICE_Mod_Sources' gen_birth, ICE_Lib_Residual's add_births). One block of 2 x 1 x 1 cells and four
!> materials, each with constant properties:
!>   material 1, family 1: an MK family that burns (q-comb 9.53 MJ/kg, d-gen 1 um)
!>   material 2, family 2: an inert MK family, where the generated mass is born
!>   material 3, family 3: an inert IG family
!>   material 4, families 4 and 5: an inert MK material with two populations
!>   G1 set_gas: the gas state lands in the cell, the other cell untouched
!>   G2 exchange_rate: what the gas receives is what the drops lose less what is born elsewhere; the energy comes from
!>      the closure's energy slot (MK slot 5, IG slot 6); the heat of reaction is the family's
!>   G3 take_removed: the record of the cell, with the heat of reaction of all of its mass, and the record emptied
!>   B1 birth_packet: m, m u, m (cp T + |u|^2/2), m / (pi/6 rho d^3)
!>   B2 add_birth: the cell's conservative state gains the packet (to round-off), the other cell untouched
!>   B3 add_birth_rate after clear_birth_rates, twice: the stage's residual gains the sum times dt in slots 1:4, 5
!>      and 6; a frozen slot stays 0; a family with no births is untouched, bit for bit
!>   B4 birth_family: the family of an inert MK material; none for a burning material, an IG family, a material with
!>      two families or an unknown name, each with its reason
!>   B5 set_generation: taken for the burning material; refused, with nothing set, for an inert one or without d-gen
!>   B6 gen_birth, time-accurate: the generation family gains dt times the burning family's gen_due, which is cleared
!>   B7 gen_birth, steady: the generation family's births_gen is the burning family's gen_due, which is cleared
!> Prints every assertion and exits non-zero if any failed.
program test_exchange
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use ICE_Global_m,         only: ngroups, ncond, nbase, ncond_max, mat_of, solid_of
  use ICE_Advanced_Types_m, only: ICE_domain_type
  use ICE_Config_Types_m,   only: obj_condensed, obj_time_scheme
  use ICE_Mod_Freeze,       only: nfrozen, frozen
  use ICE_Lib_MK,           only: prim_2_cons_MK
  use ICE_Lib_Residual,     only: add_births
  use ICE_Mod_Sources,      only: gen_birth
  use ICE_Mod_Exchange
  implicit none
  real(R8), parameter :: PI = acos(-1._R8), RHO = 4000._R8, CP = 1000._R8, QCOMB = 9.53e6_R8
  type(ICE_domain_type) :: grid
  type(exchange_t)      :: x
  character(len=:), allocatable :: why
  real(R8) :: pk(6), pk2(6), c0(6), c1(6), r0(6,2,1,1), r1(7,2,1,1), dt
  integer  :: nfail = 0, p

  call setup()

  ! G1: the gas state
  call set_gas(grid%blk(1), 2, 1, 1, 1.5_R8, [10._R8, -2._R8, 0.5_R8], 2500._R8, 287._R8, 1.2_R8, 1.8e-5_R8, 0.06_R8)
  call check(all(grid%blk(1)%gas_phase%prim(:,2,1,1) == [1.5_R8, 10._R8, -2._R8, 0.5_R8, 2500._R8]) .and. &
             grid%blk(1)%gas_phase%R(2,1,1) == 287._R8 .and. grid%blk(1)%gas_phase%gam(2,1,1) == 1.2_R8 .and. &
             grid%blk(1)%gas_phase%mu(2,1,1) == 1.8e-5_R8 .and. grid%blk(1)%gas_phase%k(2,1,1) == 0.06_R8, &
             'G1 set_gas: the state lands in the cell')
  call check(all(grid%blk(1)%gas_phase%prim(:,1,1,1) == 0._R8), 'G1 set_gas: the other cell untouched')

  ! G2: the exchange rate of the burning MK family and of the IG family
  grid%blk(1)%cond_phase(1)%source(:,1,1,1)  = [-3._R8, -0.3_R8, 0.6_R8, -0.9_R8, -7.e6_R8, 0._R8]
  grid%blk(1)%cond_phase(1)%gen_rate(:,1,1,1) = [1.5_R8, 0.1_R8, 0.2_R8, -0.4_R8, 2.e6_R8]
  grid%blk(1)%cond_phase(1)%heat(1,1,1)      = 1.1e7_R8
  x = exchange_rate(grid%blk(1), 1, 1, 1, 1)
  call check(x%mass == 3._R8 - 1.5_R8 .and. x%generated == 1.5_R8, &
             'G2 MK: the gas gets the mass lost less the generated')
  call check(all(x%momentum == [0.3_R8 - 0.1_R8, -0.6_R8 - 0.2_R8, 0.9_R8 + 0.4_R8]), &
             'G2 MK: the momentum less that of the generated mass')
  call check(x%energy == 7.e6_R8 - 2.e6_R8 .and. x%heat == 1.1e7_R8, 'G2 MK: energy from slot 5, the heat as it is')
  grid%blk(1)%cond_phase(3)%source(:,1,1,1) = [-1._R8, 2._R8, 3._R8, 4._R8, 50._R8, -60._R8, 0._R8]
  x = exchange_rate(grid%blk(1), 3, 1, 1, 1)
  call check(x%mass == 1._R8 .and. all(x%momentum == [-2._R8, -3._R8, -4._R8]) .and. x%energy == 60._R8 .and. &
             x%generated == 0._R8, 'G2 IG: the energy from slot 6, nothing generated')

  ! G3: the removals
  grid%blk(1)%cond_phase(1)%removed(:,2,1,1) = [0.4_R8, 1._R8, 2._R8, 3._R8, 5.e5_R8, 0.6_R8]
  call take_removed(grid%blk(1), 1, 2, 1, 1, x)
  call check(x%mass == 0.4_R8 .and. all(x%momentum == [1._R8, 2._R8, 3._R8]) .and. x%energy == 5.e5_R8 .and. &
             x%generated == 0.6_R8, 'G3 take_removed: the record of the cell')
  call check(near1(x%heat, QCOMB*1._R8), 'G3 take_removed: the heat of reaction of all the removed mass')
  call check(all(grid%blk(1)%cond_phase(1)%removed == 0._R8), 'G3 take_removed: the record emptied')

  ! B1: a packet
  pk = birth_packet(2, 2._R8, [3._R8, 0._R8, 4._R8], 1000._R8, 1.e-6_R8)
  call check(pk(1) == 2._R8 .and. all(pk(2:4) == [6._R8, 0._R8, 8._R8]) .and. &
             near1(pk(5), 2._R8*(CP*1000._R8 + 12.5_R8)) .and. near1(pk(6), 2._R8/(PI/6._R8*RHO*1.e-18_R8)), &
             'B1 birth_packet: mass, momentum, energy, number')

  ! B2: a birth at once in cell 1 of family 2
  c0 = prim_2_cons_MK(grid%blk(1)%cond_phase(2)%prim(:,1,1,1), obj_condensed(2))
  call add_birth(grid%blk(1), 2, 1, 1, 1, pk)
  c1 = prim_2_cons_MK(grid%blk(1)%cond_phase(2)%prim(:,1,1,1), obj_condensed(2))
  call check(all(abs(c1 - (c0 + pk)) <= 1.e-12_R8*abs(c0 + pk)), 'B2 add_birth: the state gains the packet')
  call check(all(grid%blk(1)%cond_phase(2)%prim(:,2,1,1) == drops()), 'B2 add_birth: the other cell untouched')

  ! B3: births in force
  dt = 1.e-6_R8
  grid%blk(1)%cond_phase(2)%dt = dt
  grid%blk(1)%cond_phase(3)%dt = dt
  grid%blk(1)%cond_phase(2)%residual = 0.25_R8
  grid%blk(1)%cond_phase(3)%residual = 0.75_R8
  r0 = grid%blk(1)%cond_phase(2)%residual
  r1 = grid%blk(1)%cond_phase(3)%residual
  pk2 = [1._R8, 2._R8, 3._R8, 4._R8, 5._R8, 6._R8]
  call clear_birth_rates(grid, 2)
  call add_birth_rate(grid%blk(1), 2, 1, 1, 1, pk)
  call add_birth_rate(grid%blk(1), 2, 1, 1, 1, pk2)
  call add_births(grid, 2)
  call add_births(grid, 3)
  call check(all(abs(grid%blk(1)%cond_phase(2)%residual(:,1,1,1) - (r0(:,1,1,1) + (pk + pk2)*dt)) <= &
                 4._R8*epsilon(1._R8)*abs(r0(:,1,1,1) + (pk + pk2)*dt)), 'B3 the residual gains the births times dt')
  call check(all(grid%blk(1)%cond_phase(2)%residual(:,2,1,1) == r0(:,2,1,1)), 'B3 a cell without births untouched')
  call check(all(grid%blk(1)%cond_phase(3)%residual == r1), 'B3 a family without births untouched')
  allocate(nfrozen(5), frozen(ncond_max, 5))
  nfrozen = 0
  frozen  = 0
  nfrozen(2)   = 1
  frozen(1, 2) = 3
  grid%blk(1)%cond_phase(2)%residual = r0
  call add_births(grid, 2)
  call check(grid%blk(1)%cond_phase(2)%residual(3,1,1,1) == 0._R8, 'B3 a frozen slot stays 0')
  nfrozen = 0
  call clear_birth_rates(grid, 2)
  call check(all(grid%blk(1)%cond_phase(2)%births_ext == 0._R8), 'B3 clear_birth_rates clears')

  ! B4: where a birth goes
  call check(birth_family('B', why) == 2 .and. why == '', 'B4 the family of an inert MK material')
  p = birth_family('A', why)
  call check(p == 0 .and. len(why) > 0, 'B4 none for a burning material: '//why)
  p = birth_family('C', why)
  call check(p == 0 .and. len(why) > 0, 'B4 none for an IG family: '//why)
  p = birth_family('D', why)
  call check(p == 0 .and. len(why) > 0, 'B4 none for two families: '//why)
  p = birth_family('X', why)
  call check(p == 0 .and. len(why) > 0, 'B4 none for an unknown material: '//why)

  ! B5: what a burning material generates
  call set_generation(2, 'B', 1._R8, why)
  call check(len(why) > 0 .and. obj_condensed(2)%genFam == 0, 'B5 refused for an inert material: '//why)
  obj_condensed(1)%dGen = 0._R8
  call set_generation(1, 'B', 1.89_R8, why)
  call check(len(why) > 0 .and. obj_condensed(1)%genFam == 0, 'B5 refused without d-gen: '//why)
  obj_condensed(1)%dGen = 1.e-6_R8
  call set_generation(1, 'B', 1.89_R8, why)
  call check(why == '' .and. obj_condensed(1)%genFam == 2 .and. obj_condensed(1)%nuGen == 1.89_R8, &
             'B5 taken for the burning material')

  ! B6: ICE's own generated mass, time-accurate
  obj_time_scheme%time_accurate = .true.
  grid%blk(1)%cond_phase(1)%dt = dt
  grid%blk(1)%cond_phase(2)%prim(:,2,1,1) = drops()
  c0 = prim_2_cons_MK(grid%blk(1)%cond_phase(2)%prim(:,2,1,1), obj_condensed(2))
  grid%blk(1)%cond_phase(1)%gen_due(:,2,1,1) = pk
  call gen_birth(grid)
  c1 = prim_2_cons_MK(grid%blk(1)%cond_phase(2)%prim(:,2,1,1), obj_condensed(2))
  call check(all(abs(c1 - (c0 + dt*pk)) <= 1.e-12_R8*abs(c0 + dt*pk)), 'B6 time-accurate: born at once, dt times it')
  call check(all(grid%blk(1)%cond_phase(1)%gen_due == 0._R8), 'B6 gen_due cleared')

  ! B7: ICE's own generated mass, steady
  obj_time_scheme%time_accurate = .false.
  grid%blk(1)%cond_phase(1)%gen_due(:,2,1,1) = pk
  call gen_birth(grid)
  call check(all(grid%blk(1)%cond_phase(2)%births_gen(:,2,1,1) == pk) .and. &
             all(grid%blk(1)%cond_phase(2)%births_gen(:,1,1,1) == 0._R8), 'B7 steady: the births in force')
  call check(all(grid%blk(1)%cond_phase(1)%gen_due == 0._R8), 'B7 gen_due cleared')

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_exchange: ', nfail, ' FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_exchange: PASS'

contains

  !> Four materials and five families, one block of two cells, drops in every MK cell
  subroutine setup()
    integer :: f
    ngroups = 5
    allocate(ncond(5), nbase(5), mat_of(5), solid_of(5))
    ncond    = [6, 6, 7, 6, 6]
    nbase    = ncond
    mat_of   = [1, 2, 3, 4, 4]
    solid_of = .false.
    allocate(obj_time_scheme%model(5))
    obj_time_scheme%model = ['MK', 'MK', 'IG', 'MK', 'MK']
    allocate(obj_condensed(4))
    do f = 1, 4
      obj_condensed(f)%name   = achar(iachar('A') + f - 1)
      obj_condensed(f)%rho_al = RHO
      obj_condensed(f)%cs_al  = CP
    enddo
    obj_condensed(1)%combSelect = 1
    obj_condensed(1)%qComb      = QCOMB
    obj_condensed(1)%dGen       = 1.e-6_R8
    grid%nb = 1
    allocate(grid%blk(1))
    grid%blk(1)%dim = [2, 1, 1]
    allocate(grid%blk(1)%gas_phase)
    allocate(grid%blk(1)%gas_phase%prim(5,2,1,1), grid%blk(1)%gas_phase%R(2,1,1), grid%blk(1)%gas_phase%gam(2,1,1), &
             grid%blk(1)%gas_phase%mu(2,1,1), grid%blk(1)%gas_phase%k(2,1,1))
    grid%blk(1)%gas_phase%prim = 0._R8
    allocate(grid%blk(1)%cond_phase(5))
    do f = 1, 5
      associate (cp_ => grid%blk(1)%cond_phase(f))
      allocate(cp_%prim(ncond(f),2,1,1), cp_%source(ncond(f),2,1,1), cp_%residual(ncond(f),2,1,1), cp_%dt(2,1,1), &
               cp_%heat(2,1,1), cp_%removed(6,2,1,1), cp_%gen_rate(5,2,1,1), cp_%gen_due(6,2,1,1))
      cp_%prim = 0._R8; cp_%source = 0._R8; cp_%residual = 0._R8; cp_%dt = 0._R8
      cp_%heat = 0._R8; cp_%removed = 0._R8; cp_%gen_rate = 0._R8; cp_%gen_due = 0._R8
      if (ncond(f) == 6) then
        cp_%prim(:,1,1,1) = drops()
        cp_%prim(:,2,1,1) = drops()
      endif
      end associate
    enddo
  end subroutine setup


  !> An MK state: 0.1 kg/m^3 of drops moving at (5, 1, 0) m/s at 1500 K, 1e12 drops per m^3
  pure function drops() result(q)
    real(R8) :: q(6)
    q = [0.1_R8, 5._R8, 1._R8, 0._R8, 1500._R8, 1.e12_R8]
  end function drops


  logical function near1(a, b)
    real(R8), intent(in) :: a, b
    near1 = abs(a - b) <= 4._R8*epsilon(1._R8)*abs(b)
  end function near1


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

end program test_exchange
