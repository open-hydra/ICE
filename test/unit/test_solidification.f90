!> Unit test of ICE_Lib_Solidification (IGLOO's solid-box material; a hypercooled variant with f0 > 1; a cold solid
!> whose energy is negative): the energy against IGLOO's hSolid on its three branches, the recalescence conserving
!> the energy, the nucleation threshold, melting at T-melt, the injection rule, the plateau rate, the ranges, the
!> nucleated fraction's events and its majority threshold. Prints every assertion and exits non-zero if any failed.
program test_solidification
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ICE_Lib_Solidification
  implicit none
  real(R8), parameter :: cl = 1250._R8, cs = 600._R8, hf = 1.07e6_R8, Tm = 2327._R8, hf2 = 4.e5_R8
  real(R8) :: Tn, etol, Ttol
  integer  :: nfail = 0

  Tn   = 0.8_R8*Tm
  etol = 4._R8*spacing(cl*Tm)
  Ttol = etol/min(cl, cs)

  call test_energy()
  call test_jump()
  call test_threshold()
  call test_melting()
  call test_cold_solid()
  call test_injection()
  call test_plateau_rate()
  call test_range()
  call test_memory()
  call test_majority()

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_solidification: ', nfail, ' FAILED'
    error stop 1
  endif
  write(*,'(A)') 'test_solidification: PASS'

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


  !> U1: c_l T - f L(T) is hSolid - hOff on the liquid (f = 0), plateau (T = T_m) and solid (f = 1) branches.
  subroutine test_energy()
    real(R8) :: T, f, w(3)
    character(len=160) :: msg
    integer :: i

    w = 0._R8
    do i = 0, 49
      T = 300._R8 + 2700._R8*real(i, R8)/49._R8
      w(1) = max(w(1), abs(cl*T + solid_de(T, 0._R8, cl, cs, hf, Tm) - hSolid(T, 0._R8, phLiquid, cl, cs, hf, Tm, 0._R8)))
      w(1) = max(w(1), abs(cl*T + solid_de(T, 0._R8, cl, cs, hf, Tm) - hSolid(T, 0._R8, phUndercooled, cl, cs, hf, Tm, 0._R8)))
      w(3) = max(w(3), abs(cl*T + solid_de(T, 1._R8, cl, cs, hf, Tm) - hSolid(T, 1._R8, phSolid, cl, cs, hf, Tm, 0._R8)))
    enddo
    do i = 0, 20
      f = real(i, R8)/20._R8
      w(2) = max(w(2), abs(cl*Tm + solid_de(Tm, f, cl, cs, hf, Tm) - hSolid(Tm, f, phPlateau, cl, cs, hf, Tm, 0._R8)))
    enddo
    write(msg,'(A,3ES10.2,A,ES9.2,A)') 'U1 energy = hSolid - hOff on the liquid/plateau/solid branches (worst', w, &
                                      ' J/kg; tol ', etol, ')'
    call check(maxval(w) <= etol, trim(msg))
  end subroutine test_energy


  !> U2: the state after nucleation from T_ev has the energy c_l T_ev (both branches); at f0 = 1 both give T_m.
  subroutine test_jump()
    real(R8) :: Tev, T, f, chi, h, w(2), Tev1, wc
    character(len=160) :: msg
    integer :: i, m

    w = 0._R8
    do m = 1, 2
      h = merge(hf, hf2, m == 1)
      do i = 0, 39
        Tev = Tn - 600._R8*real(i, R8)/39._R8
        call solid_state(Tev, 0._R8, cl, cs, Tm, Tn, h, T, f, chi)
        w(m) = max(w(m), abs(cl*T + solid_de(T, f, cl, cs, h, Tm) - cl*Tev))
      enddo
    enddo
    write(msg,'(A,2ES10.2,A,ES9.2,A)') 'U2 recalescence conserves the energy, plateau and hypercooled materials (worst', w, &
                                      ' J/kg; tol ', etol, ')'
    call check(maxval(w) <= etol, trim(msg))

    call solid_state(Tn, 0._R8, cl, cs, Tm, Tn, hf2, T, f, chi)
    write(msg,'(A,F10.3,A)') 'U2 hypercooled jump from T_n lands on the solid at ', T, ' K (2024.083 K)'
    call check(abs(T - (Tm - (cl*(Tm - Tn) - hf2)/cs)) <= Ttol .and. abs(T - 2024.083_R8) < 1.e-3_R8 .and. f == 1._R8, &
               trim(msg))

    Tev1 = Tm - hf/cl
    wc = 0._R8
    do i = -1, 1
      Tev = Tev1
      if (i /= 0) Tev = nearest(Tev1, real(i, R8))
      call solid_state(Tev, 0._R8, cl, cs, Tm, Tn, hf, T, f, chi)
      wc = max(wc, abs(T - Tm))
    enddo
    write(msg,'(A,F8.3,A,ES9.2,A)') 'U2 the two branches meet at f0 = 1 (T_ev = ', Tev1, ' K; worst |T - T_m| ', wc, ' K)'
    call check(wc <= Ttol, trim(msg))
  end subroutine test_jump


  !> U3: a liquid content at exactly T_nuc nucleates (IGLOO's injection rule T <= T_nuc); the next double above stays
  !> liquid.
  subroutine test_threshold()
    real(R8) :: T, f, chi, Tup
    character(len=160) :: msg

    call solid_state(Tn, 0._R8, cl, cs, Tm, Tn, hf, T, f, chi)
    write(msg,'(A,F9.6,A,F4.1)') 'U3 at T_liq = T_nuc: plateau at T_m, f = ', f, ' (f0 = 0.543692), chi = ', chi
    call check(T == Tm .and. abs(f - cl*(Tm - Tn)/hf) <= 4._R8*spacing(1._R8) .and. abs(f - 0.543692_R8) < 1.e-6_R8 &
               .and. chi == 1._R8, trim(msg))
    Tup = nearest(Tn, 1._R8)
    call solid_state(Tup, 0._R8, cl, cs, Tm, Tn, hf, T, f, chi)
    call check(T == Tup .and. f == 0._R8 .and. chi == 0._R8, 'U3 one double above T_nuc: liquid, T unchanged, f = chi = 0')
  end subroutine test_threshold


  !> U4: no content is solid above T_m. A solid heated to the energy of T_m + 50 is the plateau at T_m, f = 0.971963;
  !> a nucleated band content is the plateau; at T_liq = T_m it is liquid with chi cleared; over a sweep of nucleated
  !> contents T <= T_m, with T = T_m or f = 1.
  subroutine test_melting()
    real(R8) :: Tl, T, f, chi
    character(len=160) :: msg
    logical  :: ok
    integer  :: i

    Tl = (cl*Tm - hf + 50._R8*cs)/cl
    call solid_state(Tl, 1._R8, cl, cs, Tm, Tn, hf, T, f, chi)
    write(msg,'(A,F8.3,A,F10.4,A,F9.6,A)') 'U4 solid heated to T_m + 50 (T_liq = ', Tl, ' K): T = ', T, ' K, f = ', f, &
                                          ' (2327 K, 0.971963)'
    call check(T == Tm .and. abs(f - (1._R8 - 50._R8*cs/hf)) <= 4._R8*spacing(1._R8) .and. chi == 1._R8 .and. &
               abs(cl*T + solid_de(T, f, cl, cs, hf, Tm) - cl*Tl) <= etol, trim(msg))

    call solid_state(2000._R8, 1._R8, cl, cs, Tm, Tn, hf, T, f, chi)
    write(msg,'(A,F9.6,A)') 'U4 nucleated band content at T_liq = 2000 K: plateau at T_m, f = ', f, ' (0.382009)'
    call check(T == Tm .and. abs(f - 0.382009_R8) < 1.e-6_R8 .and. chi == 1._R8, trim(msg))

    call solid_state(Tm, 1._R8, cl, cs, Tm, Tn, hf, T, f, chi)
    call check(T == Tm .and. f == 0._R8 .and. chi == 0._R8, 'U4 at T_liq = T_m a nucleated content is liquid, chi cleared')

    ok = .true.
    do i = 0, 199
      Tl = (Tn - 600._R8) + (Tm - (Tn - 600._R8))*real(i, R8)/200._R8
      call solid_state(Tl, 1._R8, cl, cs, Tm, Tn, hf, T, f, chi)
      ok = ok .and. T <= Tm .and. (T == Tm .or. f == 1._R8) .and. f > 0._R8
    enddo
    call check(ok, 'U4 200 nucleated contents with T_liq in [T_n - 600, T_m): T <= T_m, T = T_m or f = 1, f > 0')
  end subroutine test_melting


  !> U5: a cold solid whose energy is negative keeps its temperature.
  subroutine test_cold_solid()
    real(R8), parameter :: c2 = 1380._R8, h3 = 1.09e6_R8
    real(R8) :: Tl, T, f, chi, e, t2tol
    character(len=160) :: msg

    e  = c2*Tm - h3 - c2*(Tm - 390._R8)
    Tl = e/c2
    t2tol = 4._R8*spacing(c2*Tm)/c2
    call solid_state(Tl, 1._R8, c2, c2, Tm, 0.8_R8*Tm, h3, T, f, chi)
    write(msg,'(A,ES11.4,A,F10.4,A,F10.4,A)') 'U5 cold solid, e = ', e, ' J/kg (T_liq = ', Tl, ' K): T = ', T, ' K (390 K)'
    call check(abs(T - 390._R8) <= t2tol .and. f == 1._R8 .and. chi == 1._R8, trim(msg))
  end subroutine test_cold_solid


  !> U6: the injection rule.
  subroutine test_injection()
    integer  :: ph(3)
    real(R8) :: f(3)

    call solidPhaseAtInjection(Tn, Tm, Tn, ph(1), f(1))
    call solidPhaseAtInjection(Tn + 1._R8, Tm, Tn, ph(2), f(2))
    call solidPhaseAtInjection(Tm, Tm, Tn, ph(3), f(3))
    call check(ph(1) == phSolid .and. f(1) == 1._R8 .and. ph(2) == phUndercooled .and. f(2) == 0._R8 .and. &
               ph(3) == phLiquid .and. f(3) == 0._R8, 'U6 injection: solid at T_nuc, undercooled above it, liquid at T_m')
  end subroutine test_injection


  !> U7: on the plateau a loss de of energy freezes de/h_fus (IGLOO's plateauRate in energy form).
  subroutine test_plateau_rate()
    real(R8), parameter :: de = 1.e3_R8
    real(R8) :: T1, f1, T2, f2, chi, r
    character(len=160) :: msg

    call solid_state(1700._R8, 0._R8, cl, cs, Tm, Tn, hf, T1, f1, chi)
    call solid_state((cl*1700._R8 - de)/cl, 0._R8, cl, cs, Tm, Tn, hf, T2, f2, chi)
    r = (f2 - f1)/(de/hf)
    write(msg,'(A,F12.9,A)') 'U7 plateau: df/(de/h_fus) = ', r, ' (1)'
    call check(T1 == Tm .and. T2 == Tm .and. abs(r - 1._R8) <= 1.e-9_R8, trim(msg))
  end subroutine test_plateau_rate


  !> U8: ranges over 1000 (T_liq, chi) pairs per material, a Weyl sequence over [T_n - 600, T_m + 100] x [0, 1].
  subroutine test_range()
    real(R8), parameter :: g1 = 0.6180339887498949_R8, g2 = 0.4142135623730951_R8
    real(R8) :: Tl, x, T, f, chi, f0, h
    logical  :: ok
    integer  :: i, m

    ok = .true.
    do m = 1, 2
      h  = merge(hf, hf2, m == 1)
      f0 = min(cl*(Tm - Tn)/h, 1._R8)
      do i = 1, 1000
        Tl = (Tn - 600._R8) + (Tm + 100._R8 - (Tn - 600._R8))*modulo(real(i, R8)*g1, 1._R8)
        x  = modulo(real(i, R8)*g2, 1._R8)
        call solid_state(Tl, x, cl, cs, Tm, Tn, h, T, f, chi)
        ok = ok .and. ieee_is_finite(T) .and. ieee_is_finite(f) .and. chi >= 0._R8 .and. chi <= 1._R8
        if (f > 0._R8) ok = ok .and. T <= Tm .and. (T == Tm .or. f == 1._R8) .and. f <= 1._R8
        if (Tl >= Tm)  ok = ok .and. f == 0._R8 .and. chi == 0._R8
        if (Tl <= Tn)  ok = ok .and. chi == 1._R8 .and. f >= f0*(1._R8 - 4._R8*epsilon(1._R8))
        if (Tl > Tn .and. Tl < Tm) ok = ok .and. chi == x
      enddo
    enddo
    call check(ok, 'U8 ranges over 2000 (T_liq, chi) states: finite, chi in [0,1], f > 0 only with T <= T_m, '// &
               'cleared from T_m, set to T_nuc, kept between')
  end subroutine test_range


  !> U10: the nucleated fraction changes only at nucleation (set) and at complete melting (cleared).
  subroutine test_memory()
    real(R8) :: T, f, chi
    logical  :: ok(4)
    character(len=160) :: msg

    call solid_state(Tn - 1._R8, 0._R8, cl, cs, Tm, Tn, hf, T, f, chi)
    ok(1) = T == Tm .and. f > 0._R8 .and. chi == 1._R8
    call solid_state(Tm, 1._R8, cl, cs, Tm, Tn, hf, T, f, chi)
    ok(2) = T == Tm .and. f == 0._R8 .and. chi == 0._R8
    call solid_state(2000._R8, 0.3_R8, cl, cs, Tm, Tn, hf, T, f, chi)
    ok(3) = T == 2000._R8 .and. f == 0._R8 .and. chi == 0.3_R8
    call solid_state(2000._R8, 0.7_R8, cl, cs, Tm, Tn, hf, T, f, chi)
    ok(4) = T == Tm .and. abs(f - 0.382009_R8) < 1.e-6_R8 .and. chi == 0.7_R8
    write(msg,'(A,4L2,A,F10.4,A,F9.6)') 'U10 chi set at T_nuc, cleared at T_m, kept in the band (', ok, &
                                       '); last row T = ', T, ' K, chi = ', chi
    call check(all(ok), trim(msg))
  end subroutine test_memory


  !> U11: in the band a content is nucleated iff chi >= 1/2.
  subroutine test_majority()
    real(R8) :: T1, f1, c1, T2, f2, c2
    character(len=160) :: msg

    call solid_state(2000._R8, 0.5_R8, cl, cs, Tm, Tn, hf, T1, f1, c1)
    call solid_state(2000._R8, nearest(0.5_R8, -1._R8), cl, cs, Tm, Tn, hf, T2, f2, c2)
    write(msg,'(A,F10.4,A,F10.4,A)') 'U11 majority at T_liq = 2000 K: chi = 1/2 gives T = ', T1, &
                                    ' K, one double below gives ', T2, ' K'
    call check(T1 == Tm .and. f1 > 0._R8 .and. T2 == 2000._R8 .and. f2 == 0._R8, trim(msg))
  end subroutine test_majority

end program test_solidification
