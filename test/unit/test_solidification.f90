!> Unit test of ICE_Lib_Solidification (IGLOO's solid-box material; a hypercooled variant with f0 > 1; a cold solid
!> whose energy is negative): the energy against IGLOO's hSolid on its three branches, the recalescence conserving
!> the energy, the nucleation threshold, melting at T-melt, the injection rule, the plateau rate, the ranges, the
!> nucleated fraction's events and its majority threshold; then (U9) the closure wrappers of ICE_Lib_Solid and (U12-U15)
!> the reconstruction of a solidifying family. Prints every assertion and exits non-zero if any failed.
program test_solidification
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ICE_Lib_Solidification
  use ICE_Lib_Solid
  use ICE_Lib_MK,          only: prim_2_cons_MK, cons_2_prim_MK, flux_make_MK, source_make_MK, check_prim_MK
  use ICE_Lib_IG,          only: prim_2_cons_IG, cons_2_prim_IG, flux_make_IG, source_make_IG, check_prim_IG
  use ICE_Lib_AG,          only: prim_2_cons_AG, cons_2_prim_AG, flux_make_AG, source_make_AG, check_prim_AG
  use ICE_Lib_Riemann,     only: assign_riemann, riemann
  use ICE_Lib_Model,          only: check_prim
  use ICE_Lib_Limiters,       only: assign_limiter, limiter
  use ICE_Lib_Reconstruction, only: state_reconstruction
  use ICE_Config_Types_m,  only: condensed_phase_t
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
  call test_wrappers()
  call test_reconstruction()

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


  !> U9: the wrappers. Round trips on the branches (liquid, a liquid band content with chi = 0.3, a nucleated band content,
  !> the plateau, the solid); a liquid state gives the base's slots bit for bit through the four wrappers; f and chi
  !> move with the mass flux; Saurel's variant with f = 0 is Saurel in flux(1:6); a cold solid with a negative energy
  !> comes back at its temperature.
  subroutine test_wrappers()
    type(condensed_phase_t) :: mat, cold
    real(R8) :: states(3,5), pr(14), cn(14), bk(14), fx(14), sr(14), base(12), force(6), nrm(3), fa(8), fb(8)
    real(R8) :: p1(8), p4(8), wT, wf, wx, wo
    logical  :: bit
    integer  :: k, s, nb, iT, j
    character(len=160) :: msg

    mat%cs_al = cl; mat%rho_al = 2950._R8; mat%Tmelt = Tm; mat%Tnuc = Tn; mat%hFus = hf; mat%cpSol = cs
    mat%solidSelect = 1; mat%solid = .true.
    states(:,1) = [2400._R8, 0._R8, 0._R8]
    states(:,2) = [2000._R8, 0._R8, 0.3_R8]
    states(:,3) = [Tm, 0.3_R8, 0.7_R8]
    states(:,4) = [Tm, 0.7_R8, 1._R8]
    states(:,5) = [1500._R8, 1._R8, 1._R8]
    nrm   = [1._R8, 0._R8, 0._R8]
    force = [0._R8, 0.1_R8, -0.2_R8, 0.05_R8, 3.e3_R8, 7.e-3_R8]

    wT = 0._R8; wf = 0._R8; wx = 0._R8; wo = 0._R8; bit = .true.
    do k = 1, 3
      nb = merge(6, merge(7, 12, k == 2), k == 1)
      iT = nb - 1
      do s = 1, 5
        pr = 0._R8
        pr(1:4) = [2.95_R8, 10._R8, 2._R8, -1._R8]
        if (k == 2) pr(5) = 1.e-6_R8
        if (k == 3) pr(5:10) = [1.e-6_R8, 0._R8, 0._R8, 1.e-6_R8, 0._R8, 1.e-6_R8]
        pr(iT) = states(1,s); pr(nb) = 1.e8_R8; pr(nb+1) = states(2,s); pr(nb+2) = states(3,s)
        select case (k)
        case (1); cn(1:nb+2) = prim_2_cons_MK_S(pr(1:nb+2), mat); bk(1:nb+2) = cons_2_prim_MK_S(cn(1:nb+2), mat)
        case (2); cn(1:nb+2) = prim_2_cons_IG_S(pr(1:nb+2), mat); bk(1:nb+2) = cons_2_prim_IG_S(cn(1:nb+2), mat)
        case (3); cn(1:nb+2) = prim_2_cons_AG_S(pr(1:nb+2), mat); bk(1:nb+2) = cons_2_prim_AG_S(cn(1:nb+2), mat)
        end select
        wT = max(wT, abs(bk(iT) - pr(iT)))
        wf = max(wf, abs(bk(nb+1) - pr(nb+1)))
        wx = max(wx, abs(bk(nb+2) - pr(nb+2)))
        ! every other slot is the base closure's own round trip, bit for bit
        select case (k)
        case (1); base(1:nb) = cons_2_prim_MK(prim_2_cons_MK(pr(1:nb), mat), mat)
        case (2); base(1:nb) = cons_2_prim_IG(prim_2_cons_IG(pr(1:nb), mat), mat)
        case (3); base(1:nb) = cons_2_prim_AG(prim_2_cons_AG(pr(1:nb), mat), mat)
        end select
        do j = 1, nb
          if (j /= iT .and. .not. same([bk(j)], [base(j)])) wo = wo + 1._R8
        enddo
        if (s <= 2) then
          select case (k)
          case (1)
            base(1:nb) = prim_2_cons_MK(pr(1:nb), mat); bit = bit .and. same(cn(1:nb), base(1:nb))
            base(1:nb) = cons_2_prim_MK(cn(1:nb), mat); bit = bit .and. same(bk(1:nb), base(1:nb))
            fx(1:nb+2) = flux_make_MK_S(pr(1:nb+2), nrm, mat); base(1:nb) = flux_make_MK(pr(1:nb), nrm, mat)
            bit = bit .and. same(fx(1:nb), base(1:nb))
            sr(1:nb+2) = source_make_MK_S(pr(1:nb+2), force, mat); base(1:nb) = source_make_MK(pr(1:nb), force, mat)
            bit = bit .and. same(sr(1:nb), base(1:nb))
          case (2)
            base(1:nb) = prim_2_cons_IG(pr(1:nb), mat); bit = bit .and. same(cn(1:nb), base(1:nb))
            base(1:nb) = cons_2_prim_IG(cn(1:nb), mat); bit = bit .and. same(bk(1:nb), base(1:nb))
            fx(1:nb+2) = flux_make_IG_S(pr(1:nb+2), nrm, mat); base(1:nb) = flux_make_IG(pr(1:nb), nrm, mat)
            bit = bit .and. same(fx(1:nb), base(1:nb))
            sr(1:nb+2) = source_make_IG_S(pr(1:nb+2), force, mat); base(1:nb) = source_make_IG(pr(1:nb), force, mat)
            bit = bit .and. same(sr(1:nb), base(1:nb))
          case (3)
            base(1:nb) = prim_2_cons_AG(pr(1:nb), mat); bit = bit .and. same(cn(1:nb), base(1:nb))
            base(1:nb) = cons_2_prim_AG(cn(1:nb), mat); bit = bit .and. same(bk(1:nb), base(1:nb))
            fx(1:nb+2) = flux_make_AG_S(pr(1:nb+2), nrm, mat); base(1:nb) = flux_make_AG(pr(1:nb), nrm, mat)
            bit = bit .and. same(fx(1:nb), base(1:nb))
            sr(1:nb+2) = source_make_AG_S(pr(1:nb+2), force, mat); base(1:nb) = source_make_AG(pr(1:nb), force, mat)
            bit = bit .and. same(sr(1:nb), base(1:nb))
          end select
          bit = bit .and. cn(nb+1) == 0._R8 .and. same([cn(nb+2)], [pr(1)*pr(nb+2)])
        endif
        if (s == 3) then
          select case (k)
          case (1); fx(1:nb+2) = flux_make_MK_S(pr(1:nb+2), nrm, mat)
          case (2); fx(1:nb+2) = flux_make_IG_S(pr(1:nb+2), nrm, mat)
          case (3); fx(1:nb+2) = flux_make_AG_S(pr(1:nb+2), nrm, mat)
          end select
          bit = bit .and. same(fx(nb+1:nb+2), [fx(1)*pr(nb+1), fx(1)*pr(nb+2)])
        endif
      enddo
    enddo
    write(msg,'(A,ES9.2,A,ES9.2,A,ES9.2,A,I0,A)') 'U9 round trips MK/IG/AG x 5 states: T ', wT, ' K, f ', wf, &
                                                  ', chi ', wx, '; ', int(wo), ' other slots off the base round trip'
    call check(wT <= 1.e-9_R8 .and. wf <= 1.e-12_R8 .and. wx <= 4._R8*spacing(1._R8) .and. wo == 0._R8, trim(msg))
    call check(bit, 'U9 liquid states give the base slots bit for bit through the four wrappers; f, chi move with the mass')

    ! Saurel's variant with f = 0 and chi = 0: Saurel in flux(1:6), for veln > 0, < 0 and = 0
    bit = .true.
    do s = 1, 3
      p1 = [2.95_R8, 10._R8, 0._R8, 0._R8, 2400._R8, 1.e8_R8, 0._R8, 0._R8]
      p4 = [2.00_R8,  8._R8, 1._R8, 0._R8, 2200._R8, 7.e7_R8, 0._R8, 0._R8]
      if (s == 2) then; p1(2) = -10._R8; p4(2) = -8._R8; endif
      if (s == 3) then; p1(2) = 5._R8;   p4(2) = -5._R8; endif
      call assign_riemann('Saurel');          fa = riemann(p1, p4, nrm, mat)
      call assign_riemann('Saurel', .true.);  fb = riemann(p1, p4, nrm, mat)
      bit = bit .and. same(fa(1:6), fb(1:6)) .and. fb(7) == 0._R8 .and. fb(8) == 0._R8
    enddo
    call check(bit, 'U9 Saurel with f = chi = 0 is Saurel in flux(1:6) for veln > 0, < 0, = 0; no f, chi flux')

    ! A cold solid with a negative energy (c_l = c_s = 1380, h_fus 1.09e6, 390 K) comes back at 390 K
    cold = mat; cold%cs_al = 1380._R8; cold%cpSol = 1380._R8; cold%hFus = 1.09e6_R8
    pr = 0._R8
    pr(1:8) = [2.95_R8, 10._R8, 0._R8, 0._R8, 390._R8, 1.e8_R8, 1._R8, 1._R8]
    cn(1:8) = prim_2_cons_MK_S(pr(1:8), cold); bk(1:8) = cons_2_prim_MK_S(cn(1:8), cold)
    write(msg,'(A,ES11.4,A,F10.4,A)') 'U9 cold solid, rho e = ', cn(5) - 0.5_R8*pr(1)*100._R8, ' J/m^3: back at ', bk(5), &
                                      ' K (390 K)'
    call check(abs(bk(5) - 390._R8) <= 1.e-9_R8 .and. bk(7) == 1._R8 .and. bk(8) == 1._R8, trim(msg))
  end subroutine test_wrappers


  !> U12-U15: the reconstruction of a solidifying family (van Leer, unit spacing). U12 a liquid stencil gives the faces
  !> of the ordinary reconstruction bit for bit (at temperatures where c_l T / c_l is not T); U13 across a nucleation front
  !> each face holds the energy of its reconstructed T_liq and the phase solid_state gives it; U14 a cold solid (negative
  !> T_liq) gets second-order faces at its physical temperature; U15 no check_prim reads T or f.
  subroutine test_reconstruction()
    type(condensed_phase_t) :: mat, cold
    real(R8), parameter :: T12(4) = [2400.1234567891_R8, 2390.9876543219_R8, 2000.3141592653_R8, 1990.2718281828_R8], &
                           T13(4) = [1875._R8, 1868._R8, Tm, Tm], &
                           f13(4) = [0._R8, 0._R8, 0.55_R8, 0.57_R8]
    real(R8) :: c(14,4), fl(14), fr(14), gl(14), gr(14), tl(4), el, er, xl, xr, T, f, chi
    logical  :: bit, ok
    integer  :: k, nb, iT, j
    character(len=200) :: msg

    mat%cs_al = cl; mat%rho_al = 2950._R8; mat%Tmelt = Tm; mat%Tnuc = Tn; mat%hFus = hf; mat%cpSol = cs
    mat%solidSelect = 1; mat%solid = .true.
    call assign_limiter('vanleer')

    ! U12: liquid states (f = 0; chi 0 and 0.3 in the band, or above T_m), MK, IG and AG
    bit = .true.
    do k = 1, 3
      nb = merge(6, merge(7, 12, k == 2), k == 1)
      iT = nb - 1
      select case (k)
      case (1); check_prim => check_prim_MK
      case (2); check_prim => check_prim_IG
      case (3); check_prim => check_prim_AG
      end select
      c = 0._R8
      do j = 1, 4
        c(1:4,j) = [2.95_R8 - 0.1_R8*j, 10._R8 + j, 0.5_R8*j, -0.2_R8*j]
        if (k == 2) c(5,j) = 1.e-6_R8*j
        if (k == 3) c(5:10,j) = [1.e-6_R8*j, 0._R8, 0._R8, 2.e-6_R8, 0._R8, 1.e-6_R8]
        c(iT,j) = T12(j); c(nb,j) = 1.e8_R8
        c(nb+2,j) = merge(0._R8, 0.3_R8, j <= 2)
      enddo
      call state_reconstruction(c(1:nb+2,1), c(1:nb+2,2), c(1:nb+2,3), c(1:nb+2,4), 1._R8, 1._R8, 1._R8, 0.5_R8, &
                                0.5_R8, fl(1:nb+2), fr(1:nb+2), 1._R8, mat, .true.)
      call state_reconstruction(c(1:nb+2,1), c(1:nb+2,2), c(1:nb+2,3), c(1:nb+2,4), 1._R8, 1._R8, 1._R8, 0.5_R8, &
                                0.5_R8, gl(1:nb+2), gr(1:nb+2), 1._R8, mat, .false.)
      bit = bit .and. same(fl(1:nb+2), gl(1:nb+2)) .and. same(fr(1:nb+2), gr(1:nb+2))
    enddo
    call check(bit, 'U12 a liquid stencil (MK, IG, AG) gives the ordinary faces bit for bit')

    ! U13: liquid 1875 and 1868 K, then plateau cells at T_m (f 0.55, 0.57, chi = 1), MK
    check_prim => check_prim_MK
    c = 0._R8
    do j = 1, 4
      c(1:6,j) = [2.95_R8, 10._R8, 0._R8, 0._R8, T13(j), 1.e8_R8]
      c(7,j) = f13(j); c(8,j) = merge(0._R8, 1._R8, j <= 2)
      tl(j) = c(5,j) + solid_de(c(5,j), c(7,j), cl, cs, hf, Tm)/cl
    enddo
    call state_reconstruction(c(1:8,1), c(1:8,2), c(1:8,3), c(1:8,4), 1._R8, 1._R8, 1._R8, 0.5_R8, 0.5_R8, &
                              fl(1:8), fr(1:8), 1._R8, mat, .true.)
    xl = tl(2) + limiter(tl(3) - tl(2), tl(2) - tl(1))*0.5_R8
    xr = tl(3) - limiter(tl(4) - tl(3), tl(3) - tl(2))*0.5_R8
    el = abs(cl*fl(5) + solid_de(fl(5), fl(7), cl, cs, hf, Tm) - cl*xl)
    er = abs(cl*fr(5) + solid_de(fr(5), fr(7), cl, cs, hf, Tm) - cl*xr)
    ok = el <= etol .and. er <= etol .and. min(tl(2), tl(3)) <= xl .and. xl <= max(tl(2), tl(3)) &
         .and. min(tl(2), tl(3)) <= xr .and. xr <= max(tl(2), tl(3))
    call solid_state(xl, fl(8), cl, cs, Tm, Tn, hf, T, f, chi)
    ok = ok .and. abs(fl(5) - T) <= Ttol .and. abs(fl(7) - f) <= 4._R8*spacing(1._R8)
    call solid_state(xr, fr(8), cl, cs, Tm, Tn, hf, T, f, chi)
    ok = ok .and. abs(fr(5) - T) <= Ttol .and. abs(fr(7) - f) <= 4._R8*spacing(1._R8)
    write(msg,'(A,F9.4,A,F8.6,A,F4.2,A,F9.4,A,F8.6,A,F4.2,A,ES9.2,A)') 'U13 front faces: left T ', fl(5), ' f ', fl(7), &
      ' chi ', fl(8), ', right T ', fr(5), ' f ', fr(7), ' chi ', fr(8), '; energy off ', max(el, er), &
      ' J/kg, each (T, f) = solid_state(T_liq face, chi face)'
    call check(ok .and. fl(7) == 0._R8 .and. fr(7) > 0._R8, trim(msg))

    ! U14: U5's cold solid (c_l = c_s = 1380, h_fus 1.09e6) at 380, 390, 400, 410 K: T_liq from -410 to -380 K
    cold = mat; cold%cs_al = 1380._R8; cold%cpSol = 1380._R8; cold%hFus = 1.09e6_R8
    c = 0._R8
    do j = 1, 4
      c(1:8,j) = [2.95_R8, 10._R8, 0._R8, 0._R8, 370._R8 + 10._R8*j, 1.e8_R8, 1._R8, 1._R8]
    enddo
    call state_reconstruction(c(1:8,1), c(1:8,2), c(1:8,3), c(1:8,4), 1._R8, 1._R8, 1._R8, 0.5_R8, 0.5_R8, &
                              fl(1:8), fr(1:8), 1._R8, cold, .true.)
    write(msg,'(A,F9.4,A,F9.4,A,F9.4,A)') 'U14 cold solid (T_liq ', c(5,2) - 1.09e6_R8/1380._R8, ' K): faces at ', &
      fl(5), ' and ', fr(5), ' K (second order: 395 K; first order would be 390 and 400 K), f = 1'
    call check(abs(fl(5) - 395._R8) <= 1.e-9_R8 .and. abs(fr(5) - 395._R8) <= 1.e-9_R8 .and. fl(7) == 1._R8 &
               .and. fr(7) == 1._R8, trim(msg))

    ! U15: no check_prim reads T or f, so the back-transform's place in the halving loop has no effect today
    ok = .true.
    do k = 1, 3
      nb = merge(6, merge(7, 12, k == 2), k == 1)
      iT = nb - 1
      c = 0._R8
      c(1:4,1) = [2.95_R8, 10._R8, 0._R8, 0._R8]
      if (k == 2) c(5,1) = 1.e-6_R8
      if (k == 3) c(5:10,1) = [1.e-6_R8, 0._R8, 0._R8, 1.e-6_R8, 0._R8, 1.e-6_R8]
      c(iT,1) = 1500._R8; c(nb,1) = 1.e8_R8; c(nb+1,1) = 0.5_R8; c(nb+2,1) = 1._R8
      gl(1:nb+2) = c(1:nb+2,1); gl(iT) = -1.e3_R8; gl(nb+1) = 5._R8
      select case (k)
      case (1); ok = ok .and. check_prim_MK(c(1:nb+2,1)) .and. check_prim_MK(gl(1:nb+2))
      case (2); ok = ok .and. check_prim_IG(c(1:nb+2,1)) .and. check_prim_IG(gl(1:nb+2))
      case (3); ok = ok .and. check_prim_AG(c(1:nb+2,1)) .and. check_prim_AG(gl(1:nb+2))
      end select
    enddo
    call check(ok, 'U15 check_prim (MK, IG, AG) accepts T = -1e3 K and f = 5: it reads neither')
  end subroutine test_reconstruction


  !> Bitwise equality of two arrays of reals.
  logical function same(a, b)
    real(R8), intent(in) :: a(:), b(:)
    same = all(transfer(a, 0_I8, size(a)) == transfer(b, 0_I8, size(b)))
  end function same

end program test_solidification
