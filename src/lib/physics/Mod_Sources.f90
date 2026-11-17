module ICE_Mod_Sources
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
  use ICE_Mod_MPI, only: is_local_block

  implicit none
  private
  public :: compute_source, burnout_switch, smoke_accumulate, smoke_birth, bound_burn_dt
  !> compute_source_ for the unit test; cell_surface for the unit test and the probes
  public :: compute_source_, cell_surface

  !> Families already warned of an unresolved evaporation step under boiling = ZGR
  logical, allocatable :: stiff_warned(:)

contains

  !> The sources of family p. The slots frozen at set-up (setup_freeze) exchange nothing with the gas:
  !> their source is zero, also for a coupled solver, which takes the source of the gas from this array.
  subroutine compute_source (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Mod_ThreadGroups, only: n_tgroups
    use ICE_Mod_Freeze,       only: nfrozen, frozen
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, nf
    logical          :: zgr
    real(kind=R8)    :: stiff

    if (n_tgroups > 1) then
      call compute_source_grouped(grid, p)
      return
    end if

    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)
    zgr   = obj_condensed(mat_of(p))%boilSelect > 0
    stiff = 0._R8

    !$OMP PARALLEL REDUCTION(max:stiff)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
  
      !$OMP DO COLLAPSE(3) schedule(runtime)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        call compute_source_ (grid%blk(b)%cond_phase(p)%prim(:,i,j,k),   &
                              grid%blk(b)%cond_phase(p)%tau(i,j,k),      &
                              grid%blk(b)%gas_phase%prim(:,i,j,k),       &
                              grid%blk(b)%gas_phase%R(i,j,k),            &
                              grid%blk(b)%gas_phase%gam(i,j,k),          &
                              grid%blk(b)%gas_phase%k(i,j,k),            &
                              grid%blk(b)%gas_phase%mu(i,j,k),           &
                              grid%blk(b)%cond_phase(p)%source(:,i,j,k), &
                              obj_condensed(mat_of(p)), nbase(p),        &
                              gas_heat=grid%blk(b)%cond_phase(p)%heat(i,j,k), &
                              burnt=grid%blk(b)%cond_phase(p)%burn(i,j,k))
        if (nf > 0) grid%blk(b)%cond_phase(p)%source(frozen(1:nf,p),i,j,k) = 0._R8
        if (zgr) stiff = max(stiff, step_mass_fraction(grid%blk(b)%cond_phase(p), i, j, k))
  
      enddo ; enddo ; enddo
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL

    if (zgr) call stiff_warning(grid, p, stiff)

  end subroutine compute_source
          
  
  !> compute_source with thread groups: each thread sweeps its slice of its group's blocks.
  subroutine compute_source_grouped(grid, p)
    use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, thread_group, cell_slice
    !$ use omp_lib, only: omp_get_thread_num
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Mod_Freeze,       only: nfrozen, frozen
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, n, tid, g, m, s, i1, j1, k1, i2, j2, k2, nf
    logical          :: zgr
    real(kind=R8)    :: stiff

    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)
    zgr   = obj_condensed(mat_of(p))%boilSelect > 0
    stiff = 0._R8

    !$OMP PARALLEL PRIVATE(b, i, j, k, n, tid, g, m, s, i1, j1, k1, i2, j2, k2) REDUCTION(max:stiff)
    tid = 0
    !$ tid = omp_get_thread_num()
    call thread_group(tid, g, m, s)
    do n = tg_first(g), tg_first(g+1) - 1
      b = tg_blocks(n)
      call cell_slice(grid%blk(b)%dim, m, s, i1, j1, k1, i2, j2, k2)
      do k = k1, k2
      do j = merge(j1, 1, k == k1), merge(j2, grid%blk(b)%dim(2), k == k2)
      do i = merge(i1, 1, k == k1 .and. j == j1), merge(i2, grid%blk(b)%dim(1), k == k2 .and. j == j2)
          call compute_source_ (grid%blk(b)%cond_phase(p)%prim(:,i,j,k),   &
                                grid%blk(b)%cond_phase(p)%tau(i,j,k),      &
                                grid%blk(b)%gas_phase%prim(:,i,j,k),       &
                                grid%blk(b)%gas_phase%R(i,j,k),            &
                                grid%blk(b)%gas_phase%gam(i,j,k),          &
                                grid%blk(b)%gas_phase%k(i,j,k),            &
                                grid%blk(b)%gas_phase%mu(i,j,k),           &
                                grid%blk(b)%cond_phase(p)%source(:,i,j,k), &
                                obj_condensed(mat_of(p)), nbase(p),        &
                                gas_heat=grid%blk(b)%cond_phase(p)%heat(i,j,k), &
                                burnt=grid%blk(b)%cond_phase(p)%burn(i,j,k))
          if (nf > 0) grid%blk(b)%cond_phase(p)%source(frozen(1:nf,p),i,j,k) = 0._R8
          if (zgr) stiff = max(stiff, step_mass_fraction(grid%blk(b)%cond_phase(p), i, j, k))
      enddo ; enddo ; enddo
    enddo
    !$OMP END PARALLEL

    if (zgr) call stiff_warning(grid, p, stiff)

  end subroutine compute_source_grouped


  !> Steady runs, a burning family: no cell steps further than its burn can take its drops' mass. The burn per unit
  !> mass grows without bound as the drops burn down (Beckstead 3 K X/(n d^n), direct as d^-2) and nothing else bounds
  !> the step by it (the tau cap bounds the drag, and only above rho_empty): each cell's step at most burn_dt_factor
  !> times its bulk density over its burn, the burn as compute_source_ evaluates it at the cell's state (the cell's tau
  !> kept), in every cell, whatever its density; and a cell at or below rho_empty that holds drops, the burnt-down tail
  !> whose relaxation time falls with the drops' size, takes compute_dt's tau cap too. A bound on the step moves no
  !> fixed point. After compute_dt.
  subroutine bound_burn_dt (grid, p, tau_factor)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Lib_Model,      only: assign_all, cons_2_prim
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4),      intent(in)    :: p
    real(kind=R8),         intent(in)    :: tau_factor
    integer(kind=I4) :: b, i, j, k
    real(kind=R8)    :: tl, gh, bt, s(ncond(p)), zero(ncond(p)), empty1

    if (.not. obj_condensed(mat_of(p))%combSelect > 0) return
    ! compute_source_ takes the closure's procedures: this family's, before its stages bind them
    call assign_all(p)
    ! ICE's empty state, which holds no drops
    zero   = 0._R8
    s      = cons_2_prim(zero, obj_condensed(mat_of(p)))
    empty1 = s(1)
    !$OMP PARALLEL PRIVATE(b, i, j, k, tl, gh, bt, s) SHARED(empty1)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      !$OMP DO COLLAPSE(3) schedule(runtime)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)
        tl = grid%blk(b)%cond_phase(p)%tau(i,j,k)
        s  = 0._R8
        call compute_source_ (grid%blk(b)%cond_phase(p)%prim(:,i,j,k), tl,     &
                              grid%blk(b)%gas_phase%prim(:,i,j,k),             &
                              grid%blk(b)%gas_phase%R(i,j,k),                  &
                              grid%blk(b)%gas_phase%gam(i,j,k),                &
                              grid%blk(b)%gas_phase%k(i,j,k),                  &
                              grid%blk(b)%gas_phase%mu(i,j,k),                 &
                              s, obj_condensed(mat_of(p)), nbase(p), gas_heat=gh, burnt=bt)
        if (bt > 0._R8) grid%blk(b)%cond_phase(p)%dt(i,j,k) = &
          min(grid%blk(b)%cond_phase(p)%dt(i,j,k), burn_dt_factor*grid%blk(b)%cond_phase(p)%prim(1,i,j,k)/bt)
        ! the drops of the burnt-down tail, at or below rho_empty, relax as fast as their size makes them: compute_dt's
        ! tau cap there too, with the tau of this state, wherever the cell holds drops
        if (grid%blk(b)%cond_phase(p)%prim(1,i,j,k) <= rho_empty .and.                                  &
            grid%blk(b)%cond_phase(p)%prim(1,i,j,k) > empty1 .and.                                     &
            grid%blk(b)%cond_phase(p)%prim(nbase(p),i,j,k) > 0._R8 .and. ieee_is_finite(tl) &
            .and. tl > 0._R8 .and. tau_factor > 0._R8) &
          grid%blk(b)%cond_phase(p)%dt(i,j,k) = min(grid%blk(b)%cond_phase(p)%dt(i,j,k), tau_factor*tl)
      enddo ; enddo ; enddo
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL

  end subroutine bound_burn_dt


  !> The fraction of a cell's bulk density its mass source removes in one step, dt |S_1|/rho_p (0 in an empty cell)
  pure function step_mass_fraction(phase, i, j, k) result(f)
    use ICE_Global_m,          only: rho_empty
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_cond_phase_type), intent(in) :: phase
    integer(kind=I4),          intent(in) :: i, j, k
    real(kind=R8) :: f

    f = 0._R8
    if (phase%prim(1,i,j,k) > rho_empty) f = phase%dt(i,j,k)*abs(phase%source(1,i,j,k))/phase%prim(1,i,j,k)

  end function step_mass_fraction


  !> Once per family: a step that removes more than half of a cell's drop mass leaves the switch at the boiling
  !> temperature unresolved, and the result of the run is not to be trusted
  subroutine stiff_warning(grid, p, stiff)
    use ICE_Global_m,         only: ngroups, mat_of
    use ICE_Advanced_Types_m
    use ICE_Mod_MPI,          only: mpi_rank_
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer(kind=I4),      intent(in) :: p
    real(kind=R8),         intent(in) :: stiff

    if (.not. allocated(stiff_warned)) then
      allocate(stiff_warned(ngroups))
      stiff_warned = .false.
    end if
    if (stiff <= 0.5_R8 .or. stiff_warned(p)) return
    stiff_warned(p) = .true.
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A)') ' [WARNING] boiling = ZGR, family ', p, ' (material ', &
      mat_of(p), ', rank ', mpi_rank_, '): at iteration ', grid%iter, ', time ', grid%time, &
      ', a step removes ', stiff, ' of a cell''s drop mass (more than 0.5): the step does not resolve the '// &
      'evaporation at the boiling temperature and the result is not trusted; lower dt-max'

  end subroutine stiff_warning


  !> The kinematics and the surface of a cell's drops, the only place they are formed: the source takes them, and
  !> so does the probe of the surface temperature. Ts is the mean temperature without P2T and below the boiling
  !> regime (boiling, at the boiling temperature TB).
  subroutine cell_surface(cond_prim, gas_prim, gas_R, gas_gam, gas_k, gas_mu, mat, nb, rho_mat, Rp, Re, Ma, Pr, Nu, Ts, &
                          TB, boiling)
    use ICE_Parameters_m, only: pi, I4
    use ICE_Config_Types_m, only: condensed_phase_t, obj_time_scheme
    use ICE_Lib_Heat
    use ICE_Lib_Properties, only: mat_rho
    use ICE_Lib_Drop_Surface, only: drop_surface
    use ICE_Lib_Evaporation, only: evaporation, blowingFactor, nep, iLv, icpv
    use ICE_Lib_Properties, only: mat_rho, mat_psat, mat_lv, mat_cpv, mat_emiss
    implicit none
    real(kind=R8), dimension(:), intent(in)  :: cond_prim, gas_prim
    real(kind=R8),               intent(in)  :: gas_R, gas_gam, gas_k, gas_mu
    type(condensed_phase_t),     intent(in)  :: mat
    integer(kind=I4),            intent(in)  :: nb
    real(kind=R8),               intent(out) :: rho_mat, Rp, Re, Ma, Pr, Nu, Ts, TB
    logical,                     intent(out) :: boiling

    integer(kind=I4) :: n

    !> Slots of the closure: T is n-1, the number density n
    n = nb

    !> Condensed-material density at the particle temperature
    rho_mat = mat_rho(mat, cond_prim(n-1))

    !> Particles radius 
    Rp = (0.75_R8*cond_prim(1) / (cond_prim(n)*pi*rho_mat))**(1._R8/3._R8)

    !> Reynolds number
    Re = 2._R8*gas_prim(1)*Rp*norm2(gas_prim(2:4)-cond_prim(2:4))/gas_mu
    !Re = max(Re,1.e-15_R8)
  
    !> Mach number
    Ma = norm2(gas_prim(2:4)-cond_prim(2:4)) / (gas_gam*gas_R*gas_prim(5))**0.5_R8

    !> Prandtl number
    Pr = gas_gam/(gas_gam-1._R8)*gas_R*gas_mu/gas_k

    !> Nusselt number of the heat law
    Nu = heat(Re,Pr,Ma,obj_time_scheme%heatSelect)

    !> Surface temperature
    call drop_surface(cond_prim, gas_prim, gas_R, gas_gam, gas_k, gas_mu, mat, n, rho_mat, Rp, Re, Nu, Ts, TB, boiling)

  end subroutine cell_surface


  !> The burn-out of family p, after its last stage. Its drops burn without a residue (d-stop 0): in a cell where they
  !> burn (direct at its burning temperature, Beckstead from T-ign) and are at or below the removal size (d-burnout,
  !> else a drop mass of m_burnout), or where the mass the step took is all that is left, every drop goes at once.
  !> The cell is emptied (every primitive to ICE's empty state but its velocity, which stays) and what it held is added
  !> to the amount the family hands over (gone: mass, momentum, energy), which a coupler delivers to its gas whole,
  !> over its own steps, and then clears.
  subroutine burnout_switch (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Parameters_m,   only: pi
    use ICE_Config_Types_m, only: obj_condensed, COMB_DIRECT
    use ICE_Lib_Model,      only: prim_2_cons, cons_2_prim
    use ICE_Lib_Properties, only: mat_rho
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, n
    real(kind=R8)    :: mb, dt, held(ncond(p)), zero(ncond(p)), empty(ncond(p)), Ton, sm(6)
    logical          :: gone

    n    = nbase(p)
    zero = 0._R8
    associate (mat => obj_condensed(mat_of(p)))
    empty = cons_2_prim(zero, mat)
    ! the temperature from which the drops burn: as the burning branches of compute_source_ test it
    Ton = mat%Tign
    if (mat%combSelect == COMB_DIRECT) Ton = mat%Tboil*(1._R8 - 1.e-10_R8)
    !$OMP PARALLEL PRIVATE(b, i, j, k, mb, dt, held, gone, sm)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      !$OMP DO COLLAPSE(3) schedule(runtime)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)
        associate (q => grid%blk(b)%cond_phase(p)%prim(:,i,j,k), s => grid%blk(b)%cond_phase(p)%source(:,i,j,k))
          if (.not. (q(1) > empty(1) .and. q(n) > 0._R8)) cycle
          ! a drop that does not burn is not removed: it stays an inert drop, whatever its size
          if (q(n-1) < Ton) cycle
          dt = grid%blk(b)%cond_phase(p)%dt(i,j,k)
          mb = m_burnout
          if (mat%dBurn > 0._R8) mb = pi/6._R8*mat_rho(mat, q(n-1))*mat%dBurn**3
          gone = q(1) <= mb*q(n)
          if (s(1) < 0._R8) gone = gone .or. q(1) <= -s(1)*dt
          if (.not. (gone .and. dt > 0._R8)) cycle
          held = prim_2_cons(q, mat)
          ! the remnant's metal burns at once: its smoke is born with its velocity, at its temperature, over this step
          ! (born holds rates)
          sm = 0._R8
          if (mat%smokeFam > 0) then
            sm = smoke_of(held(1), q(2:4), q(n-1), mat)
            grid%blk(b)%cond_phase(p)%born(:,i,j,k) = grid%blk(b)%cond_phase(p)%born(:,i,j,k) + sm/dt
          endif
          grid%blk(b)%cond_phase(p)%gone(1,i,j,k)   = grid%blk(b)%cond_phase(p)%gone(1,i,j,k)   + held(1)
          grid%blk(b)%cond_phase(p)%gone(2:4,i,j,k) = grid%blk(b)%cond_phase(p)%gone(2:4,i,j,k) + held(2:4) - sm(2:4)
          grid%blk(b)%cond_phase(p)%gone(5,i,j,k)   = grid%blk(b)%cond_phase(p)%gone(5,i,j,k)   + held(n-1) - sm(5)
          ! the drops are gone: ICE's empty state but for the velocity, which stays, since the MK flux upwinds at the
          ! mean of two cells' normal velocities and an emptied cell at rest halves its inflow; the temperature too goes
          ! to the empty state, below T-ign, or the empty state's one large drop would burn
          q(1)  = empty(1)
          q(5:) = empty(5:)
        end associate
      enddo ; enddo ; enddo
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL
    end associate

  end subroutine burnout_switch


  !> The smoke of a mass dm [kg/m^3] (or a rate [kg/(m^3 s)]) of a burning material, born with velocity u at
  !> temperature T: mass, momentum (3), energy (the smoke material's in ICE's datum, with the kinetic part), number
  !> (each drop at d-smoke).
  pure function smoke_of(dm, u, T, mat) result(sm)
    use ICE_Parameters_m,   only: pi
    use ICE_Config_Types_m, only: condensed_phase_t, obj_condensed
    use ICE_Global_m,       only: mat_of
    use ICE_Lib_Properties, only: mat_rho, mat_e, mat_cp_const
    real(kind=R8),           intent(in) :: dm, u(3), T
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8) :: sm(6), ms, es
    associate (smat => obj_condensed(mat_of(mat%smokeFam)))
      ms = mat%nuSmoke*dm
      if (smat%cs_varies) then
        es = mat_e(smat, T)
      else
        es = mat_cp_const(smat)*T
      endif
      sm(1)   = ms
      sm(2:4) = ms*u
      sm(5)   = ms*(es + 0.5_R8*dot_product(u, u))
      sm(6)   = ms/(pi/6._R8*mat_rho(smat, T)*mat%dSmoke**3)
    end associate
  end function smoke_of


  !> After the source of family p in stage srk (of the nrk the scheme takes): the smoke its drops make, as rates (birth, the coupler takes them
  !> out of the gas's share) and accumulated with the stage weights of the scheme (born; RK_stage on an accumulator
  !> from 0: forward Euler dt r1, RK2 dt (r1 + r2)/2, RK3 dt (r1/6 + r2/6 + 2 r3/3), the weights the state of p
  !> receives its burning with).
  subroutine smoke_accumulate (grid, p, srk)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Lib_RK,         only: RK_accumulate
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4),      intent(in)    :: p, srk
    integer(kind=I4) :: b, i, j, k, n
    real(kind=R8)    :: r(6)

    if (obj_condensed(mat_of(p))%smokeFam == 0) return
    n = nbase(p)
    associate (mat => obj_condensed(mat_of(p)))
    !$OMP PARALLEL PRIVATE(b, i, j, k, r)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      !$OMP DO COLLAPSE(3) schedule(runtime)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)
        r = 0._R8
        if (grid%blk(b)%cond_phase(p)%burn(i,j,k) > 0._R8) &
          r = smoke_of(grid%blk(b)%cond_phase(p)%burn(i,j,k), grid%blk(b)%cond_phase(p)%prim(2:4,i,j,k), &
                       grid%blk(b)%cond_phase(p)%prim(n-1,i,j,k), mat)
        ! the step's birth rates, the stages' weighted as a two-way coupled family's source is (RK_accumulate)
        grid%blk(b)%cond_phase(p)%born(:,i,j,k)  = RK_accumulate(srk, nrk, grid%blk(b)%cond_phase(p)%born(:,i,j,k), r)
        grid%blk(b)%cond_phase(p)%birth(:,i,j,k) = grid%blk(b)%cond_phase(p)%born(1:5,i,j,k)
      enddo ; enddo ; enddo
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL
    end associate

  end subroutine smoke_accumulate


  !> After every family's step: the smoke the burning families made (born, the step's weighted birth rates), and the
  !> accumulators cleared. A time-accurate run (one dt for every family): born in its family, dt times the rates of
  !> mass, momentum, energy and number added to the conservative state cell by cell, with the smoke family's own
  !> closure. A steady run (each family its own local dt): the rates themselves, the source of the smoke family in every
  !> stage of its next step (srate, add_external_source), so that its fixed point holds whatever either family's dt.
  subroutine smoke_birth (grid)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_condensed, obj_time_scheme
    use ICE_Lib_Model,      only: prim_2_cons, cons_2_prim, assign_prim_2_cons, assign_cons_2_prim
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, i, j, k, p, ps, ns
    real(kind=R8)    :: dtp
    real(kind=R8), allocatable :: c(:)
    logical :: steady

    steady = .not. obj_time_scheme%time_accurate
    ! a steady run: this step's rates replace the last step's
    if (steady) then
      do p = 1, ngroups
        ps = obj_condensed(mat_of(p))%smokeFam
        if (ps == 0) cycle
        do b = 1, grid%nb
          if (.not. is_local_block(b)) cycle
          if (.not. allocated(grid%blk(b)%cond_phase(ps)%srate)) &
            allocate(grid%blk(b)%cond_phase(ps)%srate(ncond(ps), grid%blk(b)%dim(1), grid%blk(b)%dim(2), &
                                                      grid%blk(b)%dim(3)))
          grid%blk(b)%cond_phase(ps)%srate = 0._R8
        enddo
      enddo
    endif
    do p = 1, ngroups
      ps = obj_condensed(mat_of(p))%smokeFam
      if (ps == 0) cycle
      ns = nbase(ps)
      ! the conversions of the smoke family's closure: the family loop left those of the last family it stepped
      call assign_prim_2_cons(ps)
      call assign_cons_2_prim(ps)
      associate (smat => obj_condensed(mat_of(ps)))
      !$OMP PARALLEL PRIVATE(b, i, j, k, c, dtp)
      allocate(c(ncond(ps)))
      do b = 1, grid%nb
        if (.not. is_local_block(b)) cycle
        !$OMP DO COLLAPSE(3) schedule(runtime)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)
          if (all(grid%blk(b)%cond_phase(p)%born(:,i,j,k) == 0._R8)) cycle
          if (steady) then
            grid%blk(b)%cond_phase(ps)%srate(1:4,i,j,k) = grid%blk(b)%cond_phase(ps)%srate(1:4,i,j,k) &
                                                        + grid%blk(b)%cond_phase(p)%born(1:4,i,j,k)
            grid%blk(b)%cond_phase(ps)%srate(ns-1,i,j,k) = grid%blk(b)%cond_phase(ps)%srate(ns-1,i,j,k) &
                                                         + grid%blk(b)%cond_phase(p)%born(5,i,j,k)
            grid%blk(b)%cond_phase(ps)%srate(ns,i,j,k) = grid%blk(b)%cond_phase(ps)%srate(ns,i,j,k) &
                                                       + grid%blk(b)%cond_phase(p)%born(6,i,j,k)
          else
            dtp = grid%blk(b)%cond_phase(p)%dt(i,j,k)
            c = prim_2_cons(grid%blk(b)%cond_phase(ps)%prim(:,i,j,k), smat)
            c(1:4)  = c(1:4)  + dtp*grid%blk(b)%cond_phase(p)%born(1:4,i,j,k)
            c(ns-1) = c(ns-1) + dtp*grid%blk(b)%cond_phase(p)%born(5,i,j,k)
            c(ns)   = c(ns)   + dtp*grid%blk(b)%cond_phase(p)%born(6,i,j,k)
            grid%blk(b)%cond_phase(ps)%prim(:,i,j,k) = cons_2_prim(c, smat)
          endif
          grid%blk(b)%cond_phase(p)%born(:,i,j,k)  = 0._R8
        enddo ; enddo ; enddo
        !$OMP END DO nowait
      enddo
      deallocate(c)
      !$OMP END PARALLEL
      end associate
    enddo

  end subroutine smoke_birth


  !> The source of one cell. gas_heat is the heat the gas receives that the source does not carry [W/m^3]: the heat of
  !> reaction of a burning material (combustion = direct), 0 otherwise; the gas takes -source plus gas_heat.
  subroutine compute_source_ (cond_prim, cond_tau, gas_prim, gas_R, gas_gam, gas_k, gas_mu, source, mat, nb, gas_heat, &
                              Ts_out, burnt)
    use ICE_Parameters_m, only: pi, sigma_SB, I4
    use ICE_Config_Types_m, only: condensed_phase_t, obj_time_scheme, COMB_DIRECT, COMB_BECKSTEAD
    use ICE_Lib_Model
    use ICE_Lib_Drag
    use ICE_Lib_Evaporation, only: evaporation, blowingFactor, zgrBoiling, burnRate, becksteadRate, nep, iLv, icpv
    use ICE_Lib_Properties, only: mat_psat, mat_lv, mat_cpv, mat_emiss, mat_e, mat_cp_const, mat_kl
    implicit none
    real(kind=R8), dimension(:), intent(in)    :: cond_prim, gas_prim
    real(kind=R8),               intent(inout) :: cond_tau
    real(kind=R8),               intent(in)    :: gas_R, gas_gam, gas_k, gas_mu
    real(kind=R8), dimension(:), intent(inout) :: source
    type(condensed_phase_t),     intent(in)    :: mat
    integer(kind=I4),            intent(in)    :: nb
    real(kind=R8),               intent(out)   :: gas_heat
    real(kind=R8), optional,     intent(out)   :: Ts_out
    real(kind=R8), optional,     intent(out)   :: burnt   !> the mass the drops burn per unit volume and time [kg/(m^3 s)]

    integer(kind=I4) :: n, ng
    real(kind=R8)    :: Rp, Re, Ma, Pr, Tr, Ts, TB
    real(kind=R8)    :: Cd, Nu
    real(kind=R8)    :: rho_mat, mdot, Qconv, Qevap, carry, cpg, Lv, Qburn, burnt_
    logical          :: override_Qdot, boiling
    real(kind=R8)    :: force(6)
    real(kind=R8)    :: ep(nep)

    !> Slots of the closure: T is n-1, the number density n
    n = nb
    ng = size(gas_prim)

    !> Density, radius, Re, Ma, Pr, Nu, the surface temperature Ts and the boiling regime
    call cell_surface(cond_prim, gas_prim, gas_R, gas_gam, gas_k, gas_mu, mat, n, rho_mat, Rp, Re, Ma, Pr, Nu, Ts, &
                      TB, boiling)

    !> Temperature ratio
    Tr = Ts/gas_prim(5)

    !> Drag coefficient
    Cd = drag(Re,Ma,gas_gam,Tr,obj_time_scheme%dragSelect)

    !> Slip velocity
    force(2:4) = gas_prim(2:4) - cond_prim(2:4)

    !> Convective heat exchange of a single particle [W]
    Qconv = 2._R8*Nu*gas_k*pi*Rp*(gas_prim(ng)-Ts)

    !> Mass exchange of a single particle [kg/s], negative while the particle
    !  evaporates, from the surface state. Qevap comes back in W; the latent sink is not
    !  included here, every closure's source_make adds it as -force(1)*L_v(T). A Psat
    !  column in the property table replaces the Clausius-Clapeyron curve, and Lv and Cpv
    !  columns replace the INI values at the surface temperature.
    ep = mat%ep
    if (mat%use_lv)  ep(iLv)  = mat_lv(mat, Ts)
    if (mat%use_cpv) ep(icpv) = mat_cpv(mat, Ts)
    if (boiling) then
      !> ZGR: the flash fed by the superheat (by the parabolic conduction of a P2T drop) and the heat-fed rate
      cpg = gas_gam*gas_R/(gas_gam - 1._R8)
      if (mat%liqSelect == 1) then
        call zgrBoiling(cond_prim(n-1), TB, 2._R8*Rp, gas_k, cpg, Nu, gas_prim(ng), ep(iLv), mdot, Qevap, &
                        kCond=mat_kl(mat, cond_prim(n-1)))
      else
        call zgrBoiling(cond_prim(n-1), TB, 2._R8*Rp, gas_k, cpg, Nu, gas_prim(ng), ep(iLv), mdot, Qevap)
      endif
      override_Qdot = .true.
    elseif (mat%use_psat) then
      call evaporation(gas_prim(1), gas_prim(ng), gas_gam, gas_R, gas_mu, gas_k,  &
                       Ts, 2._R8*Rp, Re,                                          &
                       mat%evapSelect, mat%intfSelect,                            &
                       ep, mdot, Qevap, override_Qdot,                  &
                       psatExt=mat_psat(mat, Ts))
    else
      call evaporation(gas_prim(1), gas_prim(ng), gas_gam, gas_R, gas_mu, gas_k,  &
                       Ts, 2._R8*Rp, Re,                                          &
                       mat%evapSelect, mat%intfSelect,                            &
                       ep, mdot, Qevap, override_Qdot)
    endif

    if (override_Qdot) then
      !> ASM and TC resolve the gas-side heat inside their own film, Stefan flow
      !  included, so their value replaces the Nusselt one rather than correcting it
      Qconv = Qevap
    elseif (obj_time_scheme%blowSelect == 1) then
      !> Outgoing vapour thickens the thermal film and cuts the convective heat
      Qconv = Qconv * blowingFactor(gas_gam, gas_R, gas_mu, gas_k, rho_mat, &
                                    2._R8*Rp, 4._R8/3._R8*pi*Rp**3*rho_mat, mdot)
    endif

    !> force(1) is the mass leaving the condensed phase per unit volume and time,
    !  force(5) the heat entering it
    force(1) = - mdot*cond_prim(n)
    force(5) = Qconv*cond_prim(n)                  
      
    !> Radiative heat exchange over the whole particle surface
    force(5) = force(5) + mat_emiss(mat, Ts)*sigma_SB * 4._R8*pi*Rp*Rp*cond_prim(n)*(gas_prim(ng)**4._I4-Ts**4._I4)

    !> P2T or boiling: the vapour leaves at the surface temperature; the closure's latent sink at the mean temperature
    !  is corrected to Lv(Ts) + e(Ts) - e(T)
    if (mat%liqSelect == 1 .or. boiling) then
      if (mat%cs_varies) then
        carry = mat_e(mat, Ts) - mat_e(mat, cond_prim(n-1))
      else
        carry = mat_cp_const(mat)*(Ts - cond_prim(n-1))
      endif
      force(5) = force(5) - force(1)*((ep(iLv) - mat_lv(mat, cond_prim(n-1))) + carry)
    endif

    !> Direct combustion: a drop at or above its burning temperature and larger than the residue burns while its rate is
    !  positive. Every heat it receives feeds the vapour, so the closure's latent sink cancels it and the drop keeps
    !  its temperature; the vapour leaves with the drop's energy and its heat of reaction goes to the gas. Otherwise
    !  the drop exchanges as an inert one. The temperature is compared within 1e-10 of boiling-temperature, so that
    !  the round-off of the energy does not switch a drop held there off.
    Qburn = 0._R8; burnt_ = 0._R8
    if (mat%combSelect == COMB_DIRECT) then
      if (cond_prim(n-1) >= mat%Tboil*(1._R8 - 1.e-10_R8) .and. 2._R8*Rp > mat%dRes) then
        cpg  = gas_gam*gas_R/(gas_gam - 1._R8)
        !> the model's latent heat is the material's at the burning temperature; the closure's sink, the material's
        !  at the drop's temperature, is cancelled as it is
        Lv   = mat_lv(mat, mat%Tboil)
        mdot = burnRate(gas_prim(ng), gas_k, cpg, cond_prim(n-1), 2._R8*Rp, Nu, Lv, mat%qComb)
        if (mdot > 0._R8) then
          force(1) = mdot*cond_prim(n)
          force(5) = force(1)*mat_lv(mat, cond_prim(n-1))
          Qburn    = force(1)*mat%qComb
          burnt_   = force(1)
        endif
      endif
    elseif (mat%combSelect == COMB_BECKSTEAD) then
      !> Beckstead: a drop at or above T-ign and larger than the residue burns on the d^n law. The law has no latent
      !  sink, so the closure's is cancelled; beta-part of the heat of reaction heats the drop, and the gas, which
      !  receives the whole of it, pays that share through the drop's energy.
      if (cond_prim(n-1) >= mat%Tign .and. 2._R8*Rp > mat%dRes) then
        mdot     = becksteadRate(2._R8*Rp, rho_mat, mat%Kburn, mat%nBurn, mat%Xeff)
        force(1) = mdot*cond_prim(n)
        force(5) = force(5) + force(1)*(mat_lv(mat, cond_prim(n-1)) + mat%betaPart*mat%qComb)
        Qburn    = force(1)*mat%qComb
        burnt_   = force(1)
      endif
    endif
    gas_heat = Qburn
    if (present(burnt)) burnt = burnt_

    !> Particles relaxation time (infinite without drag, so every relaxation term vanishes exactly)
    if (Cd == 0._R8) then
      cond_tau = ieee_value(1._R8, ieee_positive_inf)
    else
      cond_tau = 8._R8*rho_mat*Rp / (3._R8*gas_prim(1)*Cd*norm2(gas_prim(2:4)-cond_prim(2:4)) + 1e-20)
    endif
    !> Use only for Vie validation test
    !> cond_tau = 5.0_R8
    force(6) = cond_tau

    !> Source terms
    source = source_make(cond_prim,force,mat)

    if (present(Ts_out)) Ts_out = Ts

  end subroutine compute_source_
    

end module ICE_Mod_Sources