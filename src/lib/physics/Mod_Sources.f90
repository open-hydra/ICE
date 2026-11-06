module ICE_Mod_Sources
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
  use ICE_Mod_MPI, only: is_local_block

  implicit none
  private
  public :: compute_source
  !> compute_source_ for the unit test; cell_surface for the unit test and the probes
  public :: compute_source_, cell_surface

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

    if (n_tgroups > 1) then
      call compute_source_grouped(grid, p)
      return
    end if

    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)

    !$OMP PARALLEL
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
                              obj_condensed(mat_of(p)), nbase(p)         )
        if (nf > 0) grid%blk(b)%cond_phase(p)%source(frozen(1:nf,p),i,j,k) = 0._R8
  
      enddo ; enddo ; enddo
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL

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

    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)

    !$OMP PARALLEL PRIVATE(b, i, j, k, n, tid, g, m, s, i1, j1, k1, i2, j2, k2)
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
                                obj_condensed(mat_of(p)), nbase(p)         )
          if (nf > 0) grid%blk(b)%cond_phase(p)%source(frozen(1:nf,p),i,j,k) = 0._R8
      enddo ; enddo ; enddo
    enddo
    !$OMP END PARALLEL

  end subroutine compute_source_grouped


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


  subroutine compute_source_ (cond_prim, cond_tau, gas_prim, gas_R, gas_gam, gas_k, gas_mu, source, mat, nb, Ts_out)
    use ICE_Parameters_m, only: pi, sigma_SB, I4
    use ICE_Config_Types_m, only: condensed_phase_t, obj_time_scheme
    use ICE_Lib_Model
    use ICE_Lib_Drag
    use ICE_Lib_Evaporation, only: evaporation, blowingFactor, zgrBoiling, nep, iLv, icpv
    use ICE_Lib_Properties, only: mat_psat, mat_lv, mat_cpv, mat_emiss, mat_e, mat_cp_const, mat_kl
    implicit none
    real(kind=R8), dimension(:), intent(in)    :: cond_prim, gas_prim
    real(kind=R8),               intent(inout) :: cond_tau
    real(kind=R8),               intent(in)    :: gas_R, gas_gam, gas_k, gas_mu
    real(kind=R8), dimension(:), intent(inout) :: source
    type(condensed_phase_t),     intent(in)    :: mat
    integer(kind=I4),            intent(in)    :: nb
    real(kind=R8), optional,     intent(out)   :: Ts_out

    integer(kind=I4) :: n, ng
    real(kind=R8)    :: Rp, Re, Ma, Pr, Tr, Ts, TB
    real(kind=R8)    :: Cd, Nu
    real(kind=R8)    :: rho_mat, mdot, Qconv, Qevap, carry, cpg
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