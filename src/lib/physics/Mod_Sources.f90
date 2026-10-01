module ICE_Mod_Sources
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
  use ICE_Mod_MPI, only: is_local_block

  implicit none
  private
  public :: compute_source

contains

  subroutine compute_source (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k

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
  
      enddo ; enddo ; enddo
      !$OMP END DO nowait
    enddo
    !$OMP END PARALLEL

  end subroutine compute_source
          
  
  subroutine compute_source_ (cond_prim, cond_tau, gas_prim, gas_R, gas_gam, gas_k, gas_mu, source, mat, nb)
    use ICE_Parameters_m, only: pi, sigma_SB, I4
    use ICE_Config_Types_m, only: condensed_phase_t, obj_time_scheme
    use ICE_Lib_Model
    use ICE_Lib_Drag
    use ICE_Lib_Heat
    use ICE_Lib_Evaporation, only: evaporation, blowingFactor
    use ICE_Lib_Properties, only: mat_rho, mat_psat
    implicit none
    real(kind=R8), dimension(:), intent(in)    :: cond_prim, gas_prim
    real(kind=R8),               intent(inout) :: cond_tau
    real(kind=R8),               intent(in)    :: gas_R, gas_gam, gas_k, gas_mu
    real(kind=R8), dimension(:), intent(inout) :: source
    type(condensed_phase_t),     intent(in)    :: mat
    integer(kind=I4),            intent(in)    :: nb

    integer(kind=I4) :: n, ng
    real(kind=R8)    :: Rp, Re, Ma, Pr, Tr
    real(kind=R8)    :: Cd, Nu
    real(kind=R8)    :: rho_mat, mdot, Qconv, Qevap
    logical          :: override_Qdot
    real(kind=R8)    :: force(6)

    !> Slots of the closure: T is n-1, the number density n
    n = nb
    ng = size(gas_prim)

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

    !> Temperature ratio
    Tr = cond_prim(n-1)/gas_prim(5)

    !> Drag coefficient
    Cd = drag(Re,Ma,gas_gam,Tr,obj_time_scheme%dragSelect)

    !> Slip velocity
    force(2:4) = gas_prim(2:4) - cond_prim(2:4)

    !> Convective heat exchange of a single particle [W]
    Nu = heat(Re,Pr,Ma,obj_time_scheme%heatSelect)
    Qconv = 2._R8*Nu*gas_k*pi*Rp*(gas_prim(ng)-cond_prim(n-1))

    !> Mass exchange of a single particle [kg/s], negative while the particle
    !  evaporates. Qevap comes back in W; the latent sink is not included here,
    !  every closure's source_make adds it as -force(1)*lv_al. A Psat column in the
    !  property table replaces the Clausius-Clapeyron curve.
    if (mat%use_psat) then
      call evaporation(gas_prim(1), gas_prim(ng), gas_gam, gas_R, gas_mu, gas_k,  &
                       cond_prim(n-1), 2._R8*Rp, Re,                              &
                       mat%evapSelect, mat%intfSelect,                            &
                       mat%ep, mdot, Qevap, override_Qdot,              &
                       psatExt=mat_psat(mat, cond_prim(n-1)))
    else
      call evaporation(gas_prim(1), gas_prim(ng), gas_gam, gas_R, gas_mu, gas_k,  &
                       cond_prim(n-1), 2._R8*Rp, Re,                              &
                       mat%evapSelect, mat%intfSelect,                            &
                       mat%ep, mdot, Qevap, override_Qdot)
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
      
    !> Radiative heat exchange
    force(5) = force(5) + mat%emiss*sigma_SB * 2._R8*pi*Rp*Rp*cond_prim(n)*(gas_prim(ng)**4._I4-cond_prim(n-1)**4._I4)

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

  end subroutine compute_source_
    

end module ICE_Mod_Sources