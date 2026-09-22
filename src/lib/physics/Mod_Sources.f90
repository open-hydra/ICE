module ICE_Mod_Sources
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block

  implicit none
  private
  public :: compute_source

contains

  subroutine compute_source (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k

    !$OMP PARALLEL
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
  
      !$OMP DO COLLAPSE(3)
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
                              grid%blk(b)%cond_phase(p)%source(:,i,j,k)  )
  
      enddo ; enddo ; enddo
      !$OMP END DO

    enddo
    !$OMP END PARALLEL

  end subroutine compute_source
          
  
  subroutine compute_source_ (cond_prim, cond_tau, gas_prim, gas_R, gas_gam, gas_k, gas_mu, source)
    use ICE_Parameters_m, only: pi, sigma_SB, I4
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Lib_Model
    use ICE_Lib_Drag
    use ICE_Lib_Heat
    use ICE_Load_Table,     only: get_rho_al
    implicit none
    real(kind=R8), dimension(:), intent(in)    :: cond_prim, gas_prim
    real(kind=R8),               intent(inout) :: cond_tau
    real(kind=R8),               intent(in)    :: gas_R, gas_gam, gas_k, gas_mu
    real(kind=R8), dimension(:), intent(inout) :: source

    integer(kind=I4) :: n, ng
    real(kind=R8)    :: Rp, Re, Ma, Pr, Tr
    real(kind=R8)    :: B, Cd, Nu
    real(kind=R8)    :: force(6)

    !> Condensed variable number
    n = size(cond_prim)
    ng = size(gas_prim)

    !> Particles radius 
    Rp = (0.75_R8*cond_prim(1) / (cond_prim(n)*pi*get_rho_al(cond_prim(n-1))))**(1._R8/3._R8)

    !> Reynolds number
    Re = 2._R8*gas_prim(1)*Rp*norm2(gas_prim(2:4)-cond_prim(2:4))/gas_mu
    !Re = max(Re,1.e-15_R8)
  
    !> Mach number
    Ma = norm2(gas_prim(2:4)-cond_prim(2:4)) / (gas_gam*gas_R*gas_prim(5))**0.5_R8
    !Ma = max(Ma,1.e-15_R8)
    Ma = min(Ma,1.e+00_R8)

    !> Prandtl number
    Pr = gas_gam/(gas_gam-1._R8)*gas_R*gas_mu/gas_k

    !> Temperature ratio
    Tr = cond_prim(n-1)/gas_prim(5)

    !> Drag coefficient
    Cd = drag(Re,Ma,gas_gam,Tr)

    !> Slip velocity
    force(2:4) = gas_prim(2:4) - cond_prim(2:4)

    !> Mass and convective heat exchange
    !> No combustion
    Nu = heat(Re,Pr,Ma)
    force(1) = 0._R8
    force(5) = 2._R8*Nu*gas_k*pi*Rp*(gas_prim(ng)-cond_prim(n-1))*cond_prim(n)

    !> Combustion  
    !> B = (gas_gam/(gas_gam-1)*gas_R*(gas_prim(ng)-cond_prim(n-1))+obj_condensed%q_al)/obj_condensed%lv_al
    !> force(1) = 2._R8*pi*gas_mu/Pr * cond_prim(n)*Rp*log(1._R8+B) * &
    !>            2._R8*(1._R8+0.3_R8*Re**0.5_R8*Pr**0.333_R8)
    !> force(5) = force(1)*obj_condensed%lv_al                     
      
    !> Radiative heat exchange
    force(5) = force(5) + obj_condensed%emiss*sigma_SB * 2._R8*pi*Rp*Rp*cond_prim(n)*(gas_prim(ng)**4._I4-cond_prim(n-1)**4._I4)

    !> Particles relaxation time
    cond_tau = 8._R8*get_rho_al(cond_prim(n-1))*Rp / (3._R8*gas_prim(1)*Cd*norm2(gas_prim(2:4)-cond_prim(2:4)) + 1e-20)
    !> Use only for Vie validation test
    !> cond_tau = 5.0_R8
    force(6) = cond_tau

    !> Source terms
    source = source_make(cond_prim,force)

  end subroutine compute_source_
    

end module ICE_Mod_Sources