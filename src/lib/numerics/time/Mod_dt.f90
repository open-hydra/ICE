module ICE_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block

  implicit none
  private
  public :: compute_dt, set_dt_global

contains

  subroutine compute_dt (p, cfl, cfl_rampa_iter, dt_max, tau_factor, grid)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_sim_param
    implicit none
    integer(kind=I4), intent(in)  :: p
    real(kind=R8), intent(in)     :: cfl
    integer(kind=I4), intent(in)  :: cfl_rampa_iter
    real(kind=R8), intent(in)     :: dt_max
    real(kind=R8), intent(in)     :: tau_factor
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, i, j, k
    real(kind=R8)    :: dtmin, tau
    logical          :: tau_limit

    tau_limit = (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) .and. tau_factor > 0._R8

    !> Per-thread minimum, merged once below: updating grid%dtglobal inside the
    !> worksharing loop would race between threads, and so would reading it here
    !> while another thread merges (the blocks carry no barrier between them).
    dtmin = huge(dtmin)

    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$omp do collapse(3)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        grid%blk(b)%cond_phase(p)%dt(i,j,k) = 1e+5

        call compute_dt_ (grid%blk(b)%cond_phase(p)%prim(:,i,j,k), &
                          grid%blk(b)%dl(i,j,k)%c,                 &
                          grid%blk(b)%M(i,j,k)%c,                  &
                          grid%blk(b)%cond_phase(p)%dt(i,j,k))

         
        grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%blk(b)%cond_phase(p)%dt(i,j,k) * cfl

        if (grid%iter <= cfl_rampa_iter) then
          grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%blk(b)%cond_phase(p)%dt(i,j,k) * grid%iter / cfl_rampa_iter
        endif

        !> dt-max is a ceiling on the step actually taken, so it comes after the CFL
        !  factor and after the ramp. Its job is to bound the explicit source terms,
        !  whose relaxation times the CFL condition knows nothing about.
        grid%blk(b)%cond_phase(p)%dt(i,j,k) = min(grid%blk(b)%cond_phase(p)%dt(i,j,k), dt_max)

        !> Coupled runs: the explicit relaxation sources bound the step by the relaxation time
        !  of the cells that carry the phase; an empty cell has no meaningful tau
        if (tau_limit) then
          tau = grid%blk(b)%cond_phase(p)%tau(i,j,k)
          if (grid%blk(b)%cond_phase(p)%prim(1,i,j,k) > rho_empty .and. ieee_is_finite(tau) .and. tau > 0._R8) &
            grid%blk(b)%cond_phase(p)%dt(i,j,k) = min(grid%blk(b)%cond_phase(p)%dt(i,j,k), tau_factor*tau)
        endif

        dtmin = min (dtmin, grid%blk(b)%cond_phase(p)%dt(i,j,k))

      enddo ; enddo ; enddo
      !$omp end do nowait
    enddo

    !$omp critical (ice_dtglobal)
    grid%dtglobal = min (grid%dtglobal, dtmin)
    !$omp end critical (ice_dtglobal)
    !$omp barrier

  end subroutine compute_dt    
  
  
  subroutine compute_dt_ (prim, length, tensor, dt)
    use ICE_Lib_Model
    implicit none
    real(kind=R8), dimension(:),   intent(in)    :: prim
    real(kind=R8), dimension(:),   intent(in)    :: length
    real(kind=R8), dimension(:,:), intent(in)    :: tensor
    real(kind=R8),                 intent(inout) :: dt
    
    integer(kind=I4) :: d
    real(kind=R8)    :: versor(3)
    real(kind=R8)    :: speed, sound, dtd

    do d = 1, 3
      versor = tensor(d,:) / norm2 ( tensor(d,:) )
      speed  = abs( dot_product (prim(2:4), versor) )
      sound  = wavespeed_make(prim, versor)
      
      dtd = length(d) / (speed + sound)
      dt = min (dt,dtd)

    enddo    

  end subroutine compute_dt_

  
  subroutine set_dt_global (grid)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, p, i, j, k
  
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      do p = 1, ngroups

        !$omp do collapse(3)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)
        
          grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%dtglobal
        
        end do ; end do ; end do
        !$omp end do nowait
      enddo
    enddo  
  
  end subroutine set_dt_global
  

end module ICE_Mod_dt