module ICE_Mod_Explicit
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Explicit_Step

contains


  subroutine Explicit_Step(grid)
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m
    use ICE_Global_m
    use ICE_Lib_Model
    use ICE_Lib_Riemann,  only: assign_riemann
    use ICE_Mod_dt
    use ICE_Lib_Newstate
    use ICE_Mod_Sources
    use ICE_Mod_BC_Fluxes
    use ICE_Mod_Fluxes
    use ICE_Lib_Residual
    use ICE_Lib_IRS,        only: residual_smoothing
    use ICE_Mod_Diagnostic, only: Compute_Diagnostic
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, p, srk, ios
    real(R8)         :: average(5)
    logical          :: endsim, iosim

    grid%iter                    = grid%iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    !$omp parallel

    grid%dtglobal = 1e+5
    do p = 1, ngroups
      call assign_sound_make(p)
      call compute_dt(p, obj_time_scheme%cfl, obj_time_scheme%cfl_rampa_iter, grid)
    end do

    if (obj_time_scheme%time_accurate) then
      call set_dt_global(grid)
      !$omp master
      grid%time = grid%time + grid%dtglobal
      !$omp end master
    end if

    call state_copy(grid)

    !$omp end parallel

    do p = 1, ngroups

      call assign_all(p)
      select case (trim(obj_time_scheme%model(p)))
      case ('MK');     call assign_riemann('Saurel')
      case default;    call assign_riemann('Rusanov')
      end select
      read(obj_time_scheme%solver_type, *, iostat=ios) nrk

      do srk = 1, nrk

        call zero_residual(grid, p)

        if (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) &
          call compute_source(grid, p)

        call compute_ghost(grid)
        call fill_second_ghost(grid)

        call compute_bound(grid)

        call compute_flux(grid, p)

        call compute_residual(grid, p)

        if (obj_irs%enabled) call residual_smoothing(grid, p)

        call state_update(grid, p, srk)

      end do

    end do

    ! Compute global residual (L2 norm of prim - prim_old, over all blocks and groups)
    obj_sim_param%residuotot = 0._R8
    do b = 1, grid%nb
      do p = 1, ngroups
        call Compute_Diagnostic(new=grid%blk(b)%cond_phase(p)%prim,     &
                              old=grid%blk(b)%cond_phase(p)%prim_old, &
                              dt=grid%blk(b)%cond_phase(p)%dt,        &
                              n=grid%blk(b)%dim, nc=ncond(p),         &
                              average=average,                         &
                              total=obj_sim_param%residuotot)
      end do
    end do
    obj_sim_param%residuotot = sqrt(obj_sim_param%residuotot)

    iosim  = (mod(grid%iter, obj_io%sol_diter) == 0) &
         .or. (mod(grid%iter, obj_io%bck_diter) == 0) &
         .or. (grid%time >= obj_sim_param%time_from_call + obj_io%sol_dtime)

    endsim = (obj_sim_param%iter_from_call >= grid%itermax)                                        &
         .or. (obj_multigrid%MG_level == 1 .and.                                                 &
               obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold)                       &
         .or. (grid%time >= obj_sim_param%time_threshold)

    if (endsim) then
      obj_sim_param%TODO = 3
    elseif (iosim) then
      obj_sim_param%TODO = 2
    else
      obj_sim_param%TODO = 1
    end if

  end subroutine Explicit_Step


end module ICE_Mod_Explicit
