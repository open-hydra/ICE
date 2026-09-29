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
    use ICE_Mod_MPI,        only: is_local_block, mpi_allreduce_min_r8, &
                                  mpi_allreduce_sum_r8_array, mpi_bcast_integer
    use ICE_Mod_Timers,     only: timer_source_begin, timer_source_end,   &
                                  timer_flux_begin,   timer_flux_end,     &
                                  timer_halo_begin,   timer_halo_end,     &
                                  timer_sync_begin,   timer_sync_end
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, p, srk
    real(R8)         :: average(5), dtlocal
    logical          :: endsim, iosim

    grid%iter                    = grid%iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    grid%dtglobal = 1e+5

    do p = 1, ngroups
      call assign_sound_make(p)
      !$omp parallel
      call compute_dt(p, obj_time_scheme%cfl, obj_time_scheme%cfl_rampa_iter, &
                      obj_time_scheme%dt_max, obj_time_scheme%tau_factor, grid)
      !$omp end parallel
    end do

    !> Global time step: smallest over all ranks
    call timer_sync_begin()
    call mpi_allreduce_min_r8(grid%dtglobal, dtlocal)
    call timer_sync_end()
    grid%dtglobal = dtlocal
    if (obj_time_scheme%time_accurate) grid%time = grid%time + grid%dtglobal

    !$omp parallel
    if (obj_time_scheme%time_accurate) call set_dt_global(grid)
    call state_copy(grid)
    !$omp end parallel

    do p = 1, ngroups

      call assign_all(p)
      if (len_trim(obj_time_scheme%riemann) == 0) then
        select case (trim(obj_time_scheme%model(p)))
        case ('MK');     call assign_riemann('Saurel', solid_of(p))
        case default;    call assign_riemann('Rusanov')
        end select
      else
        !> Saurel upwinds on the mean normal velocity and assembles the flux from the
        !  MK variable layout, so it is meaningless for the Gaussian closures.
        if (trim(obj_time_scheme%riemann) == 'Saurel' .and. &
            trim(obj_time_scheme%model(p)) /= 'MK') then
          write(*,'(A)') '  [ICE] the Saurel solver is specific to the MK closure; '// &
                         'use Rusanov or HLLE with IG and AG.'
          error stop 'ICE: Saurel with a non-MK family'
        endif
        call assign_riemann(trim(obj_time_scheme%riemann), solid_of(p))
      endif
      select case (trim(obj_time_scheme%solver_type))
      case ('euler'); nrk = 1
      case ('RK2');   nrk = 2
      case ('RK3');   nrk = 3
      case default
        write(*,'(A)') '  [ERROR] unknown time-scheme "'//trim(obj_time_scheme%solver_type)// &
                       '"; choose euler, RK2 or RK3.'
        error stop 'ICE: unknown time-scheme'
      end select

      do srk = 1, nrk

        call zero_residual(grid, p)

        if (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) then
          call timer_source_begin()
          call compute_source(grid, p)
          call timer_source_end()
        end if

        call timer_halo_begin()
        call compute_ghost(grid)
        call fill_second_ghost(grid)

        call compute_bound(grid, p)
        call timer_halo_end()

        call timer_flux_begin()
        call compute_flux(grid, p)

        call compute_residual(grid, p)

        if (obj_irs%enabled) call residual_smoothing(grid, p)

        call state_update(grid, p, srk)
        call timer_flux_end()

      end do

    end do

    ! Compute global residual (L2 norm of prim - prim_old, over all blocks and groups)
    obj_sim_param%residuotot = 0._R8
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      do p = 1, ngroups
        call Compute_Diagnostic(new=grid%blk(b)%cond_phase(p)%prim,     &
                              old=grid%blk(b)%cond_phase(p)%prim_old, &
                              dt=grid%blk(b)%cond_phase(p)%dt,        &
                              n=grid%blk(b)%dim, nc=ncond(p),         &
                              iT=nbase(p)-1,                           &
                              average=average,                         &
                              total=obj_sim_param%residuotot)
      end do
    end do
    call timer_sync_begin()
    call mpi_allreduce_sum_r8_array(obj_sim_param%residuotot, nres)
    call timer_sync_end()
    obj_sim_param%residuotot = sqrt(obj_sim_param%residuotot)

    iosim  = (mod(grid%iter, obj_io%sol_diter) == 0) &
         .or. (grid%time >= obj_sim_param%time_from_call + obj_io%sol_dtime)

    !> res-threshold = 0 means "never stop on the residual". Without it a transient
    !> whose density happens to be stationary -- a cloud relaxing in velocity only --
    !> reports a zero density residual and is declared converged at the first iteration.
    endsim = (obj_sim_param%iter_from_call >= grid%itermax)                                        &
         .or. (obj_multigrid%MG_level == 1 .and. obj_sim_param%res_threshold > 0._R8 .and.       &
               obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold)                       &
         .or. (grid%time >= obj_sim_param%time_threshold)

    if (endsim) then
      obj_sim_param%TODO = 3
    elseif (iosim) then
      obj_sim_param%TODO = 2
    else
      obj_sim_param%TODO = 1
    end if

    !> Every rank sees the same residual and time, but take the decision from root
    !> so that all ranks leave the time loop together whatever the rounding.
    call timer_sync_begin()
    call mpi_bcast_integer(obj_sim_param%TODO)
    call timer_sync_end()

  end subroutine Explicit_Step


end module ICE_Mod_Explicit
