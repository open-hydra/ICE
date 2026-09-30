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
    use ICE_Mod_Diagnostic, only: residual_norm
    use ICE_Mod_MPI,        only: is_local_block, mpi_allreduce_min_r8
    use ICE_Mod_Timers,     only: timer_source_begin, timer_source_end,   &
                                  timer_flux_begin,   timer_flux_end,     &
                                  timer_halo_begin,   timer_halo_end,     &
                                  timer_sync_begin,   timer_sync_end,     &
                                  timer_region_begin, timer_region_end,   &
                                  TR_DT, TR_COPY, TR_ZERO, TR_SOURCE, TR_GHOST1, TR_GHOST2, &
                                  TR_BOUND, TR_FLUX, TR_RESID, TR_IRS, TR_UPDATE, TR_DIAG
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: p, srk
    real(R8)         :: dtlocal
    logical          :: endsim, iosim, need_norm

    grid%iter                    = grid%iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    grid%dtglobal = 1e+5

    call timer_region_begin(TR_DT)
    do p = 1, ngroups
      call assign_wavespeed_make(p)
      !$omp parallel
      call compute_dt(p, obj_time_scheme%cfl, obj_time_scheme%cfl_rampa_iter, &
                      obj_time_scheme%dt_max, obj_time_scheme%tau_factor, grid)
      !$omp end parallel
    end do
    call timer_region_end(TR_DT)

    !> Global time step: smallest over all ranks. Only a time-accurate run reads
    !  it (set_dt_global, the time, the shell line); a steady run steps each cell
    !  by its own dt and would pay the collective for nothing.
    if (obj_time_scheme%time_accurate) then
      call timer_sync_begin()
      call mpi_allreduce_min_r8(grid%dtglobal, dtlocal)
      call timer_sync_end()
      grid%dtglobal = dtlocal
      grid%time = grid%time + grid%dtglobal
    end if

    call timer_region_begin(TR_COPY)
    !$omp parallel
    if (obj_time_scheme%time_accurate) call set_dt_global(grid)
    call state_copy(grid)
    !$omp end parallel
    call timer_region_end(TR_COPY)

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

        call timer_region_begin(TR_ZERO)
        call zero_residual(grid, p)
        call timer_region_end(TR_ZERO)

        if (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) then
          call timer_source_begin()
          call timer_region_begin(TR_SOURCE)
          call compute_source(grid, p)
          call timer_region_end(TR_SOURCE)
          call timer_source_end()
        end if

        call timer_halo_begin()
        call timer_region_begin(TR_GHOST1)
        call compute_ghost(grid)
        call timer_region_end(TR_GHOST1)
        call timer_region_begin(TR_GHOST2)
        call fill_second_ghost(grid)
        call timer_region_end(TR_GHOST2)

        call timer_region_begin(TR_BOUND)
        call compute_bound(grid, p)
        call timer_region_end(TR_BOUND)
        call timer_halo_end()

        call timer_flux_begin()
        call timer_region_begin(TR_FLUX)
        call compute_flux(grid, p)
        call timer_region_end(TR_FLUX)

        call timer_region_begin(TR_RESID)
        call compute_residual(grid, p)
        call timer_region_end(TR_RESID)

        if (obj_irs%enabled) then
          call timer_region_begin(TR_IRS)
          call residual_smoothing(grid, p)
          call timer_region_end(TR_IRS)
        end if

        call timer_region_begin(TR_UPDATE)
        call state_update(grid, p, srk)
        call timer_region_end(TR_UPDATE)
        call timer_flux_end()

      end do

    end do

    !> The residual norm, only on a step that reads it: the stop test, the
    !  residual file, the shell line, the last step. Every key can change at run
    !  time, so the conditions are evaluated here, every step.
    need_norm = (obj_multigrid%MG_level == 1 .and. obj_sim_param%res_threshold > 0._R8) &
           .or. due(grid%iter, obj_io%res_diter) .or. due(grid%iter, obj_io%shell_diter) &
           .or. (obj_sim_param%iter_from_call >= grid%itermax)                          &
           .or. (grid%time >= obj_sim_param%time_threshold)
    if (need_norm) then
      call timer_region_begin(TR_DIAG)
      call residual_norm(grid, obj_sim_param%residuotot)
      call timer_region_end(TR_DIAG)
    end if

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

    !> Every rank took this decision from the same numbers: the iteration count,
    !  the time after the dt minimum, and a residual norm summed in one order on
    !  every rank. Nothing needs to be broadcast.

  end subroutine Explicit_Step


  !> Iteration `iter` is one at which something happens every `every` steps.
  pure logical function due(iter, every)
    integer, intent(in) :: iter, every
    due = (every > 0) .and. (mod(iter, every) == 0)
  end function due


end module ICE_Mod_Explicit
