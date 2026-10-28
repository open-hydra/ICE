module ICE_Wrap_Postprocess

  implicit none
  private
  public :: ICE_postprocess

  !> len=* so the length comes from the literal. Both _MG strings were declared one
  !> character short, which cut the closing parenthesis off the format and made every
  !> shell line on a multigrid run a "syntax error in format" abort.
  character(len=*), parameter, private :: A_shell_format    = "('ICE  | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(len=*), parameter, private :: B_shell_format    = "('ICE  | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"
  character(len=*), parameter, private :: A_shell_format_MG = "('ICE  Grid Level', i2, ' | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(len=*), parameter, private :: B_shell_format_MG = "('ICE  Grid Level', i2, ' | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"

contains

  subroutine ICE_postprocess(simulation)
    use iso_fortran_env,       only: R8 => real64
    use IR_Precision
    use ICE_Advanced_Types_m,  only: ICE_simulation_type
    use ICE_Config_Types_m
    use ICE_IO_Solution
    use ICE_IO_Probes,         only: Write_Probes_Data, Close_Probes, nprobes
    use ICE_Mod_Diagnostic
    use ICE_Mod_Multigrid,     only: Prolongation
    use ICE_Read_Ini,          only: Read_Inifile_Runtime
    use ICE_Mod_MPI,           only: mpi_is_root
    use ICE_Mod_Timers,        only: timer_report, timer_summary
    use ICE_Mod_Inflow400,     only: report_400_inflow
    implicit none
    type(ICE_simulation_type), intent(inout) :: simulation
    character(llen) :: solfile, bckfile, dgsfile
    real(R8)        :: sim_time
    integer         :: level
    logical         :: iter_based

    level = obj_multigrid%MG_level

    ! Update ODP solutiontime
    if (obj_time_scheme%time_accurate) then
      simulation%ODP(1)%solutiontime = simulation%domain(1)%time
    else
      simulation%ODP(1)%solutiontime = -real(obj_sim_param%iter_general, 8)
    end if

    if (level == 1) then
      solfile = trim(ICE_phase_prefix)//'field'
      dgsfile = trim(ICE_phase_prefix)//'residuals'
      bckfile = trim(ICE_phase_prefix)//'field'
    else
      solfile = trim(ICE_phase_prefix)//'field-level'//trim(str(.true.,level))
      dgsfile = trim(ICE_phase_prefix)//'residuals-level'//trim(str(.true.,level))
      bckfile = trim(ICE_phase_prefix)//'field-level'//trim(str(.true.,level))
    end if

    ! END OF SIMULATION
    if (obj_sim_param%TODO == 3) then

      ! The last window, before the final write so that the I/O stays out of
      ! it. Collective, so every rank calls it.
      call timer_report(simulation%domain(1)%iter)

      ! The 400 faces through which particles entered, over the run. Collective, so every
      ! rank calls it.
      call report_400_inflow(simulation%domain(1)%nb)

      if (mpi_is_root) then
        write(obj_io%unitRES, trim(obj_io%unitRES_format)) obj_sim_param%iter_general, &
              simulation%domain(1)%time, obj_sim_param%residuotot
        flush(obj_io%unitRES)

        if (obj_time_scheme%time_accurate) then
          if (obj_multigrid%MGL > 1) then
            write(*,B_shell_format_MG) level, obj_sim_param%iter_from_call, &
                  simulation%domain(1)%time, simulation%domain(1)%dtglobal
          else
            write(*,B_shell_format) obj_sim_param%iter_from_call, &
                  simulation%domain(1)%time, simulation%domain(1)%dtglobal
          end if
        else
          if (obj_multigrid%MGL > 1) then
            write(*,A_shell_format_MG) level, obj_sim_param%iter_from_call, &
                  obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          else
            write(*,A_shell_format) obj_sim_param%iter_from_call, &
                  obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          end if
        end if

        call Cpu_Time(obj_sim_param%cputime(2))
        sim_time = (obj_sim_param%cputime(2) - obj_sim_param%cputime(1)) / obj_sim_param%nthreads
        write(*,*)
        write(*,*) '  Time of operation was', sim_time/60, 'min'
      end if

      ! What an output iteration writes, when the last iteration is one: the
      ! probe rows due and the solution due (sol-diter, else sol-dtime) with the
      ! next output number; with sol-overwrite = true the final solution below
      ! is that output, written once. Every rank evaluates the same condition.
      ! A run ends on the finest level, whose files these are.
      if (level == 1) then
        if (nprobes > 0) call Write_Probes_Data(obj_sim_param%iter_general, simulation%domain(1)%time)
        iter_based = mod(simulation%domain(1)%iter, obj_io%sol_diter) == 0
        if (iter_based .or. simulation%domain(1)%time >= obj_sim_param%time_from_call + obj_io%sol_dtime) then
          obj_io%output_number = obj_io%output_number + 1
          if (.not. iter_based) obj_sim_param%time_from_call = simulation%domain(1)%time
          if (.not. obj_io%sol_overwrite) then
            if (mpi_is_root) then
              if (iter_based) then
                write(*,*) ' ... writing iter-based solution'
              else
                write(*,*) ' ... writing time-based solution'
              end if
            end if
            call Write_Solution(simulation%domain(1), simulation%ODP(1), &
                                trim(solfile)//trim(str(.true.,obj_io%output_number)), obj_io%sol_fmt)
            if (.not. obj_time_scheme%time_accurate) &
              call Write_Diagnostic(simulation%domain(1), simulation%ODP(1), &
                                    trim(dgsfile)//trim(str(.true.,obj_io%output_number)))
          end if
        end if
      end if

      call Write_Solution(simulation%domain(1), simulation%ODP(1), solfile, obj_io%sol_fmt)
      if (.not. obj_time_scheme%time_accurate) &
        call Write_Diagnostic(simulation%domain(1), simulation%ODP(1), dgsfile)

      call timer_summary()

      ! The history files hold their last rows: close them here, before MPI is
      ! finalized, not at the end of the process.
      call Close_Probes()
      if (mpi_is_root) close(obj_io%unitRES)

    ! INTERMEDIATE SOLUTION EVALUATION
    elseif (obj_sim_param%TODO == 2) then

      if (level == 1) then
        obj_io%output_number = obj_io%output_number + 1
        if (.not.obj_io%sol_overwrite) solfile = trim(solfile)//trim(str(.true.,obj_io%output_number))
        if (.not.obj_io%sol_overwrite) dgsfile = trim(dgsfile)//trim(str(.true.,obj_io%output_number))
      end if

      ! Probes: on the finest level, whose cells they read, at the global
      ! iteration count of the residual history
      if (nprobes > 0 .and. level == 1) &
        call Write_Probes_Data(obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Collective, so every rank calls it -- not inside the root guard below.
      if (mod(simulation%domain(level)%iter, obj_io%shell_diter) == 0d0) &
        call timer_report(simulation%domain(level)%iter)

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%res_diter) == 0d0) then
          write(obj_io%unitRES, trim(obj_io%unitRES_format)) obj_sim_param%iter_general, &
                simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        end if

        if (mod(simulation%domain(level)%iter, obj_io%shell_diter) == 0d0) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1) then
              write(*,B_shell_format_MG) level, obj_sim_param%iter_from_call, &
                    simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write(*,B_shell_format) obj_sim_param%iter_from_call, &
                    simulation%domain(1)%time, simulation%domain(1)%dtglobal
            end if
          else
            if (obj_multigrid%MGL > 1) then
              write(*,A_shell_format_MG) level, obj_sim_param%iter_from_call, &
                    obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write(*,A_shell_format) obj_sim_param%iter_from_call, &
                    obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            end if
          end if
        end if
      end if

      ! Iter-based solution
      if (mod(simulation%domain(level)%iter, obj_io%sol_diter) == 0d0) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based solution'
        call Write_Solution(simulation%domain(1), simulation%ODP(1), solfile, obj_io%sol_fmt)
        if (.not. obj_time_scheme%time_accurate) &
          call Write_Diagnostic(simulation%domain(1), simulation%ODP(1), dgsfile)
      ! Time-based solution
      elseif (simulation%domain(level)%time >= obj_sim_param%time_from_call + obj_io%sol_dtime) then
        obj_sim_param%time_from_call = simulation%domain(level)%time
        if (mpi_is_root) write(*,*) ' ... writing time-based solution'
        call Write_Solution(simulation%domain(1), simulation%ODP(1), solfile, obj_io%sol_fmt)
        if (.not. obj_time_scheme%time_accurate) &
          call Write_Diagnostic(simulation%domain(1), simulation%ODP(1), dgsfile)
      end if

      ! MG level switch: change_MG set by Wrap_Solve after coarse level completes
      if (obj_multigrid%change_MG) then
        obj_multigrid%MG_level       = obj_multigrid%MG_level - 1
        obj_sim_param%iter_from_call = 0
        obj_sim_param%res_reference  = 0._R8
        obj_sim_param%time_from_call = 0.0_R8
        obj_multigrid%change_MG      = .false.
      end if

    ! AUXILIARY OPERATIONS AT EACH ITERATION
    elseif (obj_sim_param%TODO == 1) then

      ! Probes, as at an output iteration
      if (nprobes > 0 .and. level == 1) &
        call Write_Probes_Data(obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Collective, so every rank calls it -- not inside the root guard below.
      if (mod(simulation%domain(level)%iter, obj_io%shell_diter) == 0d0) &
        call timer_report(simulation%domain(level)%iter)

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%res_diter) == 0d0) then
          write(obj_io%unitRES, obj_io%unitRES_format) obj_sim_param%iter_general, &
                simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        end if

        if (mod(simulation%domain(level)%iter, obj_io%shell_diter) == 0d0) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1) then
              write(*,B_shell_format_MG) level, obj_sim_param%iter_from_call, &
                    simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write(*,B_shell_format) obj_sim_param%iter_from_call, &
                    simulation%domain(1)%time, simulation%domain(1)%dtglobal
            end if
          else
            if (obj_multigrid%MGL > 1) then
              write(*,A_shell_format_MG) level, obj_sim_param%iter_from_call, &
                    obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write(*,A_shell_format) obj_sim_param%iter_from_call, &
                    obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            end if
          end if
        end if
      end if

    end if

    ! Runtime INI update. iter-threshold and level<n>-iter take effect at once: the stop
    ! test (Mod_Explicit) and the level switch (Wrap_Solve) read each level's itermax,
    ! set from them here as at set-up (Setup_Multigrid, ICE_setup).
    if (mod(simulation%domain(level)%iter, obj_io%ini_diter) == 0d0) then
      call Read_Inifile_Runtime(simulation%domain(level)%iter)
      block
        integer :: m
        do m = 2, obj_multigrid%MGL
          simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
        end do
        if (obj_multigrid%MGL > 1) then
          simulation%domain(1)%itermax = min(obj_sim_param%iter_threshold, obj_multigrid%iter_threshold(1))
        else
          simulation%domain(1)%itermax = obj_sim_param%iter_threshold
        end if
      end block
    end if

  end subroutine ICE_postprocess


end module ICE_Wrap_Postprocess
