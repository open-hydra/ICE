module ICE_Wrap_Postprocess

  implicit none
  private
  public :: ICE_postprocess

  integer,       private :: id_stampa = 0
  character(80), private :: A_shell_format    = "('ICE  | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(68), private :: B_shell_format    = "('ICE  | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"
  character(95), private :: A_shell_format_MG = "('ICE  Grid Level', i2, ' | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(83), private :: B_shell_format_MG = "('ICE  Grid Level', i2, ' | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"

contains

  subroutine ICE_postprocess(simulation)
    use iso_fortran_env,       only: R8 => real64
    use IR_Precision
    use ICE_Advanced_Types_m,  only: ICE_simulation_type
    use ICE_Config_Types_m
    use ICE_IO_Solution
    use ICE_IO_Probes,         only: Write_Probes_Data, nprobes
    use ICE_Mod_Diagnostic
    use ICE_Mod_Multigrid,     only: Prolongation
    use ICE_Read_Ini,          only: Read_Inifile_Runtime
    use ICE_Mod_MPI,           only: mpi_is_root
    implicit none
    type(ICE_simulation_type), intent(inout) :: simulation
    character(llen) :: solfile, bckfile, dgsfile
    real(R8)        :: sim_time
    integer         :: level

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

      if (mpi_is_root) then
        write(obj_io%unitRES, trim(obj_io%unitRES_format)) obj_sim_param%iter_general, &
              simulation%domain(1)%time, obj_sim_param%residuotot

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

      call Write_Solution(simulation%domain(1), simulation%ODP(1), solfile, obj_io%sol_fmt)
      if (.not. obj_time_scheme%time_accurate) &
        call Write_Diagnostic(simulation%domain(1), simulation%ODP(1), dgsfile)

    ! INTERMEDIATE SOLUTION EVALUATION
    elseif (obj_sim_param%TODO == 2) then

      if (level == 1) then
        id_stampa = id_stampa + 1
        if (.not.obj_io%sol_overwrite) solfile = trim(solfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%sol_overwrite) dgsfile = trim(dgsfile)//trim(str(.true.,id_stampa))
      end if

      ! Probes
      if (nprobes > 0) call Write_Probes_Data(simulation%domain(1)%iter, simulation%domain(1)%time)

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%res_diter) == 0d0) &
          write(obj_io%unitRES, trim(obj_io%unitRES_format)) obj_sim_param%iter_general, &
                simulation%domain(level)%time, obj_sim_param%residuotot

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
        obj_sim_param%time_from_call = 0.0_R8
        obj_multigrid%change_MG      = .false.
      end if

    ! AUXILIARY OPERATIONS AT EACH ITERATION
    elseif (obj_sim_param%TODO == 1) then

      ! Probes
      if (nprobes > 0) call Write_Probes_Data(simulation%domain(1)%iter, simulation%domain(1)%time)

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%res_diter) == 0d0) &
          write(obj_io%unitRES, obj_io%unitRES_format) obj_sim_param%iter_general, &
                simulation%domain(level)%time, obj_sim_param%residuotot

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

    ! Runtime INI update
    if (mod(simulation%domain(level)%iter, obj_io%ini_diter) == 0d0) &
      call Read_Inifile_Runtime()

  end subroutine ICE_postprocess


end module ICE_Wrap_Postprocess
