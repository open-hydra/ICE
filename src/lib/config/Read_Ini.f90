module ICE_Read_Ini
  implicit none

  logical, private :: irs_warned = .false.   ! the warning of irs ignored at a re-read is given once

contains

  subroutine Read_Inifile()
    use Finer,               only: file_ini
    use ICE_Read_Sim_Param,  only: Register_Sim_Param
    use ICE_Read_IO,         only: Register_IO_Fields, Register_Probes
    use ICE_Read_Numerics,   only: Register_Numerics, Register_Families
    use ICE_Read_Physics,    only: Register_Physics
    use ICE_Backend_INI,     only: Open_Ini, Load_Ini, Scan_Ini, Check_Unknown_Keys
    use ICE_Input_Registry
    use ICE_Global_m,        only: ngroups
    use ICE_Parameters_m,    only: clen
    implicit none
    type(file_ini)                   :: fini
    character(len=1024)              :: out
    integer                          :: nprobes, nmgl
    character(len=clen), allocatable :: probe_sections(:)

    call Open_Ini(fini)

    call Scan_Ini(fini, nprobes, probe_sections, nmgl, ngroups)

    call Register_Sim_Param()
    call Register_IO_Fields()
    if (nprobes > 0) call Register_Probes(nprobes, probe_sections)
    call Register_Numerics(nmgl)
    call Register_Physics()
    call Register_Families(ngroups)

    call Check_Unknown_Keys(fini)
    call Load_Ini(fini)

    out = Validate_Registry()
    if (trim(out) /= '') then
      write(*, '(A)') trim(out)
      error stop 1
    end if

  end subroutine Read_Inifile


  !> input.ini read again during the run, every ini-diter iterations. The values keep to the
  !  rules of the registry, as at set-up, or the run stops; so does a change of what set-up
  !  has sized and built on, the number of grid levels and the closure of a family. A
  !  time-accurate run keeps implicit residual smoothing off: irs = true is ignored, with a
  !  warning, rather than ending a long run on an edit.
  subroutine Read_Inifile_Runtime()
    use Finer,              only: file_ini
    use ICE_Backend_INI,    only: Open_Ini, Load_Ini
    use ICE_Assign_Setup,   only: Check_Riemann_Solver
    use ICE_Config_Types_m, only: obj_multigrid, obj_time_scheme, obj_irs
    use ICE_Mod_MPI,        only: mpi_is_root
    use ICE_Parameters_m,   only: codename, llen
    use IR_Precision,       only: str
    use ICE_Input_Registry
    implicit none
    type(file_ini)      :: fini
    character(len=1024) :: out
    character(len=llen), allocatable :: closure(:)
    integer :: levels, p

    levels  = obj_multigrid%MGL
    closure = obj_time_scheme%model

    call Open_Ini(fini)
    call Load_Ini(fini)
    out = Validate_Registry()
    if (trim(out) == '' .and. obj_multigrid%MGL /= levels) &
      write(out, '(A,I0,A,I0,A)') '[ERROR] input.ini, ['//trim(codename)//'-Multigrid] levels = ', &
        obj_multigrid%MGL, ': the number of grid levels cannot change during the run (', levels, ' at set-up)'
    do p = 1, size(closure)
      if (trim(out) /= '') exit
      if (trim(obj_time_scheme%model(p)) /= trim(closure(p))) &
        out = '[ERROR] input.ini, ['//trim(codename)//'-Family'//trim(str(.true., p))//'] closure = '// &
              trim(obj_time_scheme%model(p))//': the closure cannot change during the run ('// &
              trim(closure(p))//' at set-up)'
    end do
    if (trim(out) /= '') then
      write(*, '(A)') trim(out)
      error stop 1
    end if

    ! riemann-solver takes effect at the next step
    call Check_Riemann_Solver()

    if (obj_irs%enabled .and. obj_time_scheme%time_accurate) then
      if (mpi_is_root .and. .not. irs_warned) &
        write(*,'(A)') ' [WARNING] input.ini read again: irs = true ignored, implicit residual smoothing '// &
                       'needs time-accurate = false'
      irs_warned      = .true.
      obj_irs%enabled = .false.
    end if

  end subroutine Read_Inifile_Runtime

end module ICE_Read_Ini
