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


  !> input.ini read again during the run, every ini-diter iterations (ICE_postprocess, at
  !  iteration iter). The values keep to the rules of the registry, as at set-up; so does
  !  what set-up has sized and built on, the number of grid levels and the closure of a
  !  family, and riemann-solver keeps to the closures. An error of the reading (a file that
  !  cannot be used, a value that does not read or breaks its rule, such a change) stops the
  !  run with ini-error = stop, as at set-up. With ini-error = warn, the default, the reading
  !  is not taken: a warning gives the iteration and the error, and every parameter goes back
  !  to its value before the reading, those the reading had right too (Save_Values before
  !  Load_Ini, Restore_Values on the error), so that keys that go together never enter one
  !  without the other. The value of ini-error in force when the reading starts decides,
  !  also when the file read again writes another one. Every rank parses the text rank 0
  !  read (Open_Ini), so every rank takes the same decision. A time-accurate run keeps
  !  implicit residual smoothing off: irs = true is ignored, with a warning, rather than
  !  ending a long run on an edit.
  subroutine Read_Inifile_Runtime(iter)
    use Finer,              only: file_ini
    use ICE_Backend_INI,    only: Open_Ini, Load_Ini
    use ICE_Assign_Setup,   only: Check_Riemann_Solver
    use ICE_Config_Types_m, only: obj_multigrid, obj_time_scheme, obj_irs, obj_io
    use ICE_Mod_MPI,        only: mpi_is_root
    use ICE_Parameters_m,   only: codename, llen
    use IR_Precision,       only: str
    use ICE_Input_Registry
    implicit none
    integer, intent(in) :: iter
    type(file_ini)      :: fini
    character(len=1024) :: out
    character(len=512)  :: failure
    character(len=llen), allocatable :: closure(:)
    integer :: levels, p
    logical :: warn

    warn    = trim(obj_io%ini_error) == 'warn'
    levels  = obj_multigrid%MGL
    closure = obj_time_scheme%model

    if (warn) then
      call Open_Ini(fini, failure)
      if (trim(failure) /= '') then
        call Not_Taken(iter, failure)
        return
      end if
      call Save_Values()
    else
      call Open_Ini(fini)
    end if
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
    ! riemann-solver takes effect at the next step
    if (warn .and. trim(out) == '') call Check_Riemann_Solver(out)
    if (trim(out) /= '') then
      if (warn) then
        call Restore_Values()
        call Not_Taken(iter, out)
        return
      end if
      write(*, '(A)') trim(out)
      error stop 1
    end if
    if (.not. warn) call Check_Riemann_Solver()

    if (obj_irs%enabled .and. obj_time_scheme%time_accurate) then
      if (mpi_is_root .and. .not. irs_warned) &
        write(*,'(A)') ' [WARNING] input.ini read again: irs = true ignored, implicit residual smoothing '// &
                       'needs time-accurate = false'
      irs_warned      = .true.
      obj_irs%enabled = .false.
    end if

  end subroutine Read_Inifile_Runtime


  !> The warning of a reading that is not taken (ini-error = warn), without the [ERROR] in
  !  front of the error: the run goes on. Written once (rank 0) and at once (flush).
  subroutine Not_Taken(iter, error)
    use ICE_Mod_MPI,     only: mpi_is_root
    use iso_fortran_env, only: output_unit
    implicit none
    integer,          intent(in) :: iter
    character(len=*), intent(in) :: error
    character(len=:), allocatable :: reason

    if (.not. mpi_is_root) return
    reason = trim(adjustl(error))
    if (index(reason, '[ERROR] ') == 1)     reason = reason(9:)
    if (index(reason, '[ERROR ICE] ') == 1) reason = reason(13:)
    write(*,'(A,I0,A)') ' [WARNING] input.ini read again at iteration ', iter, ' is not taken: '//reason
    write(*,'(A)') '           every parameter keeps its value from before the reading (ini-error = warn)'
    flush(output_unit)

  end subroutine Not_Taken

end module ICE_Read_Ini
