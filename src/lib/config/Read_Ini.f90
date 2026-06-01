module ICE_Read_Ini
  implicit none

contains

  subroutine Read_Inifile()
    use Finer,               only: file_ini
    use ICE_Read_Sim_Param,  only: Register_Sim_Param
    use ICE_Read_IO,         only: Register_IO_Fields, Register_Probes
    use ICE_Read_Numerics,   only: Register_Numerics, Register_Families
    use ICE_Read_Physics,    only: Register_Physics
    use ICE_Backend_INI,     only: Load_Ini, Scan_Ini
    use ICE_Input_Registry
    use ICE_Global_m,        only: ngroups
    use ICE_Parameters_m,    only: clen
    implicit none
    type(file_ini)                   :: fini
    character(len=1024)              :: out
    integer                          :: nprobes, nmgl
    character(len=clen), allocatable :: probe_sections(:)

    call fini%load(filename='input.ini')

    call Scan_Ini(fini, nprobes, probe_sections, nmgl, ngroups)

    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Numerics(nmgl)
    call Register_Families(ngroups)
    call Register_Physics()
    if (nprobes > 0) call Register_Probes(nprobes, probe_sections)

    call Load_Ini(fini)

    out = Validate_Registry()
    if (trim(out) /= '') then
      write(*, '(A)') trim(out)
      error stop
    end if

  end subroutine Read_Inifile


  subroutine Read_Inifile_Runtime()
    use Finer,            only: file_ini
    use ICE_Backend_INI,  only: Load_Ini
    use ICE_Input_Registry
    implicit none
    type(file_ini)      :: fini
    character(len=1024) :: out

    call fini%load(filename='input.ini')
    call Load_Ini(fini)
    out = Validate_Registry()
    if (trim(out) /= '') then
      write(*, '(A)') trim(out)
      error stop
    end if

  end subroutine Read_Inifile_Runtime

end module ICE_Read_Ini
