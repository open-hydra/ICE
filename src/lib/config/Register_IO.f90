module ICE_Read_IO
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Input_Registry
  use ICE_Config_Types_m
  implicit none
  private
  public :: Register_IO_Fields, Register_Probes

contains

  subroutine Register_IO_Fields()
    use ICE_Parameters_m, only: codename
    implicit none
    character(len=:), allocatable :: section

    section = trim(codename)//'-IO'

    obj_io%warning_message = 'none'
    obj_io%error_message   = 'none'
    obj_io%description     = 'none'

    call reg%add(section, 'sol-format',    obj_io%sol_format,    &
                 'tecplot ascii', 'Solution output format (writer mode)', '', .false.)
    call reg%add(section, 'bck-format',    obj_io%bck_format,    &
                 'native binary',  'Backup output format (writer mode)',   '', .false.)

    call reg%add(section, 'sol-diter',     obj_io%sol_diter,     &
                 '1000000000', 'Solution output iteration frequency', '> 0', .false.)
    call reg%add(section, 'sol-dtime',     obj_io%sol_dtime,     &
                 '1e30',       'Solution output time frequency',      '> 0', .false.)
    call reg%add(section, 'sol-overwrite', obj_io%sol_overwrite, &
                 'true',       'Overwrite solution files',            'true, false', .false.)

    call reg%add(section, 'bck-diter',     obj_io%bck_diter,     &
                 '1000000000', 'Backup output iteration frequency', '> 0', .false.)
    call reg%add(section, 'bck-dtime',     obj_io%bck_dtime,     &
                 '1e30',       'Backup output time frequency',     '> 0', .false.)
    call reg%add(section, 'bck-overwrite', obj_io%bck_overwrite, &
                 'true',       'Overwrite backup files',           'true, false', .false.)

    call reg%add(section, 'shell-diter',   obj_io%shell_diter,   &
                 '10',         'Shell update iteration frequency',    '> 0', .false.)
    call reg%add(section, 'res-diter',     obj_io%res_diter,     &
                 '10',         'Residual history write frequency',    '> 0', .false.)
    call reg%add(section, 'ini-diter',     obj_io%ini_diter,     &
                 '1000000000', 'Runtime input.ini reload frequency', '> 0', .false.)
    call reg%add(section, 'gas-path',      obj_io%gaspath,       &
                 'INPUT/',     'Gas-phase solution path',         '',    .false.)

  end subroutine Register_IO_Fields


  subroutine Register_Probes(n, probes_name)
    use ICE_Parameters_m, only: codename
    use IR_Precision,     only: str
    implicit none
    integer,      intent(in) :: n
    character(*), intent(in) :: probes_name(n)
    integer :: p

    allocate(obj_io_probes(n))
    do p = 1, n
      call reg%add(trim(codename)//'-Probes', 'probe'//trim(str(.true.,p)), &
                   obj_io_probes(p)%file, 'probe-placeholder', 'Probe file name', '', .false.)
      call Register_One_Probe(p, probes_name(p))
    end do

  end subroutine Register_Probes


  subroutine Register_One_Probe(p, probe_name)
    implicit none
    integer,      intent(in) :: p
    character(*), intent(in) :: probe_name

    obj_io_probes(p)%warning_message = 'none'
    obj_io_probes(p)%error_message   = 'none'
    obj_io_probes(p)%description     = 'none'

    call reg%add(probe_name, 'variables',      obj_io_probes(p)%varnames, &
                 'none',         'Probe variables to write',          '', .false.)
    call reg%add(probe_name, 'dtime',          obj_io_probes(p)%dtime,    &
                 '1e30',         'Probe output time frequency',       '> 0', .false.)
    call reg%add(probe_name, 'diter',          obj_io_probes(p)%diter,    &
                 '1000000000',   'Probe output iter frequency',       '> 0', .false.)
    call reg%add(probe_name, 'index-position', obj_io_probes(p)%iloc,     &
                 '0 0 0 0',      'Probe location by index (b i j k)', '>= 0', .false.)
    call reg%add(probe_name, 'position',       obj_io_probes(p)%loc,      &
                 '0.0 0.0 0.0',  'Probe location by coordinates',     '', .false.)

  end subroutine Register_One_Probe

end module ICE_Read_IO
