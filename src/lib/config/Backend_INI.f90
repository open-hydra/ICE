module ICE_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini

contains

  subroutine Load_Ini(fini)
    use ICE_Input_Registry, only: reg
    use Finer,              only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error

    do i = 1, reg%size
      error = 1

      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
      else if (associated(reg%params(i)%value%iarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%iarr, error=error)
      else if (associated(reg%params(i)%value%rarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%rarr, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.
    end do

  end subroutine Load_Ini


  subroutine Scan_Ini(fini, nprobes, probe_sections, nmgl, ngroups)
    use Finer,            only: file_ini
    use ICE_Parameters_m, only: codename, clen
    use IR_Precision,     only: str
    implicit none
    type(file_ini),                          intent(in)  :: fini
    integer,                                 intent(out) :: nprobes
    character(len=clen), allocatable,        intent(out) :: probe_sections(:)
    integer,                                 intent(out) :: nmgl
    integer,                                 intent(out) :: ngroups
    integer             :: error, p
    character(len=clen) :: val

    ! --- Number of multigrid levels ---
    nmgl = 1
    call fini%get(section_name=trim(codename)//'-Multigrid', option_name='levels', &
                  val=nmgl, error=error)

    ! --- Number of condensed-phase families ---
    ngroups = 0
    do
      val = ''
      call fini%get(section_name=trim(codename)//'-Family'//trim(str(.true., ngroups+1)), &
                    option_name='model', val=val, error=error)
      if (error /= 0 .or. trim(val) == '') exit
      ngroups = ngroups + 1
    end do

    ! --- Probes ---
    nprobes = 0
    do
      val = ''
      call fini%get(section_name=trim(codename)//'-Probes', &
                    option_name='probe'//trim(str(.true., nprobes+1)), val=val, error=error)
      if (error /= 0 .or. trim(val) == '') exit
      nprobes = nprobes + 1
    end do

    if (nprobes > 0) then
      allocate(probe_sections(nprobes))
      do p = 1, nprobes
        call fini%get(section_name=trim(codename)//'-Probes', &
                      option_name='probe'//trim(str(.true.,p)), val=probe_sections(p), error=error)
      end do
    else
      allocate(probe_sections(0))
    end if

  end subroutine Scan_Ini

end module ICE_Backend_INI
