module ICE_Backend_INI
  implicit none
  private
  public :: Open_Ini, Load_Ini, Scan_Ini

contains

  !> Read input.ini and hand FiNeR the text with the full-line comments removed.
  !
  !  FiNeR's section sanitiser treats any line without the option separator as the
  !  continuation of the option above it and appends it there. A comment line after
  !  the last option of a section therefore ends up inside that option's value:
  !  `shell-diter = 100` followed by `# note` is read as the value `100 # note`.
  !  Only `;` escapes this, because it doubles as the inline-comment delimiter and is
  !  trimmed off again afterwards. Stripping the comment lines here makes `!`, `;`
  !  and `#` behave the same, and costs nothing else.
  subroutine Open_Ini(fini)
    use Finer, only: file_ini
    implicit none
    type(file_ini), intent(inout) :: fini
    character(len=*), parameter   :: comments = '!;#'
    character(len=:), allocatable :: source
    character(len=1024)           :: line
    integer :: unitfile, ios
    logical :: exists

    inquire(file='input.ini', exist=exists)
    if (.not. exists) then
      write(*,'(A)') '  [ICE] input.ini not found in the working directory.'
      error stop 'ICE: no input.ini'
    endif

    source = ''
    open(newunit=unitfile, file='input.ini', status='old', action='read')
    do
      read(unitfile, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (scan(adjustl(line), comments) == 1) cycle
      source = source//trim(line)//new_line('a')
    end do
    close(unitfile)

    call fini%load(source=source)

  end subroutine Open_Ini

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
