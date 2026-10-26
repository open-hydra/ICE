module ICE_Backend_INI
  implicit none
  private
  public :: Open_Ini, Load_Ini, Scan_Ini, Check_Unknown_Keys

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
  !
  !  Rank 0 alone reads the file and passes the text to every rank, which parses it:
  !  ranks reading it themselves could take different values from a file rewritten
  !  during the run, and a value that decides a collective (sol-diter, res-diter,
  !  ini-diter) would then make them call their collectives in different orders. A
  !  file that cannot be opened or read, that holds no key, or that FiNeR cannot load
  !  stops every rank, where the run went on with the defaults (at a re-read, with the
  !  previous values) or stopped later on the first required key.
  subroutine Open_Ini(fini)
    use Finer,           only: file_ini
    use ICE_Mod_MPI,     only: mpi_is_root, mpi_bcast_integer, mpi_bcast_chars
    use iso_fortran_env, only: output_unit
    implicit none
    type(file_ini), intent(inout) :: fini
    character(len=*), parameter   :: comments = '!;#'
    character(len=:), allocatable :: source
    character(len=1024)           :: line
    character(len=256)            :: msg
    integer :: unitfile, ios, state, error
    logical :: exists

    ! Rank 0 reads; state 0 read, 1 missing, 2 cannot be opened, 3 read error, 4 no key
    state  = 0
    source = ''
    if (mpi_is_root) then
      inquire(file='input.ini', exist=exists)
      if (.not. exists) then
        state = 1
        write(*,'(A)') '[ERROR ICE] input.ini not found in the working directory.'
      else
        open(newunit=unitfile, file='input.ini', status='old', action='read', iostat=ios, iomsg=msg)
        if (ios /= 0) then
          state = 2
          write(*,'(A)') '[ERROR] input.ini cannot be opened: '//trim(msg)
        else
          do
            read(unitfile, '(A)', iostat=ios, iomsg=msg) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (scan(adjustl(line), comments) == 1) cycle
            source = source//trim(line)//new_line('a')
          end do
          close(unitfile)
          if (ios > 0) then
            state = 3
            write(*,'(A)') '[ERROR] input.ini cannot be read: '//trim(msg)
          else if (len(source) == 0) then
            state = 4
            write(*,'(A)') '[ERROR] input.ini holds no key: it is empty or holds comment lines only'
          end if
        end if
      end if
      if (state /= 0) flush(output_unit)
    end if
    call mpi_bcast_integer(state)
    if (state /= 0) error stop 1

    ! Every rank parses the text of rank 0, so every rank finds the same values and errors
    call mpi_bcast_chars(source)
    call fini%load(source=source, error=error)
    if (error /= 0) then
      write(*,'(A,I0,A)') '[ERROR] input.ini cannot be read by the INI reader (FiNeR error ', error, &
                          '; a line FiNeR could not parse is named above, if any)'
      error stop 1
    end if

  end subroutine Open_Ini

  !> The values of input.ini into the variables of the registry. FiNeR leaves a parameter as
  !  it was and returns err_option_vals for a value that does not convert to its type (a word,
  !  NaN, Infinity, a decimal comma, a logical not starting with T or F), a list longer than
  !  the parameter holds, or a key with no value: such a key is marked unread, with its text,
  !  and Validate_Registry stops on it instead of the run going on with the default. So is an
  !  integer written with a period or an exponent (gfortran refuses 2.7, ifx reads it as 2)
  !  and a number or logical written with several values, of which FiNeR takes the first.
  subroutine Load_Ini(fini)
    use ICE_Input_Registry, only: reg
    use Finer,              only: file_ini, err_option_vals
    implicit none
    type(file_ini), intent(in) :: fini
    character(len=1024) :: buf
    character(len=:), allocatable :: text
    integer :: i, error, error_text

    do i = 1, reg%size
      error = 1

      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        reg%params(i)%multi = .false.
        if (reg%params(i)%per_material) then
          buf = ''
          call fini%get(reg%params(i)%section, reg%params(i)%name, val=buf, error=error)
          reg%params(i)%multi = error == 0 .and. count_fields(buf) > 1
        endif
        if (.not. reg%params(i)%multi) &
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

      reg%params(i)%unread = (error == err_option_vals)
      if (allocated(reg%params(i)%text)) deallocate(reg%params(i)%text)
      call fini%get_string(reg%params(i)%section, reg%params(i)%name, val=text, error=error_text)
      if (error_text == 0) reg%params(i)%text = text
      if (error == 0 .and. error_text == 0) then
        if (associated(reg%params(i)%value%i) .or. associated(reg%params(i)%value%iarr)) then
          if (scan(text, '.eEdDqQ') > 0) reg%params(i)%unread = .true.
        end if
        if (associated(reg%params(i)%value%i) .or. associated(reg%params(i)%value%l) .or. &
            (associated(reg%params(i)%value%r) .and. .not. reg%params(i)%multi)) then
          if (count_fields(text) > 1) reg%params(i)%unread = .true.
        end if
      end if
    end do

  end subroutine Load_Ini


  !> Blank-separated fields of a value.
  pure integer function count_fields(line)
    character(len=*), intent(in) :: line
    integer :: k
    logical :: inside
    count_fields = 0
    inside = .false.
    do k = 1, len_trim(line)
      if (line(k:k) == ' ') then
        inside = .false.
      elseif (.not. inside) then
        inside = .true.
        count_fields = count_fields + 1
      endif
    enddo
  end function count_fields


  !> Report keys ICE does not know, inside sections it does, and the sections it would skip.
  !
  !  An unrecognised key is not an error to FiNeR and not an error to the registry
  !  either: it simply never reaches a parameter, which takes its default instead.
  !  A renamed or mistyped key then changes nothing and says nothing, which is the
  !  worst way for an input file to be wrong. Sections ICE does not own are left
  !  alone, so a case may keep the mesh-generator sections that produced it. A
  !  section named after ICE that it does not read (a renamed one such as
  !  [ICE-Scheme], a family section without closure) would be skipped whole with
  !  its keys, and so would the second block of a section given twice, since FiNeR
  !  hands out the first: both are reported too. An [ICE-*] section without keys
  !  changes nothing and is left alone.
  subroutine Check_Unknown_Keys(fini)
    use ICE_Input_Registry, only: reg
    use ICE_Parameters_m,   only: codename
    use Finer,              only: file_ini
    implicit none
    type(file_ini), intent(inout) :: fini
    character(len=:), allocatable :: sections(:)
    character(len=:), allocatable :: pairs(:)
    integer :: s, k, nbad, nsec, nkeys
    logical :: known

    call fini%get_sections_list(sections)
    if (.not. allocated(sections)) return

    nsec = 0
    do s = 1, size(sections)
      if (any(sections(1:s-1) == sections(s))) cycle
      if (registered(sections(s))) then
        if (count(sections == sections(s)) < 2) cycle
        if (nsec == 0) write(*,'(A)') '  [ERROR] input.ini has sections ICE would skip:'
        nsec = nsec + 1
        write(*,'(A,I0,A)') '    ['//trim(sections(s))//'] given ', count(sections == sections(s)), &
                            ' times: only the first is read'
      else if (index(sections(s), trim(codename)//'-') == 1) then
        nkeys = 0
        do while (fini%loop(section_name=trim(sections(s)), option_pairs=pairs))
          nkeys = nkeys + 1
        end do
        if (nkeys == 0) cycle
        if (nsec == 0) write(*,'(A)') '  [ERROR] input.ini has sections ICE would skip:'
        nsec = nsec + 1
        write(*,'(A)') '    ['//trim(sections(s))//'] is not a section ICE reads'
      end if
    end do

    nbad = 0
    do s = 1, size(sections)
      !> Only sections the registry knows: anything else belongs to another tool.
      if (any(sections(1:s-1) == sections(s))) cycle
      known = registered(sections(s))
      if (.not. known) cycle

      do while (fini%loop(section_name=trim(sections(s)), option_pairs=pairs))
        known = .false.
        do k = 1, reg%size
          if (trim(reg%params(k)%section) == trim(sections(s)) .and. &
              trim(reg%params(k)%name)    == trim(pairs(1))) then
            known = .true. ; exit
          end if
        end do
        if (.not. known) then
          if (nbad == 0) write(*,'(A)') '  [ERROR] input.ini has keys ICE does not know:'
          nbad = nbad + 1
          write(*,'(A)') '    ['//trim(sections(s))//'] '//trim(pairs(1))
        end if
      end do
    end do

    if (nbad > 0 .or. nsec > 0) then
      write(*,'(A)') '  They would be ignored and their parameters left at the default.'
      write(*,'(A)') '  See docs/user/registry.md for the current names.'
      error stop 1
    end if

  contains

    !> A section the registry holds keys of.
    logical function registered(name)
      character(len=*), intent(in) :: name
      integer :: k
      registered = .false.
      do k = 1, reg%size
        if (trim(reg%params(k)%section) == trim(name)) then
          registered = .true. ; return
        end if
      end do
    end function registered

  end subroutine Check_Unknown_Keys


  subroutine Scan_Ini(fini, nprobes, probe_sections, nmgl, ngroups)
    use Finer,            only: file_ini, err_option
    use ICE_Parameters_m, only: codename, clen
    use IR_Precision,     only: str
    implicit none
    type(file_ini),                          intent(in)  :: fini
    integer,                                 intent(out) :: nprobes
    character(len=clen), allocatable,        intent(out) :: probe_sections(:)
    integer,                                 intent(out) :: nmgl
    integer,                                 intent(out) :: ngroups
    integer             :: error, error_text, p
    character(len=clen) :: val
    character(len=:), allocatable :: text

    ! --- Number of multigrid levels ---
    ! levels sizes the registry (one level<n>-iter per level) before the keys are matched
    ! against it, so a value that does not read as one integer stops here, as Load_Ini
    ! would stop it: gfortran would leave one level, ifx reads 2.0 as two.
    nmgl = 1
    call fini%get(section_name=trim(codename)//'-Multigrid', option_name='levels', &
                  val=nmgl, error=error)
    if (error /= err_option) then
      text = ''
      call fini%get_string(section_name=trim(codename)//'-Multigrid', option_name='levels', &
                           val=text, error=error_text)
      if (error /= 0 .or. scan(text, '.eEdDqQ') > 0 .or. count_fields(text) /= 1) then
        if (len_trim(text) == 0) then
          write(*,'(A)') '[ERROR] input.ini, ['//trim(codename)//'-Multigrid] levels =: the value is empty'
        else
          write(*,'(A)') '[ERROR] input.ini, ['//trim(codename)//'-Multigrid] levels = '//trim(text)// &
                         ': the value is not an integer'
        end if
        error stop 1
      end if
    end if

    ! --- Number of condensed-phase families ---
    ngroups = 0
    do
      val = ''
      call fini%get(section_name=trim(codename)//'-Family'//trim(str(.true., ngroups+1)), &
                    option_name='closure', val=val, error=error)
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
