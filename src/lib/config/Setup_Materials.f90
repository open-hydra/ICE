!>@brief The condensed materials. INPUT/<prefix>phase.txt names them, one line each,
!> "<name> <groups> [key=value ...]"; the families map onto them in that order, as ATLAS
!> orders the populations. Each material starts from the [ICE-Physics] values (one value
!> per material for the property keys), its tokens then override the model keys, as in
!> IGLOO, and the property table gives one zone per material.
module ICE_Setup_Materials
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: Setup_Materials

  !> The phase-file keys, the set IGLOO reads and ATLAS writes
  character(len=17), parameter :: word_keys(6) = [character(len=17) :: 'evaporation', 'liquid-conduction', &
                                  'interface', 'boiling', 'combustion', 'solidification']
  character(len=9),  parameter :: real_keys(14) = [character(len=9) :: 'alpha-e', 'k-liq', 'mu-liq', 'K-burn', &
                                  'n-burn', 'X-eff', 'beta-part', 'xi-cap', 'T-ign', 'q-comb', 'T-melt', 'h-fus', &
                                  'T-nuc', 'cp-solid']

  !> The [ICE-Physics] keys that may carry one value per material, and their range
  character(len=23), parameter :: vector_keys(10) = [character(len=23) :: 'density', 'specific-heat', &
                                  'latent-heat', 'emissivity', 'vapour-molar-mass', 'boiling-temperature', &
                                  'vapour-specific-heat', 'lewis-number', 'vapour-mass-fraction', 'evaporation-coefficient']
  logical, parameter :: vector_positive(10) = [.true., .true., .true., .false., .true., .true., .false., .true., &
                                               .false., .true.]
  !> The alias of each such key ('' = none), the name IGLOO also reads
  character(len=5),  parameter :: vector_alias(10) = [character(len=5) :: '', '', '', '', '', 'Tboil', '', '', '', '']

contains

  subroutine Setup_Materials()
    use ICE_Global_m,        only: ngroups, nmat, mat_of, ICE_phase_prefix, ncond, nbase, solid_of
    use ICE_Config_Types_m,  only: obj_condensed, ini_condensed, obj_sim_param, obj_time_scheme
    use ICE_Load_Table,      only: Load_Table
    use ICE_Lib_Evaporation, only: Ru, iMv, iLv, icpv, iLe, iYinf, iLvMvOverRu, iinvTboil, ialphaE
    implicit none
    character(len=64),   allocatable :: names(:)
    character(len=1024), allocatable :: rest(:)
    integer,             allocatable :: groups(:)
    integer :: m, p, q
    logical :: coupled

    call read_phase_file(names, groups, rest)
    nmat = max(size(names), 1)
    if (allocated(obj_condensed)) deallocate(obj_condensed)
    allocate(obj_condensed(nmat))
    do m = 1, nmat
      obj_condensed(m) = ini_condensed
      if (m <= size(names)) obj_condensed(m)%name = names(m)
    enddo

    ! Families map material-major, as ATLAS orders the populations
    if (allocated(mat_of)) deallocate(mat_of)
    allocate(mat_of(ngroups)); mat_of = 1
    if (nmat > 1) then
      if (sum(groups) /= ngroups) call refuse('INPUT/'//trim(ICE_phase_prefix)//'phase.txt declares '// &
        trim(itoa(sum(groups)))//' populations over '//trim(itoa(nmat))//' materials, but [ICE-Family*] defines '// &
        trim(itoa(ngroups))//' families')
      q = 0
      do m = 1, nmat
        do p = 1, groups(m)
          q = q + 1
          mat_of(q) = m
        enddo
      enddo
    endif

    call read_ini_vectors()

    coupled = obj_sim_param%owcoupled .or. obj_sim_param%twcoupled
    do m = 1, nmat
      ! The [ICE-Physics] models are every material's default; its tokens override them
      obj_condensed(m)%evapWord   = obj_time_scheme%evaporation
      obj_condensed(m)%intfWord   = obj_time_scheme%interface_model
      obj_condensed(m)%evapSelect = obj_time_scheme%evapSelect
      obj_condensed(m)%intfSelect = obj_time_scheme%intfSelect
      if (m <= size(rest)) call apply_tokens(m, rest(m))
      if (obj_condensed(m)%intfSelect == 1 .and. obj_condensed(m)%evapSelect == 1) &
        call refuse(material_label(m)//': interface = LK needs a gas-side evaporation model '// &
                    '(CEM, CEM-B, ASM or TC), not d2-law')
      call finalize_solidification(m)
      if (.not. coupled) obj_condensed(m)%evapSelect = 0

      ! Vapour properties, packed as Lib_Evaporation indexes them
      associate (mat => obj_condensed(m))
        mat%ep = 0._R8
        mat%ep(iMv)         = mat%Mv
        mat%ep(iLv)         = mat%lv_al
        mat%ep(icpv)        = mat%cpv
        mat%ep(iLe)         = mat%Le
        mat%ep(iYinf)       = mat%Yinf
        mat%ep(iLvMvOverRu) = mat%lv_al*mat%Mv/Ru
        mat%ep(iinvTboil)   = 1._R8/mat%Tboil
        mat%ep(ialphaE)     = mat%alphaE
      end associate
    enddo

    ! A solidifying family carries the frozen and nucleated fractions after its closure's slots
    do p = 1, ngroups
      solid_of(p) = obj_condensed(mat_of(p))%solid
      ncond(p)    = nbase(p) + merge(2, 0, solid_of(p))
    enddo

    call Load_Table()

    ! Solidification integrates the temperature with constant properties, as IGLOO does
    do m = 1, nmat
      if (.not. obj_condensed(m)%solid) cycle
      if (obj_condensed(m)%cs_varies) &
        call refuse(material_label(m)//': solidification=on requires a constant-cp material (temperature state)')
      if (obj_condensed(m)%rho_varies) &
        call refuse(material_label(m)//': solidification=on requires a constant-density material')
    enddo

  contains

    !> The material lines of the phase file: name, families, and the rest of the line (the tokens).
    subroutine read_phase_file(names, groups, rest)
      character(len=64),   allocatable, intent(out) :: names(:)
      integer,             allocatable, intent(out) :: groups(:)
      character(len=1024), allocatable, intent(out) :: rest(:)
      character(len=1024) :: line
      character(len=64)   :: f1, f2
      character(len=512)  :: file
      integer :: unit, ios, n, k1, k2
      logical :: exists

      allocate(names(0), groups(0), rest(0))
      file = 'INPUT/'//trim(ICE_phase_prefix)//'phase.txt'
      inquire(file=trim(file), exist=exists)
      if (.not. exists) return
      open(newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read(unit, '(A)', iostat=ios) line                  ! the type word
      n = 0
      do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        call next_field(line, 1, f1, k1)
        call next_field(line, k1, f2, k2)
        read(f2, *, iostat=ios) q
        if (len_trim(f1) == 0 .or. ios /= 0 .or. q < 1) &
          call refuse(trim(file)//': "'//trim(line)//'" is not a material line, "<name> <groups> [key=value ...]"')
        n = n + 1
        names  = [character(len=64)   :: names, f1]
        groups = [groups, q]
        rest   = [character(len=1024) :: rest, line(k2:)]
      enddo
      close(unit)
    end subroutine read_phase_file


    !> Applies the key=value tokens of material m.
    subroutine apply_tokens(m, tokens)
      use ICE_Lib_Evaporation, only: assign_evaporation, assign_interface
      integer,          intent(in) :: m
      character(len=*), intent(in) :: tokens
      character(len=256) :: tok, key, value
      real(R8) :: x
      integer  :: k, e, ios

      k = 1
      do
        call next_field(tokens, k, tok, k)
        if (len_trim(tok) == 0) exit
        e = index(tok, '=')
        if (e < 2) call refuse(material_label(m)//': "'//trim(tok)//'" is not a key=value token')
        key   = tok(:e-1)
        value = tok(e+1:)
        select case (trim(key))
        case ('evaporation')
          call assign_evaporation(trim(value), obj_condensed(m)%evapSelect)
          obj_condensed(m)%evapWord = value
        case ('interface')
          call assign_interface(trim(value), obj_condensed(m)%intfSelect)
          obj_condensed(m)%intfWord = value
        case ('liquid-conduction')
          if (value /= 'ITC') call refuse_word(key, value, ['ITC'])
        case ('boiling')
          if (value /= 'clamp') call refuse_word(key, value, ['clamp'])
        case ('combustion')
          if (value /= 'none') call refuse_word(key, value, ['none'])
        case ('solidification')
          select case (trim(value))
          case ('on');  obj_condensed(m)%solidSelect = 1
          case ('off'); obj_condensed(m)%solidSelect = 0
          case default; call refuse_word(key, value, ['on ', 'off'])
          end select
        case default
          if (.not. any(real_keys == key)) &
            call refuse(material_label(m)//': unknown key "'//trim(key)//'"; the keys are '//key_list())
          read(value, *, iostat=ios) x
          if (ios /= 0 .or. index(value, ',') > 0) &
            call refuse(material_label(m)//': '//trim(key)//'='//trim(value)//' is not a real number')
          select case (trim(key))
          case ('alpha-e')
            if (.not. x > 0._R8) call refuse(material_label(m)//': alpha-e='//trim(value)//' must be > 0')
            obj_condensed(m)%alphaE = x
          case ('T-melt');   obj_condensed(m)%Tmelt = x
          case ('h-fus');    obj_condensed(m)%hFus  = x
          case ('T-nuc');    obj_condensed(m)%Tnuc  = x
          case ('cp-solid'); obj_condensed(m)%cpSol = x
          end select
        end select
      enddo
    end subroutine apply_tokens


    !> The resolved solidification inputs of material m and IGLOO's refusals, before an uncoupled run drops evaporation.
    subroutine finalize_solidification(m)
      integer, intent(in) :: m
      if (obj_condensed(m)%Tnuc <= 0._R8) obj_condensed(m)%Tnuc = 0.8_R8*obj_condensed(m)%Tmelt
      obj_condensed(m)%solid = obj_condensed(m)%solidSelect > 0
      if (.not. obj_condensed(m)%solid) return
      if (.not. obj_condensed(m)%hFus > 0._R8) &
        call refuse(material_label(m)//': solidification=on requires h-fus > 0 in [GPB-PhaseX]')
      if (.not. obj_condensed(m)%cpSol > 0._R8) &
        call refuse(material_label(m)//': solidification=on requires cp-solid > 0 in [GPB-PhaseX]')
      if (obj_condensed(m)%Tnuc >= obj_condensed(m)%Tmelt) &
        call refuse(material_label(m)//': solidification=on requires T-nuc < T-melt')
      if (obj_condensed(m)%evapSelect > 0) &
        call refuse(material_label(m)//': solidification=on is exclusive with evaporation (alumina does not evaporate)')
    end subroutine finalize_solidification

  end subroutine Setup_Materials


  !> With several materials, a [ICE-Physics] property key carries one value per material; with
  !  one, the registry has already read it. The count is checked either way. A key may be given
  !  under its alias instead, never under both names.
  subroutine read_ini_vectors()
    use Finer,              only: file_ini
    use ICE_Backend_INI,    only: Open_Ini
    use ICE_Parameters_m,   only: codename
    use ICE_Global_m,       only: nmat
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    type(file_ini)      :: fini
    character(len=1024) :: buf, abuf
    character(len=256)  :: tok
    character(len=23)   :: key
    real(R8) :: x
    integer  :: i, m, n, k, err, aerr, ios

    call Open_Ini(fini)
    do i = 1, size(vector_keys)
      key = vector_keys(i)
      buf = ''
      call fini%get(section_name=trim(codename)//'-Physics', option_name=trim(key), val=buf, error=err)
      if (len_trim(vector_alias(i)) > 0) then
        abuf = ''
        call fini%get(section_name=trim(codename)//'-Physics', option_name=trim(vector_alias(i)), val=abuf, error=aerr)
        if (aerr == 0) then
          if (err == 0) call refuse('['//trim(codename)//'-Physics]: give '//trim(key)//' or its alias '// &
                                    trim(vector_alias(i))//', not both')
          key = vector_alias(i); buf = abuf; err = 0
        endif
      endif
      if (err /= 0 .or. len_trim(buf) == 0) cycle
      n = count_fields(buf)
      if (n /= nmat) call refuse('['//trim(codename)//'-Physics] '//trim(key)//' carries '// &
        trim(itoa(n))//' values for '//trim(itoa(nmat))//' material(s): give one value per material')
      if (nmat == 1) cycle
      k = 1
      do m = 1, nmat
        call next_field(buf, k, tok, k)
        read(tok, *, iostat=ios) x
        if (ios /= 0 .or. index(tok, ',') > 0) call refuse('['//trim(codename)//'-Physics] '// &
          trim(key)//': "'//trim(tok)//'" is not a real number')
        if ((vector_positive(i) .and. .not. x > 0._R8) .or. (.not. vector_positive(i) .and. .not. x >= 0._R8)) &
          call refuse('['//trim(codename)//'-Physics] '//trim(key)//' = '//trim(tok)// &
                      merge(' must be > 0 ', ' must be >= 0', vector_positive(i)))
        select case (i)
        case (1);  obj_condensed(m)%rho_al = x
        case (2);  obj_condensed(m)%cs_al  = x
        case (3);  obj_condensed(m)%lv_al  = x
        case (4);  obj_condensed(m)%emiss  = x
        case (5);  obj_condensed(m)%Mv     = x
        case (6);  obj_condensed(m)%Tboil  = x
        case (7);  obj_condensed(m)%cpv    = x
        case (8);  obj_condensed(m)%Le     = x
        case (9);  obj_condensed(m)%Yinf   = x
        case (10); obj_condensed(m)%alphaE = x
        end select
      enddo
    enddo
    call fini%free()
  end subroutine read_ini_vectors


  !> The next blank-separated field of line from position k; k returns just after it.
  subroutine next_field(line, k, field, knext)
    character(len=*), intent(in)  :: line
    integer,          intent(in)  :: k
    character(len=*), intent(out) :: field
    integer,          intent(out) :: knext
    integer :: a, b
    field = ''
    a = k
    do while (a <= len(line))
      if (line(a:a) /= ' ' .and. line(a:a) /= char(9)) exit
      a = a + 1
    enddo
    b = a
    do while (b <= len(line))
      if (line(b:b) == ' ' .or. line(b:b) == char(9)) exit
      b = b + 1
    enddo
    if (a <= len(line)) field = line(a:b-1)
    knext = b
  end subroutine next_field


  integer function count_fields(line)
    character(len=*), intent(in) :: line
    character(len=256) :: f
    integer :: k
    count_fields = 0
    k = 1
    do
      call next_field(line, k, f, k)
      if (len_trim(f) == 0) exit
      count_fields = count_fields + 1
    enddo
  end function count_fields


  function material_label(m) result(txt)
    use ICE_Global_m,       only: ICE_phase_prefix
    use ICE_Config_Types_m, only: obj_condensed
    integer, intent(in) :: m
    character(len=:), allocatable :: txt
    txt = 'INPUT/'//trim(ICE_phase_prefix)//'phase.txt, material '//trim(itoa(m))//' ('// &
          trim(obj_condensed(m)%name)//')'
  end function material_label


  function key_list() result(txt)
    character(len=:), allocatable :: txt
    integer :: i
    txt = trim(word_keys(1))
    do i = 2, size(word_keys)
      txt = txt//', '//trim(word_keys(i))
    enddo
    do i = 1, size(real_keys)
      txt = txt//', '//trim(real_keys(i))
    enddo
  end function key_list


  !> A model value ICE does not implement or recognise, in the format of the INI model keys.
  subroutine refuse_word(key, value, menu)
    character(len=*), intent(in) :: key, value, menu(:)
    integer :: i
    write(*,*)
    write(*,'(A)') ' Wrong '//trim(key)//' input ---> '//trim(value)
    write(*,'(A)') ' Choose one of the following :'
    do i = 1, size(menu)
      write(*,'(A)') ' - '//trim(menu(i))
    enddo
    write(*,*)
    error stop 1
  end subroutine refuse_word


  !> Every rank prints, as the BC reader does.
  subroutine refuse(msg)
    character(len=*), intent(in) :: msg
    write(*,'(A)') ' [ERROR] [ICE::Setup_Materials] '//msg
    error stop 1
  end subroutine refuse


  pure function itoa(i) result(txt)
    integer, intent(in) :: i
    character(len=16) :: txt
    write(txt, '(I0)') i
  end function itoa

end module ICE_Setup_Materials
