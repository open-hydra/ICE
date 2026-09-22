module ICE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Global_m,         only: ICE_phase_prefix
  use ICE_Advanced_Types_m, only: ICE_domain_type, ICE_bc_type

  implicit none
  private
  public :: Setup_BC, Print_BC_Summary

  integer :: nconn = 0, nsym = 0, nio = 0, next = 0, nchim = 0

contains


  !> Read INPUT/<prefix>bc.txt.
  !
  !  TWO schemas are accepted, distinguished by the column count of the first record:
  !
  !    ATLAS "401" schema (6 columns)   -- what ATLAS BCB emits today. Header is
  !        b  i  j  k  f  code
  !      with a three-digit code, and ONE record per boundary CELL (not per cell-and-group).
  !      This is the schema MOSE (MOSE_IO_BC) and IGLOO (read_cdp_bc_file) already read.
  !
  !    Legacy schema (8 columns)        -- header b i j k f mat p type with single-digit type,
  !      one record per (cell, group). Retained ONLY because ICE's own regression data
  !      (test/Doisneau/*/INPUT/part-bc.txt) is still in it. Regenerate those with ATLAS BCB
  !      and this branch can go.
  !
  !  The ATLAS codes are TRANSLATED to ICE's internal bc%type values on read, so every
  !  downstream `select case (bc%type)` (Lib_Ghost, Mod_BC_Fluxes, Mod_GhostExchange) is
  !  untouched by this migration. The mapping was derived by regenerating one case both ways
  !  and matching the per-code census exactly (15660 faces either way):
  !
  !      ATLAS 200 (axisymmetry)  -> 2   ICE skips these
  !      ATLAS 300 (symmetry)     -> 3
  !      ATLAS 400 (extrapolation)-> 11
  !      ATLAS 401 (inlet)        -> 4   with injtype = 0
  !      ATLAS 402/403 (inlet)    -> 4   with injtype = 1/2
  !      ATLAS 101/201 (connect)  -> 1
  !      ATLAS 102 (chimera)      -> 102 no legacy equivalent, so the ATLAS code is kept;
  !                                      filled by Ghost_Chimera in Lib_Ghost
  subroutine Setup_BC(grid)
    use ICE_Global_m, only: ngroups
    implicit none
    type(ICE_domain_type), intent(inout), target :: grid
    type(ICE_bc_type), dimension(:), pointer     :: bc
    integer(I4), dimension(:,:), pointer         :: nbc
    integer(I4)         :: unitfile, ios, ntok
    character(len=512)  :: firstline
    logical             :: legacy

    bc  => grid%bc
    nbc => grid%n_bf
    nconn = 0; nsym = 0; nio = 0; next = 0; nchim = 0

    open(newunit=unitfile, file='INPUT/'//trim(ICE_phase_prefix)//'bc.txt', &
         status='old', iostat=ios)
    if (ios /= 0) error stop 'BC file not found'

    !> Schema probe: count blank-separated tokens on the first record.
    read(unitfile,'(A)',iostat=ios) firstline
    if (ios /= 0) error stop 'BC file is empty'
    ntok = Count_Tokens(firstline)
    rewind(unitfile)

    select case (ntok)
      case (6); legacy = .false.
      case (8); legacy = .true.
      case default
        write(*,'(A,I0,A)') '  [ICE::Setup_BC] first BC record has ', ntok, &
                            ' columns; expected 6 (ATLAS) or 8 (legacy).'
        error stop
    end select

    nbc = 0
    if (legacy) then
      call Read_BC_Legacy(bc, nbc, unitfile)
    else
      call Read_BC_ATLAS(bc, nbc, unitfile, ngroups)
    endif

    close(unitfile)
    call Print_BC_Summary()

  end subroutine Setup_BC


  !> Blank-separated token count, used only for the schema probe.
  pure function Count_Tokens(line) result(n)
    implicit none
    character(len=*), intent(in) :: line
    integer(I4) :: n, i
    logical     :: inside

    n = 0; inside = .false.
    do i = 1, len_trim(line)
      if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
        inside = .false.
      else
        if (.not. inside) n = n + 1
        inside = .true.
      endif
    enddo

  end function Count_Tokens


  !> ATLAS "401" schema. ONE record per boundary cell; ICE's bc array holds one entry per
  !  (cell, group), so each record is fanned out over the groups with bc%p set accordingly.
  subroutine Read_BC_ATLAS(bc, nbc, unitfile, ngroups)
    implicit none
    type(ICE_bc_type), dimension(:), intent(inout) :: bc
    integer(I4), dimension(:,:),     intent(inout) :: nbc
    integer(I4),                     intent(in)    :: unitfile, ngroups
    integer(I4) :: n, ncell, p, idx, ios, code
    integer(I4) :: hb, hi, hj, hk, hf
    integer(I4) :: cs(9)
    real(R8)    :: massflux, velocity, temperature, radius, alpha, beta
    integer(I4) :: injtype
    character(len=32) :: alpha_tok, beta_tok
    integer(I4) :: nchi(2), s, nzero_chim
    integer(I4),  allocatable :: donorID(:,:)
    real(R8),     allocatable :: weight(:)

    if (ngroups <= 0) error stop 'Read_BC_ATLAS: ngroups not set'
    if (mod(size(bc), ngroups) /= 0) error stop 'Read_BC_ATLAS: bc array is not a multiple of ngroups'
    ncell = size(bc) / ngroups
    nzero_chim = 0

    do n = 1, ncell

      read(unitfile,*,iostat=ios) hb, hi, hj, hk, hf, code
      if (ios /= 0) then
        write(*,'(A,I0)') '  [ICE::Setup_BC] error reading BC header record ', n
        error stop
      endif

      !> Same gate IGLOO uses: a legacy record, or a property line swallowed as a header,
      !  lands outside {0, 100..999}. Lenient list-directed reads would otherwise mistype
      !  every face and corrupt the setup silently.
      if (code /= 0 .and. (code < 100 .or. code > 999)) then
        write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] bad BC code ', code, ' at record ', n
        write(*,'(A)') '  expected 6 columns with a 3-digit code (or 0) in column 6.'
        error stop
      endif

      alpha = 0._R8; beta = 0._R8
      massflux = 0._R8; velocity = 0._R8; temperature = 0._R8; radius = 0._R8
      injtype = 0
      cs = 0

      select case (code)

        case (101, 201)          ! connection / periodic -> ICE type 1
          nconn = nconn + 1
          read(unitfile,*,iostat=ios) cs(1:9)

        case (102)               ! chimera -> ICE type 102
          !> Payload: `nchi_g1 nchi_g2`, then nchi_g1 + nchi_g2 donor lines `b i j k weight`,
          !  ghost layer 1 first. Weights are normalised per layer by ATLAS.
          nchim = nchim + 1
          read(unitfile,*,iostat=ios) nchi
          if (ios == 0) then
            if (allocated(donorID)) deallocate(donorID, weight)
            allocate(donorID(sum(nchi),4), weight(sum(nchi)))
            do s = 1, sum(nchi)
              read(unitfile,*,iostat=ios) donorID(s,1:4), weight(s)
              if (ios /= 0) exit
            enddo
            if (ios == 0) then
              if (sum(weight(1:nchi(1))) < 0.5_R8 .or. &
                  sum(weight(nchi(1)+1:sum(nchi))) < 0.5_R8) nzero_chim = nzero_chim + 1
            endif
          endif

        case (200)               ! axisymmetry -> ICE type 2 (skipped downstream), no payload
          continue

        case (300)               ! symmetry -> ICE type 3, no payload
          nsym = nsym + 1

        case (400)               ! extrapolation -> ICE type 11, no payload
          next = next + 1

        case (401)               ! inlet, krho-style loading -> ICE type 4, injtype 0
          nio = nio + 1
          !> Condensed payload, 9 fields; ICE consumes the first six:
          !    krho, velocity-ratio, alpha, beta, temperature-ratio, radius,
          !    [ distribution-width, distribution-law, per-cell ds ]  <- IGLOO-only, ignored here.
          !  Field-for-field identical to the legacy 7-field payload minus its leading injtype,
          !  verified by regenerating the same case in both schemas.
          read(unitfile,*,iostat=ios) massflux, velocity, alpha_tok, beta_tok, temperature, radius
          alpha = Parse_Dir_Token(alpha_tok)
          beta  = Parse_Dir_Token(beta_tok)
          injtype = 0        ! krho loading: velocity and temperature are RATIOS to the gas

        case (402, 403)
          !> Prescribed-mass-flux inlets. Column meanings differ from 401 and from each other;
          !  the contract is stated at utils/ATLAS/src/BCB/builder_400dp.f90:24-56 and restated
          !  in IGLOO's obj_particles.f90:214-221, whose initializePart applies it:
          !
          !    col 1  gp  [kg/(s m^2)] mass flux          (both)   -- NOT krho
          !    col 2  402: |v| absolute [m/s]   403: kV velocity scaling on the local gas speed
          !    col 5  402: Tp absolute [K]      403: Tp absolute [K]
          !
          !  That lines up one-to-one with ICE's injtypes (Lib_Ghost.f90:123-160):
          !    injtype 1  massflux, velocity and temperature all absolute   <- 402
          !    injtype 2  massflux and temperature absolute, velocity a ratio of |v_gas|  <- 403
          !  In both, ICE forms the density as massflux/|v.n|, and gp is already kg/(s m^2),
          !  so the units match with no conversion.
          nio = nio + 1
          read(unitfile,*,iostat=ios) massflux, velocity, alpha_tok, beta_tok, temperature, radius
          alpha = Parse_Dir_Token(alpha_tok)
          beta  = Parse_Dir_Token(beta_tok)
          if (code == 402) then; injtype = 1; else; injtype = 2; endif

        case (0)
          continue

        case default
          write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] unsupported ATLAS BC code ', code, &
                                 ' at record ', n
          error stop

      end select

      if (ios /= 0) then
        write(*,'(A,I0)') '  [ICE::Setup_BC] error reading payload for record ', n
        error stop
      endif

      !> One boundary FACE per record -- counted once, not once per group.
      nbc(hb, hf) = nbc(hb, hf) + 1

      do p = 1, ngroups
        idx = (n-1)*ngroups + p
        bc(idx)%b = hb; bc(idx)%i = hi; bc(idx)%j = hj; bc(idx)%k = hk; bc(idx)%f = hf
        bc(idx)%p    = p
        bc(idx)%mat  = p          ! legacy field; unused outside this reader
        bc(idx)%type = ICE_BC_Type_From_ATLAS(code)
        bc(idx)%injtype    = injtype
        bc(idx)%massflux   = massflux
        bc(idx)%velocity   = velocity
        bc(idx)%alpha      = alpha
        bc(idx)%beta       = beta
        bc(idx)%temperature = temperature
        bc(idx)%radius     = radius
        if (code == 101 .or. code == 201) then
          bc(idx)%bs = cs(1); bc(idx)%is = cs(2); bc(idx)%js = cs(3)
          bc(idx)%ks = cs(4); bc(idx)%fs = cs(5)
          bc(idx)%d11 = cs(6); bc(idx)%d12 = cs(7)
          bc(idx)%d21 = cs(8); bc(idx)%d22 = cs(9)
        endif
        if (code == 102) then
          bc(idx)%ni              = nchi
          bc(idx)%donorID         = donorID
          bc(idx)%volume_fraction = weight
        endif
      enddo

    end do

    if (nzero_chim > 0) then
      write(*,'(A,I0,A)') '  [ICE::Setup_BC] ', nzero_chim, ' chimera records have a ghost layer '// &
        'whose donor weights sum to ~0 (no donor found by ATLAS BCB). Fix the BC file before running.'
      error stop
    endif

  end subroutine Read_BC_ATLAS


  !> ATLAS code -> ICE internal bc%type. See the table in Setup_BC's header.
  pure function ICE_BC_Type_From_ATLAS(code) result(t)
    implicit none
    integer(I4), intent(in) :: code
    integer(I4) :: t

    select case (code)
      case (101, 201); t = 1
      case (200);      t = 2
      case (300);      t = 3
      case (400);      t = 11
      case (102);      t = 102
      case (401:403);  t = 4
      case default;    t = 0
    end select

  end function ICE_BC_Type_From_ATLAS


  !> Direction token. ATLAS writes the literal `normal` when the injection direction is the
  !  face normal. ICE's convention for that is alpha = beta = 0, which Lib_Ghost.f90:116-120
  !  detects and replaces with the actual face-normal angles -- so `normal` maps to 0, NOT to
  !  MOSE's huge() sentinel. (Confirmed against the legacy file for the same case, which wrote
  !  0.0 in both slots.)
  pure function Parse_Dir_Token(tok) result(val)
    implicit none
    character(len=*), intent(in) :: tok
    real(R8) :: val
    integer  :: ios_loc
    character(len=len(tok)) :: t

    t = adjustl(tok)
    if (index(trim(t), 'normal') > 0) then
      val = 0._R8
    else
      read(t, *, iostat=ios_loc) val
      if (ios_loc /= 0) val = 0._R8
    endif

  end function Parse_Dir_Token


  !> Legacy 8-column schema, verbatim from before the ATLAS migration. Kept only for ICE's
  !  own regression data; see Setup_BC's header.
  subroutine Read_BC_Legacy(bc, nbc, unitfile)
    implicit none
    type(ICE_bc_type), dimension(:), intent(inout) :: bc
    integer(I4), dimension(:,:),     intent(inout) :: nbc
    integer(I4),                     intent(in)    :: unitfile
    integer(I4) :: i, ios
    real(R8)    :: injtype_r

    do i = 1, size(bc)

      read(unitfile,*,iostat=ios) &
        bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%mat, bc(i)%p, bc(i)%type
      if (ios /= 0) write(*,*) '  Error reading BC file at entry', i

      nbc(bc(i)%b, bc(i)%f) = nbc(bc(i)%b, bc(i)%f) + 1

      select case (bc(i)%type)
        case (1)
          nconn = nconn + 1
          read(unitfile,*,iostat=ios) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
            bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
        case (3)
          nsym = nsym + 1
        case (4, 14)
          nio = nio + 1
          read(unitfile,*,iostat=ios) &
            injtype_r, bc(i)%massflux, bc(i)%velocity, &
            bc(i)%alpha, bc(i)%beta, bc(i)%temperature, bc(i)%radius
          bc(i)%injtype = nint(injtype_r)
        case (11)
          next = next + 1
      end select

      if (ios /= 0) then
        write(*,*) '  Error in BC file'
        error stop
      end if

    end do

  end subroutine Read_BC_Legacy


  subroutine Print_BC_Summary()
    implicit none
    write(*,*)
    write(*,'(A)') ' Boundary Conditions:'
    if (nconn > 0) write(*,'(A,T35,I0)') '   Connection',    nconn
    if (nsym  > 0) write(*,'(A,T35,I0)') '   Symmetry',      nsym
    if (nio   > 0) write(*,'(A,T35,I0)') '   Inflow',        nio
    if (next  > 0) write(*,'(A,T35,I0)') '   Extrapolation', next
    if (nchim > 0) write(*,'(A,T35,I0)') '   Chimera',       nchim
    write(*,*)
  end subroutine Print_BC_Summary


end module ICE_IO_BC
