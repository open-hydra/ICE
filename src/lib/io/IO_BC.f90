module ICE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Global_m,         only: ICE_phase_prefix
  use ICE_Advanced_Types_m, only: ICE_domain_type, ICE_bc_type

  implicit none
  private
  public :: Setup_BC, Print_BC_Summary

  integer :: nconn = 0, nsym = 0, nio = 0, next = 0, nchim = 0, naxi = 0

  !> One record of the bc file: the header and its payload
  type :: bc_record
    integer(I4) :: b = 0, i = 0, j = 0, k = 0, f = 0, code = 0
    real(R8)    :: massflux = 0._R8, velocity = 0._R8, temperature = 0._R8, radius = 0._R8
    real(R8)    :: alpha = 0._R8, beta = 0._R8
    integer(I4) :: cs(9) = 0, nchi(2) = 0
    integer(I4), allocatable :: donorID(:,:)
    real(R8),    allocatable :: weight(:)
  end type bc_record

contains


  !> Read INPUT/<prefix>bc.txt, written by ATLAS BCB. Header per boundary cell:
  !
  !      b  i  j  k  f  code
  !
  !  followed by a code-dependent payload. ICE's bc array holds one entry per (cell, group).
  !  The file holds either ONE record per boundary cell, fanned out over the groups, or one
  !  block of records per group inside each mesh block, in ATLAS order (mesh block, group,
  !  faces), each group repeating the first one's faces. bc%type stores the ATLAS code as is:
  !
  !      0        planar 2-D face        ghost = copy of the interior cell, no boundary flux
  !      101/201  connection / periodic
  !      102      chimera
  !      200      wedge side face        ghosts mirrored; boundary flux for IG/AG (the hoop
  !                                      pressure), none for MK
  !      300      symmetry
  !      301      dispersed-phase wall   treated as 300: the symmetry ghost already absorbs
  !                                      particles moving toward the face and mirrors receding ones
  !      400      extrapolation
  !      401-403  inlet (see Read payload below)
  subroutine Setup_BC(grid)
    use ICE_Global_m, only: ngroups
    use ICE_Mod_MPI, only: mpi_is_root
    implicit none
    type(ICE_domain_type), intent(inout), target :: grid
    type(ICE_bc_type), dimension(:), pointer     :: bc
    integer(I4), dimension(:,:), pointer         :: nbc
    type(bc_record), allocatable :: recs(:)
    integer(I4) :: unitfile, ios
    integer(I4) :: n, ncell, nrec, p, code, pos, nb, c, r, n0
    integer(I4) :: hb, hi, hj, hk, hf, nzero_chim

    bc  => grid%bc
    nbc => grid%n_bf
    nconn = 0; nsym = 0; nio = 0; next = 0; nchim = 0; naxi = 0

    if (ngroups <= 0) error stop 'Setup_BC: ngroups not set'
    if (mod(size(bc), ngroups) /= 0) error stop 'Setup_BC: bc array is not a multiple of ngroups'
    ncell = size(bc) / ngroups
    nzero_chim = 0

    open(newunit=unitfile, file='INPUT/'//trim(ICE_phase_prefix)//'bc.txt', &
         status='old', iostat=ios)
    if (ios /= 0) error stop 'BC file not found'

    ! Every record of the file
    allocate(recs(ngroups*ncell))
    nrec = 0
    do
      read(unitfile,*,iostat=ios) hb, hi, hj, hk, hf, code
      if (ios < 0) exit
      if (ios /= 0) then
        write(*,'(A,I0)') '  [ICE::Setup_BC] error reading BC header record ', nrec + 1
        error stop
      endif
      if (nrec == size(recs)) then
        write(*,'(A,I0,A,I0,A)') '  [ICE::Setup_BC] the BC file holds more than ', size(recs), &
          ' records (', ngroups, ' group(s) times the boundary faces)'
        error stop
      endif
      nrec = nrec + 1

      !> Same gate IGLOO uses: a pre-ATLAS 8-column record, or a property line swallowed as
      !  a header, lands outside {0, 100..999}. Lenient list-directed reads would otherwise
      !  mistype every face and corrupt the setup silently.
      if (code /= 0 .and. (code < 100 .or. code > 999)) then
        write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] bad BC code ', code, ' at record ', nrec
        write(*,'(A)') '  expected the ATLAS schema: 6 columns with a 3-digit code (or 0) in column 6.'
        error stop
      endif
      recs(nrec)%b = hb; recs(nrec)%i = hi; recs(nrec)%j = hj; recs(nrec)%k = hk; recs(nrec)%f = hf
      recs(nrec)%code = code
      call read_payload(unitfile, recs(nrec), nrec, nzero_chim)
    end do
    close(unitfile)

    nbc = 0
    if (nrec == ncell) then
      !> One record per boundary face, fanned out over the groups
      do n = 1, ncell
        call count_face(recs(n), nbc)
        do p = 1, ngroups
          call store_record(bc((n-1)*ngroups + p), recs(n), p)
        enddo
      enddo
    elseif (ngroups > 1 .and. nrec == ngroups*ncell) then
      !> One block per group inside each mesh block
      pos = 1; n0 = 0
      do while (pos <= nrec)
        hb = recs(pos)%b
        nb = 2*(grid%blk(hb)%dim(2)*grid%blk(hb)%dim(3) + grid%blk(hb)%dim(1)*grid%blk(hb)%dim(3) + &
                grid%blk(hb)%dim(1)*grid%blk(hb)%dim(2))
        if (pos + ngroups*nb - 1 > nrec) then
          write(*,'(A,I0,A,I0,A)') '  [ICE::Setup_BC] mesh block ', hb, ' needs ', ngroups*nb, &
            ' records (one block per group) but the file ends before'
          error stop
        endif
        do c = 1, nb
          call count_face(recs(pos + c - 1), nbc)
        enddo
        do p = 1, ngroups
          do c = 1, nb
            r = pos + (p-1)*nb + c - 1
            if (recs(r)%b /= recs(pos+c-1)%b .or. recs(r)%i /= recs(pos+c-1)%i .or. recs(r)%j /= recs(pos+c-1)%j &
                .or. recs(r)%k /= recs(pos+c-1)%k .or. recs(r)%f /= recs(pos+c-1)%f) then
              write(*,'(A,I0,A,I0,A,I0)') '  [ICE::Setup_BC] group ', p, ' of mesh block ', hb, &
                ' does not repeat the faces of group 1, at record ', r
              error stop
            endif
            call store_record(bc((n0 + c - 1)*ngroups + p), recs(r), p)
          enddo
        enddo
        n0 = n0 + nb
        pos = pos + ngroups*nb
      enddo
    else
      write(*,'(A,I0,A,I0,A,I0,A)') '  [ICE::Setup_BC] the BC file holds ', nrec, ' records: expected ', ncell, &
        ' (one per boundary face) or ', ngroups*ncell, ' (one block per group)'
      error stop
    endif

    if (nzero_chim > 0) then
      write(*,'(A,I0,A)') '  [ICE::Setup_BC] ', nzero_chim, ' chimera records have a ghost layer '// &
        'whose donor weights sum to ~0 (no donor found by ATLAS BCB). Fix the BC file before running.'
      error stop
    endif

    if (mpi_is_root) then
      call Print_BC_Summary()
      call Check_Phase_Groups(ngroups)
    endif

  end subroutine Setup_BC


  !> The payload of one record, by its code.
  subroutine read_payload(unitfile, rec, n, nzero_chim)
    implicit none
    integer(I4),     intent(in)    :: unitfile, n
    type(bc_record), intent(inout) :: rec
    integer(I4),     intent(inout) :: nzero_chim
    character(len=32) :: alpha_tok, beta_tok
    integer(I4) :: ios, s

    ios = 0
    select case (rec%code)

      case (101, 201)          ! connection / periodic
        read(unitfile,*,iostat=ios) rec%cs(1:9)

      case (102)               ! chimera
        !> Payload: `nchi_g1 nchi_g2`, then nchi_g1 + nchi_g2 donor lines `b i j k weight`,
        !  ghost layer 1 first. Weights are normalised per layer by ATLAS.
        read(unitfile,*,iostat=ios) rec%nchi
        if (ios == 0) then
          allocate(rec%donorID(sum(rec%nchi),4), rec%weight(sum(rec%nchi)))
          do s = 1, sum(rec%nchi)
            read(unitfile,*,iostat=ios) rec%donorID(s,1:4), rec%weight(s)
            if (ios /= 0) exit
          enddo
          if (ios == 0) then
            if (sum(rec%weight(1:rec%nchi(1))) < 0.5_R8 .or. &
                sum(rec%weight(rec%nchi(1)+1:sum(rec%nchi))) < 0.5_R8) nzero_chim = nzero_chim + 1
          endif
        endif

      case (200, 300, 301, 400, 0)   ! no payload
        continue

      case (401:403)           ! inlet
        !> Payload, 9 fields; ICE consumes the first six. Column meanings per code, as in
        !  ATLAS builder_400dp.f90 and IGLOO's obj_particles.f90:
        !
        !    401  krho, kV, alpha, beta, kT, radius     loading, velocity and temperature
        !                                               are RATIOS to the local gas
        !    402  gp,   |v|, alpha, beta, Tp, radius    mass flux [kg/(s m^2)], speed and
        !                                               temperature all absolute
        !    403  gp,   kV,  alpha, beta, Tp, radius    as 402, but speed is a ratio of |v_gas|
        !
        !  [ sigmap, distribution law, per-cell ds ] follow and are IGLOO-only.
        read(unitfile,*,iostat=ios) rec%massflux, rec%velocity, alpha_tok, beta_tok, rec%temperature, rec%radius
        rec%alpha = Parse_Dir_Token(alpha_tok)
        rec%beta  = Parse_Dir_Token(beta_tok)

      case default
        write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] unsupported ATLAS BC code ', rec%code, &
                               ' at record ', n
        error stop

    end select

    if (ios /= 0) then
      write(*,'(A,I0)') '  [ICE::Setup_BC] error reading payload for record ', n
      error stop
    endif

  end subroutine read_payload


  !> One boundary FACE: counted once, not once per group.
  subroutine count_face(rec, nbc)
    implicit none
    type(bc_record),  intent(in)    :: rec
    integer(I4),      intent(inout) :: nbc(:,:)
    nbc(rec%b, rec%f) = nbc(rec%b, rec%f) + 1
    select case (rec%code)
      case (101, 201);   nconn = nconn + 1
      case (102);        nchim = nchim + 1
      case (200);        naxi  = naxi + 1
      case (300, 301);   nsym  = nsym + 1
      case (400);        next  = next + 1
      case (401:403);    nio   = nio + 1
    end select
  end subroutine count_face


  subroutine store_record(bcp, rec, p)
    implicit none
    type(ICE_bc_type), intent(inout) :: bcp
    type(bc_record),   intent(in)    :: rec
    integer(I4),       intent(in)    :: p
    bcp%b = rec%b; bcp%i = rec%i; bcp%j = rec%j; bcp%k = rec%k; bcp%f = rec%f
    bcp%p    = p
    bcp%type = rec%code
    bcp%massflux    = rec%massflux
    bcp%velocity    = rec%velocity
    bcp%alpha       = rec%alpha
    bcp%beta        = rec%beta
    bcp%temperature = rec%temperature
    bcp%radius      = rec%radius
    if (rec%code == 101 .or. rec%code == 201) then
      bcp%bs = rec%cs(1); bcp%is = rec%cs(2); bcp%js = rec%cs(3)
      bcp%ks = rec%cs(4); bcp%fs = rec%cs(5)
      bcp%d11 = rec%cs(6); bcp%d12 = rec%cs(7)
      bcp%d21 = rec%cs(8); bcp%d22 = rec%cs(9)
    endif
    if (rec%code == 102) then
      bcp%ni              = rec%nchi
      bcp%donorID         = rec%donorID
      bcp%volume_fraction = rec%weight
    endif
  end subroutine store_record


  !> ATLAS writes INPUT/<prefix>phase.txt for the condensed phase ICE solves: one line per material,
  !  "<name> <groups> [key=value ...]" (the tokens are IGLOO's; the list-directed read below stops after
  !  the first two items). ICE takes its group count from [ICE-Family*] and never opened this file, so a
  !  mismatch went unnoticed -- and the bc.txt carries one block of records per (material, population),
  !  of which Setup_BC consumes the first and fans it out over ngroups. Warn when the two disagree.
  !  The file is optional (standalone ICE cases may not carry it): absent = silent.
  subroutine Check_Phase_Groups(ngroups)
    implicit none
    integer(I4), intent(in) :: ngroups
    integer(I4)         :: u, ios, ng, nsum, nmat
    character(len=512)  :: line
    character(len=64)   :: name

    open(newunit=u, file='INPUT/'//trim(ICE_phase_prefix)//'phase.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    read(u, '(A)', iostat=ios) line          ! header: "<type word> [modeling=...]"
    nsum = 0; nmat = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios) name, ng
      if (ios /= 0) cycle
      nmat = nmat + 1; nsum = nsum + ng
    enddo
    close(u)
    if (nmat == 0) return
    if (nsum /= ngroups) then
      write(*,'(A,I0,A,I0,A)') '  [WARNING] INPUT/'//trim(ICE_phase_prefix)//'phase.txt declares ', nsum, &
        ' population(s) over ', nmat, ' material(s), but [ICE-Family*] defines '
      write(*,'(A,I0,A)') '            ', ngroups, ' group(s). The bc.txt carries one record block per (material,'// &
        ' population); ICE reads the FIRST block and fans it out over its groups, so the two counts should agree.'
    endif

  end subroutine Check_Phase_Groups


  !> Direction token. ATLAS writes the literal `normal` when the injection direction is the
  !  face normal. ICE's convention for that is alpha = beta = 0, which Lib_Ghost detects and
  !  replaces with the actual face-normal angles -- so `normal` maps to 0, NOT to MOSE's
  !  huge() sentinel.
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


  subroutine Print_BC_Summary()
    implicit none
    write(*,*)
    write(*,'(A)') ' Boundary Conditions:'
    if (nconn > 0) write(*,'(A,T35,I0)') '   Connection',    nconn
    if (nsym  > 0) write(*,'(A,T35,I0)') '   Symmetry',      nsym
    if (naxi  > 0) write(*,'(A,T35,I0)') '   Axisymmetry',   naxi
    if (nio   > 0) write(*,'(A,T35,I0)') '   Inflow',        nio
    if (next  > 0) write(*,'(A,T35,I0)') '   Extrapolation', next
    if (nchim > 0) write(*,'(A,T35,I0)') '   Chimera',       nchim
    write(*,*)
  end subroutine Print_BC_Summary


end module ICE_IO_BC
