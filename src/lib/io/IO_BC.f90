module ICE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Global_m,         only: ICE_phase_prefix
  use ICE_Advanced_Types_m, only: ICE_domain_type, ICE_bc_type
  implicit none
  private

  public :: Setup_BC, Print_BC_Summary

  integer :: nconn = 0, nsym = 0, nio = 0, next = 0, nchim = 0, naxi = 0

contains

  !> Read the BC file of one grid level, written by ATLAS BCB: INPUT/<prefix>bc.txt
  !  for the fine level and INPUT/<prefix>bc<level>.txt for a coarse multigrid level.
  !  Header per boundary cell:
  !
  !      b  i  j  k  f  code
  !
  !> Header per boundary cell:
  !
  !      b  i  j  k  f  code
  !
  !  followed by a code-dependent payload. ICE's bc array holds one entry per
  !  (cell, group). The file has either:
  !
  !    one copy   a single-population phase, or a file written before populations
  !               existed. The record is fanned out over every group -- one set of
  !               injection data for all families.
  !
  !    one copy   per (material, population) pair, in ATLAS order, block by block:
  !    per group  ATLAS_BCB::write_dp_bc loops materials and populations inside its
  !               block loop. Copy c then belongs to family c and nothing is shared.
  !
  !  With multiple copies, Setup_BC therefore walks block by block and expects
  !  the same boundary-cell order in every copy. Any other count is a mismatch
  !  between the phase file used to write the BC file and this run's groups; it
  !  stops setup rather than binding a family's data to the wrong copy.

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
  subroutine Setup_BC(grid, level)
    use ICE_Global_m, only: ngroups
    use ICE_Mod_MPI, only: mpi_is_root
    implicit none
    type(ICE_domain_type), intent(inout), target :: grid
    integer(I4), intent(in)                      :: level  !< Multigrid level, 1 = fine
    type(ICE_bc_type), dimension(:), pointer     :: bc
    integer(I4), dimension(:,:), pointer         :: nbc
    integer(I4) :: unitfile, ios
    integer(I4) :: n, ncell, p, idx, code
    integer(I4) :: b, c, r, nbb, nbmax, ncopy, nrec, cellbase, pfirst, plast
    integer(I4) :: hb, hi, hj, hk, hf
    integer(I4), allocatable :: key(:,:)   !< (i,j,k,f) of copy 1, to check the others against
    integer(I4) :: cs(9)
    real(R8)    :: massflux, velocity, temperature, radius, alpha, beta
    character(len=32) :: alpha_tok, beta_tok
    integer(I4) :: nchi(2), s, nzero_chim
    integer(I4), allocatable :: donorID(:,:)
    real(R8), allocatable :: weight(:)

    bc  => grid%bc
    nbc => grid%n_bf

    nconn = 0; nsym = 0; nio = 0; next = 0; nchim = 0; naxi = 0

    if (ngroups <= 0) error stop 'Setup_BC: ngroups not set'
    if (mod(size(bc), ngroups) /= 0) error stop 'Setup_BC: bc array is not a multiple of ngroups'

    ncell = size(bc) / ngroups
    nzero_chim = 0
    nrec = Count_BC_Records(level)

    if (mod(nrec, ncell) /= 0) then
      write(*,'(A,I0,A,I0,A)') '  [ICE::Setup_BC] the BC file holds ', nrec, &
        ' records, which is not a whole number of copies of the ', ncell, &
        ' boundary cells of this mesh'
      error stop
    endif

    ncopy = nrec / ncell

    if (ncopy /= 1 .and. ncopy /= ngroups) then
      write(*,'(A,I0,A,I0,A)') '  [ICE::Setup_BC] the BC file carries ', ncopy, &
        ' copies of the boundary table but this run has ', ngroups, ' particle families'
      write(*,'(A)') '  ATLAS writes one copy per (material, population) pair, counted from'
      write(*,'(A)') '  the <material> <npCP> lines of the phase file: declare one'
      write(*,'(A)') '  [ICE-Family*] section per pair, or use a single-population phase.'
      error stop
    endif

    open(newunit=unitfile, file=trim(bc_filename(level)), status='old', iostat=ios)

    if (ios /= 0) then
      write(*,'(A)') '  [ICE::Setup_BC] boundary-condition file not found: '//trim(bc_filename(level))
      error stop 'BC file not found'
    endif

    nbmax = 0

    do b = 1, grid%nb
      nbmax = max(nbmax, 2*(grid%blk(b)%dim(2)*grid%blk(b)%dim(3)    &
                          + grid%blk(b)%dim(1)*grid%blk(b)%dim(3)    &
                          + grid%blk(b)%dim(1)*grid%blk(b)%dim(2)))
    enddo

    allocate(key(4, nbmax))

    nbc = 0
    nrec     = 0
    cellbase = 0

    !> The file is walked the way ATLAS wrote it: block by block, and within a block
    !  one run of its boundary cells per copy. Reading it in that shape is what lets a
    !  record be attributed to a family without trusting its position in the file.
    do b = 1, grid%nb
      nbb = 2*(grid%blk(b)%dim(2)*grid%blk(b)%dim(3)    &
             + grid%blk(b)%dim(1)*grid%blk(b)%dim(3)    &
             + grid%blk(b)%dim(1)*grid%blk(b)%dim(2))
      do c = 1, ncopy
        do r = 1, nbb
          n    = cellbase + r
          nrec = nrec + 1
          read(unitfile,*,iostat=ios) hb, hi, hj, hk, hf, code
          if (ios /= 0) then
            write(*,'(A,I0)') '  [ICE::Setup_BC] error reading BC header record ', nrec
            error stop
          endif

          !> With one copy a record stands on its own -- every entry carries its own
          !  (b,i,j,k,f) -- so the file may be in any order, as it always could. With
          !  several, the order is the only thing that says which copy a record is in,
          !  so it has to be the order ATLAS writes: block by block, copies together,
          !  each copy over the same cells.
          if (ncopy > 1) then
            if (hb /= b) then
              write(*,'(A,I0,A,I0,A,I0)') '  [ICE::Setup_BC] record ', nrec, ' belongs to block ', &
                hb, ' where the mesh expects block ', b
              write(*,'(A)') '  a multi-population file must list the blocks in mesh order, with'
              write(*,'(A)') '  every copy of one block together, the way ATLAS writes them.'
              error stop
            endif
            if (c == 1) then
              key(1:4, r) = [hi, hj, hk, hf]
            elseif (any(key(1:4, r) /= [hi, hj, hk, hf])) then
              write(*,'(A,I0,A,I0,A)') '  [ICE::Setup_BC] copy ', c, ' of block ', b, &
                ' does not repeat the cells of copy 1 in the same order'
              error stop
            endif
          endif

          !> Same gate IGLOO uses: a pre-ATLAS 8-column record, or a property line swallowed as
          !  a header, lands outside {0, 100..999}. Lenient list-directed reads would otherwise
          !  mistype every face and corrupt the setup silently.
          if (code /= 0 .and. (code < 100 .or. code > 999)) then
            write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] bad BC code ', code, ' at record ', nrec
            write(*,'(A)') '  expected the ATLAS schema: 6 columns with a 3-digit code (or 0) in column 6.'
            error stop
          endif

          alpha = 0._R8; beta = 0._R8
          massflux = 0._R8; velocity = 0._R8; temperature = 0._R8; radius = 0._R8
          cs = 0

          !> Payload dispatch. Count_BC_Records walks the same rules to work out how
          !  many copies the file holds: a code added here needs adding there too.
          select case (code)
            case (101, 201)          ! connection / periodic
              nconn = nconn + 1
              read(unitfile,*,iostat=ios) cs(1:9)

            case (102)               ! chimera
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

            case (200)               ! axisymmetry, no payload
              naxi = naxi + 1

            case (300, 301)          ! symmetry / dispersed-phase wall, no payload
              nsym = nsym + 1

            case (400)               ! extrapolation, no payload
              next = next + 1

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
              nio = nio + 1
              read(unitfile,*,iostat=ios) massflux, velocity, alpha_tok, beta_tok, temperature, radius
              alpha = Parse_Dir_Token(alpha_tok)
              beta  = Parse_Dir_Token(beta_tok)

            case (0)
              continue

            case default
              write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] unsupported ATLAS BC code ', code, &
                                     ' at record ', nrec
              error stop
          end select

          if (ios /= 0) then
            write(*,'(A,I0)') '  [ICE::Setup_BC] error reading payload for record ', nrec
            error stop
          endif

          !> One boundary FACE per cell -- counted once, not once per group and not
          !  once per copy.
          if (c == 1) nbc(hb, hf) = nbc(hb, hf) + 1

          !> One copy: the record speaks for every family. One copy per family: it
          !  speaks for its own.
          if (ncopy == 1) then
            pfirst = 1 ; plast = ngroups
          else
            pfirst = c ; plast = c
          endif
          do p = pfirst, plast
            idx = (n-1)*ngroups + p
            bc(idx)%b = hb; bc(idx)%i = hi; bc(idx)%j = hj; bc(idx)%k = hk; bc(idx)%f = hf
            bc(idx)%p    = p
            bc(idx)%type = code
            bc(idx)%massflux    = massflux
            bc(idx)%velocity    = velocity
            bc(idx)%alpha       = alpha
            bc(idx)%beta        = beta
            bc(idx)%temperature = temperature
            bc(idx)%radius      = radius
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
        enddo
      enddo
      cellbase = cellbase + nbb
    enddo

    close(unitfile)

    if (nzero_chim > 0) then
      write(*,'(A,I0,A)') '  [ICE::Setup_BC] ', nzero_chim, &
        ' chimera records have a ghost layer whose donor weights sum to ~0 (no donor found by ATLAS BCB). '// &
        'Fix the BC file before running.'
      error stop
    endif

    if (mpi_is_root) then
      call Print_BC_Summary()
      call Check_Phase_Groups(ngroups)
    endif

  end subroutine Setup_BC


  !> How many records does the BC file hold? Walk the payload with the same rules as
  !  Setup_BC, storing nothing. The count divided by the mesh's boundary cells is the
  !  number of copies of the table, which has to be known before records can be
  !  attributed to families.
  integer(I4) function Count_BC_Records(level) result(nrec)
    implicit none
    integer(I4), intent(in) :: level
    integer(I4) :: unitfile, ios, s, nskip
    integer(I4) :: hb, hi, hj, hk, hf, code, nchi(2)
    real(R8)    :: first

    nrec = 0

    open(newunit=unitfile, file=trim(bc_filename(level)), status='old', iostat=ios)
    if (ios /= 0) error stop 'BC file not found'

    do
      read(unitfile,*,iostat=ios) hb, hi, hj, hk, hf, code
      if (ios /= 0) exit

      nrec = nrec + 1
      nskip = 0

      select case (code)
      case (101, 201, 401:403)
        nskip = 1
      case (102)
        read(unitfile,*,iostat=ios) nchi
        if (ios /= 0) exit
        nskip = sum(nchi)
      end select

      do s = 1, nskip
        read(unitfile,*,iostat=ios) first
        if (ios /= 0) exit
      enddo

      if (ios /= 0) exit
    enddo

    close(unitfile)

  end function Count_BC_Records

  !> ATLAS writes INPUT/<prefix>phase.txt for the condensed phase ICE solves: one line per material,
  !  "<name> <groups> [key=value ...]". ICE takes its group count from [ICE-Family*] and never
  !  opened this file, so a mismatch can go unnoticed. Warn when the two disagree.
  !  The file is optional (standalone ICE cases may not carry it): absent = silent.
  subroutine Check_Phase_Groups(ngroups)

    implicit none

    integer(I4), intent(in) :: ngroups
    integer(I4) :: u, ios, ng, nsum, nmat
    character(len=512) :: line
    character(len=64)  :: name

    open(newunit=u, file='INPUT/'//trim(ICE_phase_prefix)//'phase.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) return

    read(u, '(A)', iostat=ios) line
    nsum = 0
    nmat = 0

    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios) name, ng
      if (ios /= 0) cycle
      nmat = nmat + 1
      nsum = nsum + ng
    enddo

    close(u)

    if (nmat == 0) return

    if (nsum /= ngroups) then
      write(*,'(A,I0,A,I0,A)') '  [WARNING] INPUT/'//trim(ICE_phase_prefix)//'phase.txt declares ', nsum, &
        ' population(s) over ', nmat, ' material(s), but [ICE-Family*] defines '
      write(*,'(A,I0,A)') '            ', ngroups, ' group(s). The bc.txt carries one record block per (material,'// &
        ' population); the counts should agree.'
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

  !> Which file a level reads. ATLAS BCB writes <prefix>bc.txt for the fine level and
  !  <prefix>bc<level>.txt for each coarse one -- write_ig_bc, write_sp_bc and
  !  write_dp_bc all take the level and branch the same way -- and MOSE reads the same
  !  two names. A coarse level is a different mesh, so it needs its own boundary table:
  !  the records are per cell, and a coarse cell is not a fine one.
  function bc_filename(level) result(fname)
    implicit none
    integer(I4), intent(in) :: level
    character(len=256)      :: fname
    character(len=16)       :: lvl

      if (level <= 1) then
        fname = 'INPUT/'//trim(ICE_phase_prefix)//'bc.txt'
      else
        write(lvl,'(I0)') level
        fname = 'INPUT/'//trim(ICE_phase_prefix)//'bc'//trim(lvl)//'.txt'
      endif
  end function bc_filename

end module ICE_IO_BC
