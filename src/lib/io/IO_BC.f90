module ICE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ICE_Global_m,         only: ICE_phase_prefix
  use ICE_Advanced_Types_m, only: ICE_domain_type, ICE_bc_type
  implicit none
  private

  public :: Setup_BC, Print_BC_Summary

  integer :: nconn = 0, nsym = 0, nio = 0, next = 0, nchim = 0, naxi = 0

  !> A payload field still holding its sentinel after the read of its line was left empty:
  !  a null value (`,,`) leaves a list-directed item unchanged.
  real(R8),    parameter :: EMPTY  = -huge(1._R8)
  integer(I4), parameter :: IEMPTY = -huge(1_I4)

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
  !  Every header and payload line is read whole and on its own (next_line), so a
  !  short line cannot take its missing values from the next one. A line that is
  !  short or holds a value that is not a number, an empty value (`,,`), a value
  !  that is not finite, a direction that is neither a number nor `normal`, a cell
  !  outside the mesh or off its face stops the setup, naming the file, the line
  !  and the record.

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
    integer(I4) :: nline                     !< lines of the file read so far (next_line)
    character(len=:), allocatable :: line

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
    nline    = 0

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
          call next_line(unitfile, line, ios, nline)
          if (ios == 0) read(line,*,iostat=ios) hb, hi, hj, hk, hf, code
          if (ios /= 0) call line_error(level, nline, nrec, -1, line)

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
            write(*,'(A,I0,A,I0,A,I0)') '  [ICE::Setup_BC] bad BC code ', code, ' at record ', nrec, &
              ', line ', nline
            write(*,'(A)') '  expected the ATLAS schema: 6 columns with a 3-digit code (or 0) in column 6.'
            error stop
          endif

          !> The record's cell is a cell of the mesh, against the face it names: every
          !  later use indexes the block, its cells and their ghosts with it.
          call check_cell(grid, level, nline, nrec, code, 'its cell', hb, hi, hj, hk, hf)

          alpha = 0._R8; beta = 0._R8
          massflux = 0._R8; velocity = 0._R8; temperature = 0._R8; radius = 0._R8
          cs = 0

          !> Payload dispatch. Count_BC_Records walks the same rules to work out how
          !  many copies the file holds: a code added here needs adding there too.
          !  Every field is first set to a sentinel (EMPTY, IEMPTY, blank): one still
          !  holding it after the read was left empty.
          select case (code)
            case (101, 201)          ! connection / periodic
              nconn = nconn + 1
              cs = IEMPTY
              call next_line(unitfile, line, ios, nline)
              if (ios == 0) read(line,*,iostat=ios) cs(1:9)
              if (ios == 0) call need_integers(level, nline, nrec, code, cs, &
                [character(len=3) :: 'bs', 'is', 'js', 'ks', 'fs', 'd11', 'd12', 'd21', 'd22'])

            case (102)               ! chimera
              !> Payload: `nchi_g1 nchi_g2`, then nchi_g1 + nchi_g2 donor lines `b i j k weight`,
              !  ghost layer 1 first. Weights are normalised per layer by ATLAS.
              nchim = nchim + 1
              nchi = IEMPTY
              call next_line(unitfile, line, ios, nline)
              if (ios == 0) read(line,*,iostat=ios) nchi
              if (ios == 0) call need_integers(level, nline, nrec, code, nchi, [character(len=7) :: 'nchi_g1', 'nchi_g2'])
              if (ios == 0) then
                if (allocated(donorID)) deallocate(donorID, weight)
                allocate(donorID(sum(nchi),4), weight(sum(nchi)))
                do s = 1, sum(nchi)
                  donorID(s,1:4) = IEMPTY ; weight(s) = EMPTY
                  call next_line(unitfile, line, ios, nline)
                  if (ios == 0) read(line,*,iostat=ios) donorID(s,1:4), weight(s)
                  if (ios /= 0) exit
                  call need_integers(level, nline, nrec, code, donorID(s,1:4), [character(len=7) :: 'donor b', &
                                     'donor i', 'donor j', 'donor k'])
                  call need_reals(level, nline, nrec, code, weight(s:s), ['weight'])
                  call check_cell(grid, level, nline, nrec, code, 'its donor', donorID(s,1), donorID(s,2), &
                                  donorID(s,3), donorID(s,4), 0)
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
              massflux = EMPTY ; velocity = EMPTY ; temperature = EMPTY ; radius = EMPTY
              alpha_tok = ' ' ; beta_tok = ' '
              call next_line(unitfile, line, ios, nline)
              if (ios == 0) read(line,*,iostat=ios) massflux, velocity, alpha_tok, beta_tok, temperature, radius
              if (ios == 0) call check_inlet(level, nline, nrec, code, massflux, velocity, alpha_tok, beta_tok, &
                                             temperature, radius, alpha, beta)

            case (0)
              continue

            case default
              write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] unsupported ATLAS BC code ', code, &
                                     ' at record ', nrec
              error stop
          end select

          if (ios /= 0) call line_error(level, nline, nrec, code, line)

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
    integer(I4) :: unitfile, ios, s, nskip, nline
    integer(I4) :: hb, hi, hj, hk, hf, code, nchi(2)
    character(len=:), allocatable :: line

    nrec  = 0
    nline = 0

    open(newunit=unitfile, file=trim(bc_filename(level)), status='old', iostat=ios)
    if (ios /= 0) error stop 'BC file not found'

    !> Line by line, as Setup_BC reads it: the end of the file ends the count, a line
    !  that is not a header where one is due, or a payload cut short, stops the setup.
    do
      call next_line(unitfile, line, ios, nline)
      if (ios /= 0) exit
      read(line,*,iostat=ios) hb, hi, hj, hk, hf, code
      if (ios /= 0) call line_error(level, nline, nrec + 1, -1, line)

      nrec = nrec + 1
      nskip = 0

      select case (code)
      case (101, 201, 401:403)
        nskip = 1
      case (102)
        nchi = IEMPTY
        call next_line(unitfile, line, ios, nline)
        if (ios == 0) read(line,*,iostat=ios) nchi
        if (ios /= 0) call line_error(level, nline, nrec, code, line)
        call need_integers(level, nline, nrec, code, nchi, [character(len=7) :: 'nchi_g1', 'nchi_g2'])
        nskip = sum(nchi)
      end select

      do s = 1, nskip
        call next_line(unitfile, line, ios, nline)
        if (ios /= 0) call line_error(level, nline, nrec, code, line)
      enddo
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
  !  huge() sentinel. A token that is neither `normal` nor a number gives ios /= 0.
  pure subroutine Parse_Dir_Token(tok, val, ios)
    implicit none
    character(len=*), intent(in)  :: tok
    real(R8),         intent(out) :: val
    integer(I4),      intent(out) :: ios

    if (trim(adjustl(tok)) == 'normal') then
      val = 0._R8
      ios = 0
    else
      read(tok, *, iostat=ios) val
      if (ios /= 0) val = 0._R8
    endif

  end subroutine Parse_Dir_Token


  !> The fields of an inlet line (401-403), read: each number finite and given, each direction
  !  `normal` or a finite number. With a mass flux, the ghost divides by the speed and by the
  !  radius cubed (Lib_Ghost), so both must be non-zero there and the radius positive; with
  !  none the inlet is a symmetry and neither is used.
  subroutine check_inlet(level, nline, nrec, code, massflux, velocity, alpha_tok, beta_tok, temperature, radius, &
                         alpha, beta)
    implicit none
    integer(I4),      intent(in)  :: level, nline, nrec, code
    real(R8),         intent(in)  :: massflux, velocity, temperature, radius
    character(len=*), intent(in)  :: alpha_tok, beta_tok
    real(R8),         intent(out) :: alpha, beta
    integer(I4) :: ios

    call need_reals(level, nline, nrec, code, [massflux, velocity, temperature, radius], &
                    [character(len=11) :: 'massflux', 'velocity', 'temperature', 'radius'])
    if (len_trim(alpha_tok) == 0) call field_error(level, nline, nrec, code, 'the value of alpha is empty')
    if (len_trim(beta_tok)  == 0) call field_error(level, nline, nrec, code, 'the value of beta is empty')
    call Parse_Dir_Token(alpha_tok, alpha, ios)
    if (ios /= 0) call field_error(level, nline, nrec, code, &
      "alpha is neither a number nor normal: '"//trim(adjustl(alpha_tok))//"'")
    call Parse_Dir_Token(beta_tok, beta, ios)
    if (ios /= 0) call field_error(level, nline, nrec, code, &
      "beta is neither a number nor normal: '"//trim(adjustl(beta_tok))//"'")
    call need_reals(level, nline, nrec, code, [alpha, beta], [character(len=5) :: 'alpha', 'beta'])
    if (massflux /= 0._R8) then
      if (velocity == 0._R8) call field_error(level, nline, nrec, code, &
        'the velocity is 0 with a mass flux: the injected density would be infinite')
      if (.not. (radius > 0._R8)) call field_error(level, nline, nrec, code, &
        'the radius is not positive with a mass flux: the injected number density would be infinite or negative')
    endif

  end subroutine check_inlet


  !> The reals v of a payload line, in the order of names: each one finite, then given (not
  !  left at EMPTY by a null value). Finite first: a NaN fails every comparison, so the EMPTY
  !  test alone would let one through.
  subroutine need_reals(level, nline, nrec, code, v, names)
    implicit none
    integer(I4),      intent(in) :: level, nline, nrec, code
    real(R8),         intent(in) :: v(:)
    character(len=*), intent(in) :: names(:)
    integer(I4) :: k

    do k = 1, size(v)
      if (.not. ieee_is_finite(v(k))) call field_error(level, nline, nrec, code, &
        'the value of '//trim(names(k))//' is not finite (NaN or infinity)')
    enddo
    do k = 1, size(v)
      if (v(k) == EMPTY) call field_error(level, nline, nrec, code, 'the value of '//trim(names(k))//' is empty')
    enddo

  end subroutine need_reals


  !> The integers v of a payload line, in the order of names: each one given (not left at
  !  IEMPTY by a null value).
  subroutine need_integers(level, nline, nrec, code, v, names)
    implicit none
    integer(I4),      intent(in) :: level, nline, nrec, code
    integer(I4),      intent(in) :: v(:)
    character(len=*), intent(in) :: names(:)
    integer(I4) :: k

    do k = 1, size(v)
      if (v(k) == IEMPTY) call field_error(level, nline, nrec, code, 'the value of '//trim(names(k))//' is empty')
    enddo

  end subroutine need_integers


  !> Stop unless (b,i,j,k) is a cell of the mesh and, for f /= 0, lies against its face f:
  !  i = 1 for face 1, i = dim(1) for face 2, j = 1 for face 3, and so on. what names the
  !  cell in the message.
  subroutine check_cell(grid, level, nline, nrec, code, what, b, i, j, k, f)
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer(I4),           intent(in) :: level, nline, nrec, code, b, i, j, k, f
    character(len=*),      intent(in) :: what
    character(len=160) :: msg
    logical :: inside

    inside = b >= 1 .and. b <= grid%nb
    if (inside) inside = i >= 1 .and. i <= grid%blk(b)%dim(1) .and. j >= 1 .and. j <= grid%blk(b)%dim(2) &
                   .and. k >= 1 .and. k <= grid%blk(b)%dim(3)
    if (.not. inside) then
      write(msg,'(A,4(1X,I0),A)') what//' (block i j k)', b, i, j, k, ' is not a cell of the mesh'
      call field_error(level, nline, nrec, code, trim(msg))
    endif
    if (f == 0) return

    select case (f)
      case (1) ; inside = i == 1
      case (2) ; inside = i == grid%blk(b)%dim(1)
      case (3) ; inside = j == 1
      case (4) ; inside = j == grid%blk(b)%dim(2)
      case (5) ; inside = k == 1
      case (6) ; inside = k == grid%blk(b)%dim(3)
      case default
        write(msg,'(A,I0,A)') 'face ', f, ' is not a face (1 to 6)'
        call field_error(level, nline, nrec, code, trim(msg))
    end select
    if (.not. inside) then
      write(msg,'(A,4(1X,I0),A,I0)') what//' (block i j k)', b, i, j, k, ' does not lie against face ', f
      call field_error(level, nline, nrec, code, trim(msg))
    endif

  end subroutine check_cell


  !> The next line of the BC file that is not blank, whole, whatever its length (read in
  !  pieces without advancing), as a list-directed read skips blank lines. nline counts
  !  every line read, blank ones included, so it is the line number in the file. ios < 0
  !  at the end of the file; a last line without its end-of-line mark is still a line.
  subroutine next_line(unit, line, ios, nline)
    implicit none
    integer(I4),                   intent(in)    :: unit
    character(len=:), allocatable, intent(out)   :: line
    integer(I4),                   intent(out)   :: ios
    integer(I4),                   intent(inout) :: nline
    character(len=256) :: piece
    integer(I4)        :: nread

    do
      line = ''
      do
        read(unit, '(A)', advance='no', iostat=ios, size=nread) piece
        line = line//piece(1:nread)
        if (ios /= 0) exit
      enddo
      if (is_iostat_eor(ios)) ios = 0
      if (is_iostat_end(ios) .and. len(line) > 0) ios = 0
      if (ios /= 0) exit
      nline = nline + 1
      if (len_trim(line) > 0) exit
    enddo

  end subroutine next_line


  !> Stop on line nline of the BC file of this level, in record nrec, which cannot be read as
  !  the values it must hold: its header (code < 0) or a payload line of that code. A blank
  !  line is the end of the file.
  subroutine line_error(level, nline, nrec, code, line)
    implicit none
    integer(I4),      intent(in) :: level, nline, nrec, code
    character(len=*), intent(in) :: line

    if (len_trim(line) == 0) then
      write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] '//trim(bc_filename(level))//' ends at line ', nline, &
        ', before the end of record ', nrec
    elseif (code < 0) then
      write(*,'(A,I0,A,I0,A)') '  [ICE::Setup_BC] '//trim(bc_filename(level))//' line ', nline, ' (record ', nrec, &
        '): the line is not block, i, j, k, face and code'
    else
      write(*,'(A,I0,A,I0,A,I0,A)') '  [ICE::Setup_BC] '//trim(bc_filename(level))//' line ', nline, &
        ' (record ', nrec, ', code ', code, '): the line is short or holds a value that is not a number'
    endif
    if (len_trim(line) > 0) write(*,'(A)') '  "'//trim(line(1:min(len(line), 200)))//'"'
    error stop 1

  end subroutine line_error


  !> Stop on line nline of the BC file of this level, in record nrec of that code, whose
  !  values were read but one of them cannot be used: why says which and what is wrong.
  subroutine field_error(level, nline, nrec, code, why)
    implicit none
    integer(I4),      intent(in) :: level, nline, nrec, code
    character(len=*), intent(in) :: why

    write(*,'(A,I0,A,I0,A,I0,A)') '  [ICE::Setup_BC] '//trim(bc_filename(level))//' line ', nline, &
      ' (record ', nrec, ', code ', code, '): '//why
    error stop 1

  end subroutine field_error

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
