module ICE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ICE_Global_m,         only: ICE_phase_prefix
  use ICE_Advanced_Types_m, only: ICE_domain_type, ICE_bc_type
  use ICE_Config_Types_m,   only: obj_sim_param
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
    logical     :: normal_dir
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
              alpha_tok = ' ' ; beta_tok = ' ' ; normal_dir = .false.
              call next_line(unitfile, line, ios, nline)
              if (ios == 0) read(line,*,iostat=ios) massflux, velocity, alpha_tok, beta_tok, temperature, radius
              if (ios == 0) call check_inlet(level, nline, nrec, code, massflux, velocity, alpha_tok, beta_tok, &
                                             temperature, radius, alpha, beta, normal_dir)
              !> 401 and 403 take the speed of what they inject from the gas (401 also its density and
              !  temperature): a run without a gas phase has none to read
              if (ios == 0 .and. code /= 402 .and. massflux /= 0._R8 .and. &
                  .not. (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled)) &
                call field_error(level, nline, nrec, code, 'the inlet takes its speed from the gas, and this run '// &
                  'has no gas phase: use 402, which takes the speed in m/s')

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
            bc(idx)%normal_dir  = normal_dir
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

    call Pair_Connections(grid, level)

    if (mpi_is_root) then
      call Print_BC_Summary()
      call Check_Phase_Groups(ngroups)
    endif

  end subroutine Setup_BC


  !> The two records of each connection face (101, 201) of a family, one per block, must
  !  point at each other: each is the other's partner. One of them owns the face: the record
  !  whose face is the high-index one of its cell (2, 4, 6) when the other's is a low-index
  !  one, as on an interior face, whose flux takes the cell below; between faces of the same
  !  parity (blocks that meet with their index directions permuted) the record with the
  !  smaller (block, face, k, j, i), a key each record has of its own, so that the owner does
  !  not depend on the order of the file. A record that points at a cell and face where no
  !  connection record of its family stands, at its own cell and face, or at a record that
  !  does not point back at it stops the setup. So does one whose block across has fewer than
  !  gc cells across its face at this level: the ghost cells of the face are the first gc
  !  cells of that block, and the second would be one of its ghost cells. The search goes
  !  through a map, per block face, from the two indices of a cell along the face to its
  !  record: one pass over the records, in any order.
  subroutine Pair_Connections(grid, level)
    use ICE_Global_m, only: ngroups, gc
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(I4),           intent(in)    :: level  !< Multigrid level, 1 = fine
    type :: face_records_type
      integer(I4), allocatable :: r(:,:)
    end type face_records_type
    type(face_records_type), allocatable :: rec(:,:)
    integer(I4) :: p, n, q, b, f, t(2), nt(2), across, nf
    character(len=40)  :: at_level
    character(len=200) :: msg

    at_level = ''
    if (level > 1) write(at_level,'(A,I0)') ' at multigrid level ', level

    allocate(rec(grid%nb, 6))
    if (allocated(grid%conn)) deallocate(grid%conn)
    allocate(grid%conn(ngroups))

    do p = 1, ngroups

      !> Map of the connection records of family p, per block face
      do n = p, size(grid%bc), ngroups
        if (grid%bc(n)%type /= 101 .and. grid%bc(n)%type /= 201) cycle
        b = grid%bc(n)%b ; f = grid%bc(n)%f
        if (.not. allocated(rec(b,f)%r)) then
          nt = tangential(f, grid%blk(b)%dim(1), grid%blk(b)%dim(2), grid%blk(b)%dim(3))
          allocate(rec(b,f)%r(nt(1), nt(2)))
          rec(b,f)%r = 0
        endif
        t = tangential(f, grid%bc(n)%i, grid%bc(n)%j, grid%bc(n)%k)
        rec(b,f)%r(t(1), t(2)) = n
      enddo

      !> Partner and owner of each
      do n = p, size(grid%bc), ngroups
        if (grid%bc(n)%type /= 101 .and. grid%bc(n)%type /= 201) cycle
        b = grid%bc(n)%bs ; f = grid%bc(n)%fs
        q = 0
        if (b >= 1 .and. b <= grid%nb .and. f >= 1 .and. f <= 6) then
          if (allocated(rec(b,f)%r)) then
            t = tangential(f, grid%bc(n)%is, grid%bc(n)%js, grid%bc(n)%ks)
            if (all(t >= 1) .and. all(t <= shape(rec(b,f)%r))) q = rec(b,f)%r(t(1), t(2))
          endif
        endif
        if (q /= 0) then
          if (grid%bc(q)%i /= grid%bc(n)%is .or. grid%bc(q)%j /= grid%bc(n)%js .or. &
              grid%bc(q)%k /= grid%bc(n)%ks) q = 0
        endif
        if (q == 0) then
          call pair_error(n, 'points at '//trim(cell_face(grid%bc(n)%bs, grid%bc(n)%is, grid%bc(n)%js, &
            grid%bc(n)%ks, grid%bc(n)%fs))//', where no connection record of its family stands')
        endif
        if (q == n) call pair_error(n, 'points at its own cell and face')
        if (grid%bc(q)%type /= grid%bc(n)%type .or. &
            grid%bc(q)%bs /= grid%bc(n)%b .or. grid%bc(q)%fs /= grid%bc(n)%f .or. grid%bc(q)%is /= grid%bc(n)%i .or. &
            grid%bc(q)%js /= grid%bc(n)%j .or. grid%bc(q)%ks /= grid%bc(n)%k) then
          write(msg,'(A,I0,A)') ', whose record (code ', grid%bc(q)%type, ') points at '
          call pair_error(n, 'points at '//trim(cell_face(b, grid%bc(q)%i, grid%bc(q)%j, grid%bc(q)%k, f))// &
            trim(msg)//' '//trim(cell_face(grid%bc(q)%bs, grid%bc(q)%is, grid%bc(q)%js, grid%bc(q)%ks, &
            grid%bc(q)%fs))//': the two records of a face must have one code and point at each other')
        endif
        across = grid%blk(b)%dim((f+1)/2)
        if (across < gc) then
          write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') 'has across it block ', b, ', with ', across, ' cell across its face ', f, &
            trim(at_level)//', fewer than the ', gc, ' a connection needs'
          call pair_error(n, trim(msg))
        endif

        grid%bc(n)%partner = q
        if (mod(grid%bc(n)%f, 2) == 0 .and. mod(grid%bc(q)%f, 2) == 1) then
          grid%bc(n)%owner = .true.
        elseif (mod(grid%bc(n)%f, 2) == 1 .and. mod(grid%bc(q)%f, 2) == 0) then
          grid%bc(n)%owner = .false.
        else
          grid%bc(n)%owner = key_less([grid%bc(n)%b, grid%bc(n)%f, grid%bc(n)%k, grid%bc(n)%j, grid%bc(n)%i], &
                                      [grid%bc(q)%b, grid%bc(q)%f, grid%bc(q)%k, grid%bc(q)%j, grid%bc(q)%i])
        endif
      enddo

      !> The faces of family p, numbered by their owner records in table order: every rank reads the whole
      !  table, so every rank numbers them alike
      nf = count([( (grid%bc(n)%type == 101 .or. grid%bc(n)%type == 201) .and. grid%bc(n)%owner, &
                    n = p, size(grid%bc), ngroups )])
      grid%conn(p)%n = nf
      allocate(grid%conn(p)%prt(max(nf, 1)))
      nf = 0
      do n = p, size(grid%bc), ngroups
        if (grid%bc(n)%type /= 101 .and. grid%bc(n)%type /= 201) cycle
        if (.not. grid%bc(n)%owner) cycle
        nf = nf + 1
        grid%bc(n)%cf = nf
        grid%bc(grid%bc(n)%partner)%cf = nf
        grid%conn(p)%prt(nf) = grid%bc(n)%partner
      enddo

      do b = 1, grid%nb
        do f = 1, 6
          if (allocated(rec(b,f)%r)) deallocate(rec(b,f)%r)
        enddo
      enddo
    enddo

  contains

    !> The two indices of a cell (or the two dimensions of a block) along face f
    pure function tangential(f, i, j, k) result(t)
      integer(I4), intent(in) :: f, i, j, k
      integer(I4) :: t(2)
      select case (f)
        case (1:2)    ; t = [j, k]
        case (3:4)    ; t = [i, k]
        case default  ; t = [i, j]
      end select
    end function tangential

    !> Whether key1 comes before key2, comparing (block, face, k, j, i) in that order
    pure logical function key_less(key1, key2)
      integer(I4), intent(in) :: key1(5), key2(5)
      integer(I4) :: m
      key_less = .false.
      do m = 1, 5
        if (key1(m) /= key2(m)) then
          key_less = key1(m) < key2(m)
          return
        endif
      enddo
    end function key_less

    !> "block b, cell i j k, face f"
    function cell_face(b, i, j, k, f) result(s)
      integer(I4), intent(in) :: b, i, j, k, f
      character(len=80) :: s
      write(s,'(A,I0,A,3(1X,I0),A,I0)') 'block ', b, ', cell', i, j, k, ', face ', f
    end function cell_face

    !> Stop on the connection record n of this level: why says what is wrong with it
    subroutine pair_error(n, why)
      integer(I4),      intent(in) :: n
      character(len=*), intent(in) :: why
      character(len=24) :: family
      family = ''
      if (ngroups > 1) write(family,'(A,I0)') ', family ', grid%bc(n)%p
      write(*,'(A,I0,A)') '  [ICE::Setup_BC] '//trim(bc_filename(level))//': the connection record of '// &
        trim(cell_face(grid%bc(n)%b, grid%bc(n)%i, grid%bc(n)%j, grid%bc(n)%k, grid%bc(n)%f))//trim(family)// &
        ' (code ', grid%bc(n)%type, ') '//why
      error stop 1
    end subroutine pair_error

  end subroutine Pair_Connections


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
  !  face normal: it maps to 0 here, and check_inlet sets the record's normal_dir flag, which
  !  Lib_Ghost reads (numbers are angles, 0 and 0 included). A token that is neither `normal`
  !  nor a number gives ios /= 0.
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
                         alpha, beta, normal_dir)
    implicit none
    integer(I4),      intent(in)  :: level, nline, nrec, code
    real(R8),         intent(in)  :: massflux, velocity, temperature, radius
    character(len=*), intent(in)  :: alpha_tok, beta_tok
    real(R8),         intent(out) :: alpha, beta
    logical,          intent(out) :: normal_dir
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
    !> The face normal fixes both angles: `normal` for one and a number for the other says nothing
    normal_dir = (trim(adjustl(alpha_tok)) == 'normal')
    if (normal_dir .neqv. (trim(adjustl(beta_tok)) == 'normal')) call field_error(level, nline, nrec, code, &
      "alpha and beta are both normal or both numbers: '"//trim(adjustl(alpha_tok))//"', '"// &
      trim(adjustl(beta_tok))//"'")
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
