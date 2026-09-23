module ICE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Global_m,         only: ICE_phase_prefix
  use ICE_Advanced_Types_m, only: ICE_domain_type, ICE_bc_type

  implicit none
  private
  public :: Setup_BC, Print_BC_Summary

  integer :: nconn = 0, nsym = 0, nio = 0, next = 0, nchim = 0

contains


  !> Read INPUT/<prefix>bc.txt, written by ATLAS BCB. Header per boundary cell:
  !
  !      b  i  j  k  f  code
  !
  !  followed by a code-dependent payload. ONE record per boundary cell; ICE's bc array
  !  holds one entry per (cell, group), so each record is fanned out over the groups with
  !  bc%p set accordingly. bc%type stores the ATLAS code as is:
  !
  !      0        null                   no ghost, no boundary flux
  !      101/201  connection / periodic
  !      102      chimera
  !      200      axisymmetry            no ghost, no boundary flux
  !      300      symmetry
  !      400      extrapolation
  !      401-403  inlet (see Read payload below)
  subroutine Setup_BC(grid)
    use ICE_Global_m, only: ngroups
    use ICE_Mod_MPI, only: mpi_is_root
    implicit none
    type(ICE_domain_type), intent(inout), target :: grid
    type(ICE_bc_type), dimension(:), pointer     :: bc
    integer(I4), dimension(:,:), pointer         :: nbc
    integer(I4) :: unitfile, ios
    integer(I4) :: n, ncell, p, idx, code
    integer(I4) :: hb, hi, hj, hk, hf
    integer(I4) :: cs(9)
    real(R8)    :: massflux, velocity, temperature, radius, alpha, beta
    character(len=32) :: alpha_tok, beta_tok
    integer(I4) :: nchi(2), s, nzero_chim
    integer(I4),  allocatable :: donorID(:,:)
    real(R8),     allocatable :: weight(:)

    bc  => grid%bc
    nbc => grid%n_bf
    nconn = 0; nsym = 0; nio = 0; next = 0; nchim = 0

    if (ngroups <= 0) error stop 'Setup_BC: ngroups not set'
    if (mod(size(bc), ngroups) /= 0) error stop 'Setup_BC: bc array is not a multiple of ngroups'
    ncell = size(bc) / ngroups
    nzero_chim = 0

    open(newunit=unitfile, file='INPUT/'//trim(ICE_phase_prefix)//'bc.txt', &
         status='old', iostat=ios)
    if (ios /= 0) error stop 'BC file not found'

    nbc = 0
    do n = 1, ncell

      read(unitfile,*,iostat=ios) hb, hi, hj, hk, hf, code
      if (ios /= 0) then
        write(*,'(A,I0)') '  [ICE::Setup_BC] error reading BC header record ', n
        error stop
      endif

      !> Same gate IGLOO uses: a pre-ATLAS 8-column record, or a property line swallowed as
      !  a header, lands outside {0, 100..999}. Lenient list-directed reads would otherwise
      !  mistype every face and corrupt the setup silently.
      if (code /= 0 .and. (code < 100 .or. code > 999)) then
        write(*,'(A,I0,A,I0)') '  [ICE::Setup_BC] bad BC code ', code, ' at record ', n
        write(*,'(A)') '  expected the ATLAS schema: 6 columns with a 3-digit code (or 0) in column 6.'
        error stop
      endif

      alpha = 0._R8; beta = 0._R8
      massflux = 0._R8; velocity = 0._R8; temperature = 0._R8; radius = 0._R8
      cs = 0

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
          continue

        case (300)               ! symmetry, no payload
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

    end do

    close(unitfile)

    if (nzero_chim > 0) then
      write(*,'(A,I0,A)') '  [ICE::Setup_BC] ', nzero_chim, ' chimera records have a ghost layer '// &
        'whose donor weights sum to ~0 (no donor found by ATLAS BCB). Fix the BC file before running.'
      error stop
    endif

    if (mpi_is_root) call Print_BC_Summary()

  end subroutine Setup_BC


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
    if (nio   > 0) write(*,'(A,T35,I0)') '   Inflow',        nio
    if (next  > 0) write(*,'(A,T35,I0)') '   Extrapolation', next
    if (nchim > 0) write(*,'(A,T35,I0)') '   Chimera',       nchim
    write(*,*)
  end subroutine Print_BC_Summary


end module ICE_IO_BC
