!> The components of the condensed state that are odd across mirror faces one cell apart, set to zero
!> at set-up and not integrated.
module ICE_Mod_Freeze
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: setup_freeze, nfrozen, frozen

  integer(kind=I4), allocatable :: nfrozen(:)   ! slots of family p not integrated: frozen(1:nfrozen(p),p)
  integer(kind=I4), allocatable :: frozen(:,:)  ! the same slot in the primitive and the conservative state

contains

  !> Between two faces that mirror the flow one cell apart (the k faces of a 2-D mesh, the j and k faces
  !> of a 1-D one) the components of the state odd across them, w across the k faces (with AG also P13
  !> and P23, whose conservative slots are rho u w + P13 and rho v w + P23), have the mode
  !> dw/dt = -2 a w / dz, a the signal speed and dz the width of the cell across the faces: the flux
  !> between the cell and its mirror image damps them. The time step follows the directions the mesh
  !> resolves (compute_dt), and on a layer thin against its other cells that mode leaves the real
  !> stability interval of the scheme, [-2, 0] for Heun and [-2.51, 0] for RK3: round-off in w grows
  !> until the state is unphysical. The flow between two mirror planes has those components zero, so
  !> they are set to zero here, compute_source zeroes their source (the gas receives none either) and
  !> compute_residual their residual:
  !>  - a 2-D plane mesh whose k faces are all 200, 300 or 301 and whose k node planes are planes normal
  !>    to z: w;
  !>  - a 2-D axisymmetric wedge whose k faces are all 200 and whose k node planes are the planes at
  !>    -/+ delthe/2 about z = 0: w, the swirl of the cell centres, which lie on z = 0. Any other wedge
  !>    stops the setup;
  !>  - a 1-D mesh: w as on a plane mesh, and v likewise between mirror j faces normal to y (with AG
  !>    also P12 and P23).
  !> A plane to within mesh_fit of the thickness of its block (mirror_planes). Mirror faces that are not
  !> planes normal to their axis, as on a layer between parallel planes inclined to z, ask the flow to
  !> follow them across the layer: the component is then integrated, with a warning, and a thin layer
  !> can still fail. The codes are those of every family on the fine level, and every family is frozen
  !> alike. Called once the initial field is in place, before the coarse levels take it.
  subroutine setup_freeze (grid)
    use ICE_Advanced_Types_m
    use ICE_Global_m,    only: ndir, ngroups, nbase, ncond_max
    use ICE_Mod_Metrics, only: delthe, mirror_planes
    use ICE_Mod_MPI,     only: mpi_is_root
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    logical          :: freeze_v, freeze_w, mirror, ag, odd(ncond_max)
    integer(kind=I4) :: n, b, p, s, nj, nk, nmj, nmk, other
    real(kind=R8)    :: fit
    character(len=64) :: what, faces

    !> The j and k faces, and how many of them mirror
    nj = 0 ; nk = 0 ; nmj = 0 ; nmk = 0 ; other = 0
    do n = 1, size(grid%bc)
      select case (grid%bc(n)%f)
      case (3, 4)
        nj = nj + 1
        if (any(grid%bc(n)%type == [200, 300, 301])) nmj = nmj + 1
      case (5, 6)
        nk = nk + 1
        if (any(grid%bc(n)%type == [200, 300, 301])) nmk = nmk + 1
        if (grid%bc(n)%type /= 200 .and. other == 0) other = n
      end select
    end do

    freeze_v = .false.
    freeze_w = .false.
    if (ndir == 2 .and. delthe /= 0._R8) then
      if (other > 0) then
        write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] [ICE::setup_freeze] face ', grid%bc(other)%f, &
          ' of the cell (block i j)', grid%bc(other)%b, grid%bc(other)%i, grid%bc(other)%j, ' carries ', &
          grid%bc(other)%type, ': the k faces of an axisymmetric wedge must carry 200, the mirror that keeps '// &
          'the flow axisymmetric'
        error stop 1
      end if
      call mirror_planes(grid, 3, mirror, fit)
      if (.not. mirror) then
        write(*,'(A,ES8.2,A)') ' [ERROR] [ICE::setup_freeze] the k node planes of the axisymmetric wedge are not '// &
          'the planes at -/+ delthe/2 about z = 0 (a node off its plane by ', fit, ' of the thickness of its '// &
          'block): the wedge must be symmetric about the plane z = 0'
        error stop 1
      end if
      freeze_w = .true.
    elseif (ndir <= 2) then
      if (nk > 0 .and. nmk == nk) then
        call mirror_planes(grid, 3, mirror, fit)
        freeze_w = mirror
        if (.not. mirror .and. mpi_is_root) write(*,'(A,ES8.2,A)') ' [WARNING] [ICE::setup_freeze] mirror k '// &
          'faces that are not planes normal to z (a node off by ', fit, ' of the thickness of its block): w is '// &
          'integrated, and its mode across the cells can make a thin layer fail'
      end if
      if (ndir == 1 .and. nj > 0 .and. nmj == nj) then
        call mirror_planes(grid, 2, mirror, fit)
        freeze_v = mirror
        if (.not. mirror .and. mpi_is_root) write(*,'(A,ES8.2,A)') ' [WARNING] [ICE::setup_freeze] mirror j '// &
          'faces that are not planes normal to y (a node off by ', fit, ' of the thickness of its block): v is '// &
          'integrated, and its mode across the cells can make thin cells fail'
      end if
    end if

    !> The slots of each family: v and w, and with AG the tensor components odd across the faces
    if (allocated(nfrozen)) deallocate(nfrozen, frozen)
    allocate(nfrozen(ngroups), frozen(ncond_max, ngroups))
    nfrozen = 0
    frozen  = 0
    do p = 1, ngroups
      odd    = .false.
      odd(3) = freeze_v
      odd(4) = freeze_w
      if (nbase(p) == 12) then
        odd(6) = freeze_v
        odd(7) = freeze_w
        odd(9) = freeze_v .or. freeze_w
      end if
      nfrozen(p) = count(odd)
      frozen(1:nfrozen(p),p) = pack([(s, s = 1, ncond_max)], odd)
      do b = 1, grid%nb
        do s = 1, nfrozen(p)
          grid%blk(b)%cond_phase(p)%prim(frozen(s,p),:,:,:) = 0._R8
        end do
      end do
    end do

    if (.not. (freeze_v .or. freeze_w) .or. .not. mpi_is_root) return
    ag = any(nbase(1:ngroups) == 12)
    if (freeze_v .and. freeze_w) then
      what = 'v and w' ; faces = 'j and k'
      if (ag) what = 'v and w (with AG also P12, P13 and P23)'
    elseif (freeze_v) then
      what = 'v' ; faces = 'j'
      if (ag) what = 'v (with AG also P12 and P23)'
    else
      what = 'w' ; faces = 'k'
      if (ag) what = 'w (with AG also P13 and P23)'
    end if
    if (ndir == 2 .and. delthe /= 0._R8) then
      write(*,'(A)') ' [WARNING] [ICE::setup_freeze] axisymmetric wedge: the swirl '//trim(what)// &
        ' set to zero and not integrated'
    else
      write(*,'(A)') ' [WARNING] [ICE::setup_freeze] '//trim(merge('2-D plane', '1-D      ', ndir == 2))// &
        ' mesh with mirror '//trim(faces)//' faces: '//trim(what)//' set to zero and not integrated'
    end if

  end subroutine setup_freeze


end module ICE_Mod_Freeze
