module ICE_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: setup_metrics, delthe, meshType, ncollapsed_corners

  integer(kind=I4) :: meshType ! 1 -> 1D, 2 -> 2D, 3 -> 3D
  real(kind=R8)    :: delthe = 0._R8  ! angle of the axisymmetric wedge of a 2-D mesh, 0 on any other mesh
  integer(kind=I4) :: ncollapsed_corners = 0 ! cell corners of the last setup_metrics on a collapsed edge, not checked
                                             ! for folding (check_folded_cells)
  real(kind=R8), parameter :: mesh_fit   = 1d-2  ! k node planes match a turn or a translation to this fraction of
                                                 ! the thickness of the block (check_mesh_type)
  real(kind=R8), parameter :: mesh_angle = 1d-6  ! an azimuth (rad) turned by a node, or a distance from the x axis
                                                 ! relative to the largest one, is zero up to this

contains

  subroutine setup_metrics (grid, IOfield, level)
    use Lib_ORION_data
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    type(orion_data),      intent(inout) :: IOfield
    integer(kind=I4), intent(in), optional :: level   !< Multigrid level, for the messages (default 1, the fine one)
    integer(kind=I4) :: b, lev, ncollapsed

    lev = 1
    if (present(level)) lev = level

    call import_nodes (grid, IOfield)

    ncollapsed_corners = 0
    do b = 1, grid%nb
      call check_folded_cells (grid%blk(b)%dim(1), &
                               grid%blk(b)%dim(2), &
                               grid%blk(b)%dim(3), &
                               grid%blk(b)%node, b, lev, ncollapsed)
      ncollapsed_corners = ncollapsed_corners + ncollapsed
    end do

    ! The type of the mesh is that of its fine level: a coarse level, with fewer cells in a
    ! direction, keeps it
    if (lev == 1) call check_mesh_type (grid)

    do b = 1, grid%nb
      call extrapolate_nodes (grid%blk(b)%dim(1), &
                              grid%blk(b)%dim(2), &
                              grid%blk(b)%dim(3), &
                              grid%blk(b)%node)

      call compute_metrics (grid%blk(b))
    end do

    do b = 1, grid%nb
      call compute_metric_tensor (grid%blk(b), b, lev)
    end do

  end subroutine setup_metrics


  subroutine import_nodes (grid, IOfield)
    use Lib_ORION_data
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    type(orion_data),      intent(inout) :: IOfield
    integer(kind=I4) :: b, d

    grid%nb = size(IOfield%block)

    do b = 1, grid%nb
      associate( ni => grid%blk(b)%dim(1), &
                 nj => grid%blk(b)%dim(2), &
                 nk => grid%blk(b)%dim(3)  )
        do d = 1, 3
          grid%blk(b)%node(0:ni, 0:nj, 0:nk)%c(d) = &
          IOfield%block(b)%mesh(d, 0:, 0:, 0:)
        enddo
      end associate
    enddo

  end subroutine import_nodes


  !> The type of the mesh, from the nodes of every block.
  !>  - meshType, and ndir, count the index directions with more than one cell: 3 if a block has
  !>    more than one cell in k, otherwise 2 if a block has more than one cell in j, otherwise 1.
  !>  - delthe is the angle of the axisymmetric wedge of a 2-D mesh, 0 on any other mesh. A 2-D
  !>    mesh is a wedge when, in every block, the k = 1 node plane is the k = 0 plane turned about
  !>    the x axis by the same angle (x and r = sqrt(y^2 + z^2) kept), and plane when, in every
  !>    block, it is the k = 0 plane translated; node_planes checks every node. Blocks of both
  !>    kinds, wedges whose angles differ by more than mesh_fit of the angle, and a block of neither
  !>    kind (nearer one than the other but not to mesh_fit, or as near both) stop the setup.
  !>  - The misfits are lengths, measured against the thickness of the block, not against the size
  !>    of the coordinates: on a mesh written with a finite number of digits, what a turn or a
  !>    translation leaves over is the round-off of the nodes, a fraction of their move from one
  !>    plane to the other. A plane slab far from the axis, which a small turn nearly fits, stays
  !>    plane: the translation fits it better. A node on the axis needs no exception, unlike the
  !>    corners of its collapsed edge in check_folded_cells: its misfit is at most the length of
  !>    its k edge plus its distance from the axis times the angle, a few 1e-10 m where the axis row
  !>    lies at r = 1e-8, against mesh_fit of a thickness of the order of r times the angle.
  subroutine check_mesh_type (grid)
    use ICE_Advanced_Types_m
    use ICE_Global_m, only: ndir
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer(kind=I4) :: b, nturn, first, flat, worst_turn(2), worst_move(2)
    real(kind=R8)    :: sr, ang, ang_b, fit, fit_move
    logical          :: turned, moved, nearer_turn, nearer_move

    meshType = 1
    do b = 1, grid%nb
      if (grid%blk(b)%dim(2) > 1) meshType = max(meshType, 2)
      if (grid%blk(b)%dim(3) > 1) meshType = 3
    end do
    ndir = meshType

    delthe = 0._R8
    if (meshType /= 2) return

    ! The largest distance of a node from the x axis, to tell a probe node on the axis
    sr = 0._R8
    do b = 1, grid%nb
      sr = max(sr, axis_distance(grid%blk(b)%node, grid%blk(b)%dim(1), grid%blk(b)%dim(2)))
    end do

    nturn = 0
    first = 0
    flat  = 0
    ang   = 0._R8
    do b = 1, grid%nb
      call node_planes (grid%blk(b)%node, grid%blk(b)%dim(1), grid%blk(b)%dim(2), sr, turned, ang_b, moved, fit, &
                        fit_move, nearer_turn, nearer_move, worst_turn, worst_move)
      if (nearer_turn .and. .not. turned) then
        write(*,'(A,I0,A,ES12.5,A,ES9.2,A,2(1X,I0),A,ES8.1,A)') ' [ERROR] [ICE::check_mesh_type] block ', b, &
          ': nearer a wedge of ', ang_b, ' rad than a plane block, but its nodes match the turn only to ', fit, &
          ' of its thickness (node', worst_turn, '; ', mesh_fit, ' needed): more significant digits are needed'
        error stop 1
      end if
      if (nearer_move .and. .not. moved) then
        write(*,'(A,I0,A,ES9.2,A,2(1X,I0),A,ES8.1,A)') ' [ERROR] [ICE::check_mesh_type] block ', b, &
          ': nearer a plane block than a wedge, but its k = 1 node plane is the k = 0 plane translated only to ', &
          fit_move, ' of its thickness (node', worst_move, '; ', mesh_fit, ' needed): the k node planes must be '// &
          'parallel, or written with more significant digits'
        error stop 1
      end if
      if (.not. (turned .or. moved)) then
        write(*,'(A,I0,A,ES9.2,A,ES9.2,A)') ' [ERROR] [ICE::check_mesh_type] block ', b, ': its k = 1 node plane '// &
          'is neither the k = 0 plane turned about the x axis nor translated (the turn fits to ', fit, &
          ' of its thickness, the translation to ', fit_move, ')'
        error stop 1
      end if
      if (.not. turned) then
        if (flat == 0) flat = b
        cycle
      end if
      nturn = nturn + 1
      if (nturn == 1) then
        ang   = ang_b
        first = b
      elseif (abs(ang_b - ang) > mesh_fit*abs(ang)) then
        write(*,'(A,I0,A,ES12.5,A,I0,A,ES12.5,A)') ' [ERROR] [ICE::check_mesh_type] block ', b, ' is a wedge of ', &
          ang_b, ' rad, block ', first, ' one of ', ang, ' rad'
        error stop 1
      end if
    end do
    if (nturn == grid%nb) then
      delthe = ang
    elseif (nturn > 0) then
      write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' [ERROR] [ICE::check_mesh_type] ', nturn, ' of the ', grid%nb, &
        ' blocks of the 2-D mesh are axisymmetric wedges (block ', first, '), the others plane (block ', flat, ')'
      error stop 1
    end if

  end subroutine check_mesh_type


  !> The largest distance from the x axis of a node of the planes k = 0 and k = 1 of a 2-D block.
  pure function axis_distance (node, im, jm) result(rmax)
    use ICE_Base_Types_m, only: ICE_vector_3D_type
    implicit none
    type(ICE_vector_3D_type), intent(in) :: node(-1:,-1:,-1:)
    integer(kind=I4),         intent(in) :: im, jm
    real(kind=R8)                        :: rmax
    rmax = maxval(sqrt(node(0:im,0:jm,0:1)%c(2)**2 + node(0:im,0:jm,0:1)%c(3)**2))
  end function axis_distance


  !> Compare the node planes k = 0 and k = 1 of a 2-D block of im x jm cells. Its thickness is the
  !> largest move of a node from one plane to the other (the largest component). turned: every node
  !> of the plane 1 is the node of the plane 0 turned about the x axis by ang, with |ang| >
  !> mesh_angle, to within mesh_fit of the thickness, and closer than as a translation. moved: every
  !> node is the node of the plane 0 translated by the move of the node (0,0), to within mesh_fit of
  !> the thickness, and closer than as a turn or with |ang| <= mesh_angle (a turn by a null angle is a
  !> translation). ang is the azimuth turned by the node (0,jm), or by the node of the plane 0
  !> farthest from the axis when (0,jm) lies on it (no farther than mesh_angle of sr, the largest
  !> distance of a node of the mesh from the axis). fit, fit_move: the misfits of the turn and of the
  !> translation over the thickness, the largest at the nodes worst_turn and worst_move; nearer_turn,
  !> nearer_move: the turn, the translation, fits better.
  !> The nodes come in as a dummy with the block's lower bounds: read through blk%node (or an
  !> associate name for it), ifx 2023.1 at -O2 loaded the nodes (0,0,0) and (0,0,1) as zeros here, so
  !> the move of the node (0,0) was 0 and every plane block was refused, while gfortran and the reads
  !> of the loop below were right. The move of the node (0,0) the loop reads is also checked against
  !> that shift: a build that misreads the nodes stops instead of taking the wrong type.
  subroutine node_planes (node, im, jm, sr, turned, ang, moved, fit, fit_move, nearer_turn, nearer_move, &
                          worst_turn, worst_move)
    use ICE_Base_Types_m, only: ICE_vector_3D_type
    use ICE_Parameters_m, only: pi
    implicit none
    type(ICE_vector_3D_type), intent(in)  :: node(-1:,-1:,-1:)
    integer(kind=I4),         intent(in)  :: im, jm
    real(kind=R8),            intent(in)  :: sr
    logical,                  intent(out) :: turned, moved, nearer_turn, nearer_move
    real(kind=R8),            intent(out) :: ang, fit, fit_move
    integer(kind=I4),         intent(out) :: worst_turn(2), worst_move(2)
    integer(kind=I4) :: i, j, i0, j0
    real(kind=R8)    :: ca, sa, r, rmax, a(3), c(3), shift(3), d, d00, turn_off, move_off, thick

    i0 = 0
    j0 = jm
    if (sqrt(node(0,jm,0)%c(2)**2 + node(0,jm,0)%c(3)**2) <= mesh_angle*sr) then
      rmax = -1._R8
      do j = 0, jm ; do i = 0, im
        r = sqrt(node(i,j,0)%c(2)**2 + node(i,j,0)%c(3)**2)
        if (r > rmax) then
          rmax = r
          i0 = i
          j0 = j
        end if
      end do ; end do
    end if
    ang = atan2(node(i0,j0,1)%c(3), node(i0,j0,1)%c(2)) - atan2(node(i0,j0,0)%c(3), node(i0,j0,0)%c(2))
    if (ang >  pi) ang = ang - 2._R8*pi
    if (ang < -pi) ang = ang + 2._R8*pi

    ca = cos(ang)
    sa = sin(ang)
    shift = node(0,0,1)%c - node(0,0,0)%c
    turn_off = 0._R8
    move_off = 0._R8
    thick    = 0._R8
    d00      = 0._R8
    worst_turn = 0
    worst_move = 0
    do j = 0, jm ; do i = 0, im
      a = node(i,j,0)%c
      c = node(i,j,1)%c
      d = max(abs(c(1) - a(1)), abs(c(2) - (a(2)*ca - a(3)*sa)), abs(c(3) - (a(2)*sa + a(3)*ca)))
      if (d > turn_off) then
        turn_off   = d
        worst_turn = [i, j]
      end if
      d = maxval(abs(c - a - shift))
      if (i == 0 .and. j == 0) d00 = d
      if (d > move_off) then
        move_off   = d
        worst_move = [i, j]
      end if
      thick = max(thick, maxval(abs(c - a)))
    end do ; end do

    if (d00 > mesh_fit*thick) then
      write(*,'(A,ES9.2,A)') ' [ERROR] [ICE::node_planes] the node (0,0) moves from the plane k = 0 to k = 1 by ', &
        d00/thick, ' of the thickness more when read in the loop than when read alone: this build misreads the '// &
        'node array'
      error stop 1
    end if

    nearer_turn = abs(ang) > mesh_angle .and. turn_off < move_off
    nearer_move = move_off < turn_off .or. abs(ang) <= mesh_angle
    turned = nearer_turn .and. turn_off <= mesh_fit*thick
    moved  = nearer_move .and. move_off <= mesh_fit*thick
    fit      = huge(1._R8)
    fit_move = huge(1._R8)
    if (thick > 0._R8) then
      fit      = turn_off/thick
      fit_move = move_off/thick
    end if

  end subroutine node_planes


  !> A cell is the image of the unit cube by the trilinear map of its eight nodes, and a finite
  !> volume needs that map one-to-one: its Jacobian keeps the sign of the block in the whole cell.
  !> At a corner the Jacobian is J = e_i . (e_j x e_k), with e_i, e_j and e_k the three edges from
  !> that corner along i, j and k. The setup stops when, at a corner of a cell of block b,
  !> s J <= 1e-12 |e_i| |e_j| |e_k| (Hadamard's bound, as for the degenerate cells of
  !> compute_metric_tensor), with s the sign of the determinant of the mean edges of the block's
  !> first cell; a degenerate first cell gives the block no orientation and stops it there. A
  !> corner on a collapsed edge, one no longer than 1e-3 of the longest of the four edges of its
  !> cell parallel to it (where the nodes of the cell meet on an axis), is not checked: its
  !> Jacobian says how the nodes of that edge were placed, not how the cell is shaped. Such
  !> corners are counted in ncollapsed. On a 2-D mesh, plane or a wedge whose nodes are off the
  !> axis, the corners decide for the whole cell; in 3-D they are a necessary condition. Every rank
  !> holds every block and checks it, so every rank stops alike; setup_metrics runs on every
  !> multigrid level, and a stop names a coarse one.
  subroutine check_folded_cells (im, jm, km, node, b, level, ncollapsed)
    use ICE_Base_Types_m
    implicit none
    integer(kind=I4),         intent(in)  :: im, jm, km, b, level
    type(ICE_vector_3D_type), intent(in)  :: node(-1:im+1, -1:jm+1, -1:km+1)
    integer(kind=I4),         intent(out) :: ncollapsed
    real(kind=R8), parameter :: collapsed = 1d-3   ! an edge this short against its longest parallel is collapsed
    integer(kind=I4)  :: i, j, k, ca, cb, cc
    real(kind=R8)     :: A(3,3), det, hm, s, h, jac
    real(kind=R8)     :: ei(3,0:1,0:1), ej(3,0:1,0:1), ek(3,0:1,0:1), li(0:1,0:1), lj(0:1,0:1), lk(0:1,0:1)
    logical           :: folded
    character(len=32) :: lev

    lev = ''
    if (level > 1) write(lev,'(A,I0,A)') ' (multigrid level ', level, ')'

    ! Orientation: the determinant of the mean edges of cell 1 1 1, as compute_metric_tensor forms them
    A(1,:) = node(1,1,1)%c - node(0,1,1)%c + node(1,0,1)%c - node(0,0,1)%c + &
             node(1,1,0)%c - node(0,1,0)%c + node(1,0,0)%c - node(0,0,0)%c
    A(2,:) = node(1,1,1)%c - node(1,0,1)%c + node(0,1,1)%c - node(0,0,1)%c + &
             node(1,1,0)%c - node(1,0,0)%c + node(0,1,0)%c - node(0,0,0)%c
    A(3,:) = node(1,1,1)%c - node(1,1,0)%c + node(0,1,1)%c - node(0,1,0)%c + &
             node(1,0,1)%c - node(1,0,0)%c + node(0,0,1)%c - node(0,0,0)%c
    A = 0.25d0 * A
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
    hm = norm2(A(1,:)) * norm2(A(2,:)) * norm2(A(3,:))
    if (.not. (abs(det) > 1d-12*hm)) then
      write(*,'(A,I0,A,ES12.4,A,ES12.4)') ' [ERROR] [ICE::check_folded_cells] block ', b, ', cell 1 1 1'//trim(lev)// &
        ': degenerate, so the block has no orientation; det = ', det, ', |A1| |A2| |A3| = ', hm
      error stop 1
    endif
    s = sign(1.d0, det)

    folded = .false.
    ncollapsed = 0
    !$omp parallel do collapse(3) private(i,j,k,ca,cb,cc,ei,ej,ek,li,lj,lk,h,jac) &
    !$omp shared(node,im,jm,km,s,b,lev,folded) reduction(+:ncollapsed)
    do k = 1, km ; do j = 1, jm ; do i = 1, im
      ! The twelve edges of the cell, e_i(cb,cc) from node (i-1, j-1+cb, k-1+cc), and so on
      do cc = 0, 1 ; do cb = 0, 1
        ei(:,cb,cc) = node(i, j-1+cb, k-1+cc)%c - node(i-1, j-1+cb, k-1+cc)%c
        li(cb,cc)   = norm2(ei(:,cb,cc))
      end do ; end do
      do cc = 0, 1 ; do ca = 0, 1
        ej(:,ca,cc) = node(i-1+ca, j, k-1+cc)%c - node(i-1+ca, j-1, k-1+cc)%c
        lj(ca,cc)   = norm2(ej(:,ca,cc))
      end do ; end do
      do cb = 0, 1 ; do ca = 0, 1
        ek(:,ca,cb) = node(i-1+ca, j-1+cb, k)%c - node(i-1+ca, j-1+cb, k-1)%c
        lk(ca,cb)   = norm2(ek(:,ca,cb))
      end do ; end do
      do cc = 0, 1 ; do cb = 0, 1 ; do ca = 0, 1
        if (li(cb,cc) <= collapsed*maxval(li) .or. lj(ca,cc) <= collapsed*maxval(lj) .or. &
            lk(ca,cb) <= collapsed*maxval(lk)) then
          ncollapsed = ncollapsed + 1
          cycle
        end if
        h   = li(cb,cc) * lj(ca,cc) * lk(ca,cb)
        jac = ei(1,cb,cc)*(ej(2,ca,cc)*ek(3,ca,cb) - ej(3,ca,cc)*ek(2,ca,cb)) &
            + ei(2,cb,cc)*(ej(3,ca,cc)*ek(1,ca,cb) - ej(1,ca,cc)*ek(3,ca,cb)) &
            + ei(3,cb,cc)*(ej(1,ca,cc)*ek(2,ca,cb) - ej(2,ca,cc)*ek(1,ca,cb))
        if (s*jac > 1d-12*h) cycle
        !$omp critical (ICE_metrics_report)
        if (.not. folded) then
          write(*,'(A,I0,A,3(1X,I0),A,3(1X,I0),A)') ' [ERROR] [ICE::check_folded_cells] block ', b, ', cell', i, j, k, &
            trim(lev)//', corner node', i-1+ca, j-1+cb, k-1+cc, ': the cell is folded or turned against its block'
          write(*,'(A,ES12.4,A,ES12.4,A,I0)') '   e_i . (e_j x e_k) = ', jac, ', |e_i| |e_j| |e_k| = ', h, &
            ', sign of the block (from its first cell): ', nint(s)
          write(*,'(A,3(1X,ES14.6))') '   corner node:', node(i-1+ca, j-1+cb, k-1+cc)%c
        endif
        folded = .true.
        !$omp end critical (ICE_metrics_report)
      end do ; end do ; end do
    end do ; end do ; end do
    !$omp end parallel do
    if (folded) then
      write(*,'(A)') ' [ERROR] [ICE::check_folded_cells] folded cell in the mesh'
      error stop 1
    endif

  end subroutine check_folded_cells


  !> Extrapolate ghost cell nodes with 2nd order accuracy.
  subroutine extrapolate_nodes (im, jm, km, node)
    use ICE_Base_Types_m
    implicit none
    integer(kind=I4),         intent(in)    :: im, jm, km
    type(ICE_vector_3D_type), intent(inout) :: node(-1:im+1, -1:jm+1, -1:km+1)
    integer(kind=I4) :: i, j, k

    ! i-faces
    do k = 0, km ; do j = 0, jm
      node(-1,j,k)%c   = 3._R8*node(0,j,k)%c  - 3._R8*node(1,j,k)%c    + node(2,j,k)%c
      node(im+1,j,k)%c = 3._R8*node(im,j,k)%c - 3._R8*node(im-1,j,k)%c + node(im-2,j,k)%c
    enddo ; enddo

    ! j-faces
    if ( jm > 2 ) then
      do k = 0, km ; do i = -1, im+1
        node(i,-1,k)%c   = 3._R8*node(i,0,k)%c  - 3._R8*node(i,1,k)%c    + node(i,2,k)%c
        node(i,jm+1,k)%c = 3._R8*node(i,jm,k)%c - 3._R8*node(i,jm-1,k)%c + node(i,jm-2,k)%c
      enddo ; enddo
    else
      do k = 0, km ; do i = -1, im+1
        node(i,-1,k)%c   = 2._R8*node(i,0,k)%c  - node(i,1,k)%c
        node(i,jm+1,k)%c = 2._R8*node(i,jm,k)%c - node(i,jm-1,k)%c
      end do ; end do
    endif

    ! k-faces: cubic needs two cells, which a block of a 3-D mesh, or its coarse level, may lack
    if ( meshType == 2 .and. delthe == 0d0 .or. meshType == 1 .or. meshType == 3 .and. km < 2 ) then
      do j = -1, jm+1 ; do i = -1, im+1
        node(i,j,-1)%c   = 2._R8*node(i,j,0)%c  - node(i,j,1)%c
        node(i,j,km+1)%c = 2._R8*node(i,j,km)%c - node(i,j,km-1)%c
      end do ; end do
    elseif (meshType == 2 .and. delthe /= 0d0) then
      do j = -1, jm+1 ; do i = -1, im+1
        node(i,j,-1)%c(1)   = node(i,j,0)%c(1)
        node(i,j,km+1)%c(1) = node(i,j,km)%c(1)

        node(i,j,-1)%c(2) = node(i,j,0)%c(2)/cos(delthe*0.5_R8)*cos(delthe*(0.5_R8+1._R8))
        node(i,j,-1)%c(3) = node(i,j,0)%c(3)/sin(delthe*0.5_R8)*sin(delthe*(0.5_R8+1._R8))

        node(i,j,km+1)%c(2) =  node(i,j,-1)%c(2)
        node(i,j,km+1)%c(3) = -node(i,j,-1)%c(3)
      end do ; end do
    elseif (meshType == 3) then
      do j = -1, jm+1 ; do i = -1, im+1
        node(i,j,-1)%c   = 3._R8*node(i,j,0)%c  - 3._R8*node(i,j,1)%c    + node(i,j,2)%c
        node(i,j,km+1)%c = 3._R8*node(i,j,km)%c - 3._R8*node(i,j,km-1)%c + node(i,j,km-2)%c
      end do ; end do
    endif

  end subroutine extrapolate_nodes


  subroutine compute_metrics (b)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_block_type), intent(inout) :: b
    real(kind=R8)    :: Ai, snix, sniy, sniz, Aj, snjx, snjy, snjz, Ak, snkx, snky, snkz, vol
    integer(kind=I4) :: i, j, k, im, jm, km
    real(kind=R8)    :: d1(3), d2(3), d3(3), vx(8), vy(8), vz(8)
    real(kind=R8)    :: snixx, sniyy, snizz, snjxx, snjyy, snjzz, snkxx, snkyy, snkzz
    real(kind=R8)    :: scal, signi, signj, signk
    logical          :: bad_vol

    im = b%dim(1) ; jm = b%dim(2) ; km = b%dim(3)

    ! ------ Sign of normal vectors ------

    ! i direction
    d1 = b%node(1,1,1)%c - b%node(1,0,0)%c
    d2 = b%node(1,0,1)%c - b%node(1,1,0)%c
    d3(1) = d1(2)*d2(3) - d1(3)*d2(2)
    d3(2) = d1(3)*d2(1) - d1(1)*d2(3)
    d3(3) = d1(1)*d2(2) - d1(2)*d2(1)
    snixx = d3(1) ; sniyy = d3(2) ; snizz = d3(3)
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(1,0,1)%c + b%node(1,0,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(0,1,0)%c )
    scal = (d1(1)-d2(1))*snixx + (d1(2)-d2(2))*sniyy + (d1(3)-d2(3))*snizz
    signi = sign(1.d0, scal)

    ! j direction
    d1 = b%node(0,1,1)%c - b%node(1,1,0)%c
    d2 = b%node(1,1,1)%c - b%node(0,1,0)%c
    d3(1) = d1(2)*d2(3) - d1(3)*d2(2)
    d3(2) = d1(3)*d2(1) - d1(1)*d2(3)
    d3(3) = d1(1)*d2(2) - d1(2)*d2(1)
    snjxx = d3(1) ; snjyy = d3(2) ; snjzz = d3(3)
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,1,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(1,0,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
    scal = (d1(1)-d2(1))*snjxx + (d1(2)-d2(2))*snjyy + (d1(3)-d2(3))*snjzz
    signj = sign(1.d0, scal)

    ! k direction
    d1 = b%node(0,1,1)%c - b%node(1,0,1)%c
    d2 = b%node(0,0,1)%c - b%node(1,1,1)%c
    d3(1) = d1(2)*d2(3) - d1(3)*d2(2)
    d3(2) = d1(3)*d2(1) - d1(1)*d2(3)
    d3(3) = d1(1)*d2(2) - d1(2)*d2(1)
    snkxx = d3(1) ; snkyy = d3(2) ; snkzz = d3(3)
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(1,0,1)%c )
    d2 = 0.25d0*( b%node(1,1,0)%c + b%node(0,1,0)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
    scal = (d1(1)-d2(1))*snkxx + (d1(2)-d2(2))*snkyy + (d1(3)-d2(3))*snkzz
    signk = sign(1.d0, scal)

    ! ------ Compute metrics ------

    bad_vol = .false.
    !$omp parallel private (d1,d2,d3,i,j,k,snix,sniy,sniz,Ai,Aj,snjx,snjy,snjz,Ak,snkx,snky,snkz), &
    !$omp private (vx,vy,vz,vol)

    ! i direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 1, jm ; do i = 0, im
      d1 = b%node(i,j,k)%c   - b%node(i,j-1,k-1)%c
      d2 = b%node(i,j-1,k)%c - b%node(i,j,k-1)%c
      d3 = [ d1(2)*d2(3)-d1(3)*d2(2), d1(3)*d2(1)-d1(1)*d2(3), d1(1)*d2(2)-d1(2)*d2(1) ] * 0.5d0
      Ai = sqrt(d3(1)**2 + d3(2)**2 + d3(3)**2)
      if (Ai == 0d0) then
        snix = 0d0 ; sniy = 0d0 ; sniz = 0d0
      else
        snix = d3(1)/Ai*signi ; sniy = d3(2)/Ai*signi ; sniz = d3(3)/Ai*signi
      end if
      b%dir(1)%f(i,j,k)%A = Ai
      b%dir(1)%f(i,j,k)%N = [ snix, sniy, sniz ]
    end do ; end do ; end do

    ! j direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 0, jm ; do i = 1, im
      d1 = b%node(i-1,j,k)%c - b%node(i,j,k-1)%c
      d2 = b%node(i,j,k)%c   - b%node(i-1,j,k-1)%c
      d3 = [ d1(2)*d2(3)-d1(3)*d2(2), d1(3)*d2(1)-d1(1)*d2(3), d1(1)*d2(2)-d1(2)*d2(1) ] * 0.5d0
      Aj = sqrt(d3(1)**2 + d3(2)**2 + d3(3)**2)
      if (Aj == 0d0) then
        snjx = 0d0 ; snjy = 0d0 ; snjz = 0d0
      else
        snjx = d3(1)/Aj*signj ; snjy = d3(2)/Aj*signj ; snjz = d3(3)/Aj*signj
      end if
      b%dir(2)%f(i,j,k)%A = Aj
      b%dir(2)%f(i,j,k)%N = [ snjx, snjy, snjz ]
    end do ; end do ; end do

    ! k direction
    !$omp do collapse(3)
    do k = 0, km ; do j = 1, jm ; do i = 1, im
      d1 = b%node(i-1,j,k)%c   - b%node(i,j-1,k)%c
      d2 = b%node(i-1,j-1,k)%c - b%node(i,j,k)%c
      d3 = [ d1(2)*d2(3)-d1(3)*d2(2), d1(3)*d2(1)-d1(1)*d2(3), d1(1)*d2(2)-d1(2)*d2(1) ] * 0.5d0
      Ak = sqrt(d3(1)**2 + d3(2)**2 + d3(3)**2)
      if (Ak == 0d0) then
        snkx = 0d0 ; snky = 0d0 ; snkz = 0d0
      else
        snkx = d3(1)/Ak*signk ; snky = d3(2)/Ak*signk ; snkz = d3(3)/Ak*signk
      end if
      b%dir(3)%f(i,j,k)%A = Ak
      b%dir(3)%f(i,j,k)%N = [ snkx, snky, snkz ]
    end do ; end do ; end do

    ! Cell volumes
    !$omp do collapse(3)
    do k = 1, km ; do j = 1, jm ; do i = 1, im
      vx = [ b%node(i-1,j-1,k-1)%c(1), b%node(i,j-1,k-1)%c(1), b%node(i-1,j,k-1)%c(1), b%node(i,j,k-1)%c(1), &
             b%node(i-1,j-1,k  )%c(1), b%node(i,j-1,k  )%c(1), b%node(i-1,j,k  )%c(1), b%node(i,j,k  )%c(1) ]
      vy = [ b%node(i-1,j-1,k-1)%c(2), b%node(i,j-1,k-1)%c(2), b%node(i-1,j,k-1)%c(2), b%node(i,j,k-1)%c(2), &
             b%node(i-1,j-1,k  )%c(2), b%node(i,j-1,k  )%c(2), b%node(i-1,j,k  )%c(2), b%node(i,j,k  )%c(2) ]
      vz = [ b%node(i-1,j-1,k-1)%c(3), b%node(i,j-1,k-1)%c(3), b%node(i-1,j,k-1)%c(3), b%node(i,j,k-1)%c(3), &
             b%node(i-1,j-1,k  )%c(3), b%node(i,j-1,k  )%c(3), b%node(i-1,j,k  )%c(3), b%node(i,j,k  )%c(3) ]
      vol = tvol(vx,vy,vz,1,2,3,5) + tvol(vx,vy,vz,2,4,3,8) &
          + tvol(vx,vy,vz,5,8,6,2) + tvol(vx,vy,vz,5,7,8,3) &
          + tvol(vx,vy,vz,5,8,2,3)
      if (vol <= 0d0) then
        !$omp critical (ICE_metrics_report)
        write(*,*) 'Negative volume in i,j,k', i, j, k
        bad_vol = .true.
        !$omp end critical (ICE_metrics_report)
      endif
      b%vol(i,j,k) = vol
    end do ; end do ; end do
    !$omp end parallel
    if (bad_vol) then
      write(*,'(A)') ' [ERROR] [ICE::compute_metrics] non-positive cell volume'
      error stop 1
    endif

  end subroutine compute_metrics


  pure function tvol(vx,vy,vz,ind1,ind2,ind3,ind4) result(volume)
    implicit none
    real(kind=R8),    intent(in) :: vx(8), vy(8), vz(8)
    integer(kind=I4), intent(in) :: ind1, ind2, ind3, ind4
    real(kind=R8) :: volume

    volume = abs( ( (vx(ind2)-vx(ind1)) * &
        ( (vy(ind3)-vy(ind1))*(vz(ind4)-vz(ind1)) - (vy(ind4)-vy(ind1))*(vz(ind3)-vz(ind1)) ) + &
                    (vy(ind2)-vy(ind1)) * &
        ( (vx(ind4)-vx(ind1))*(vz(ind3)-vz(ind1)) - (vx(ind3)-vx(ind1))*(vz(ind4)-vz(ind1)) ) + &
                    (vz(ind2)-vz(ind1)) * &
        ( (vx(ind3)-vx(ind1))*(vy(ind4)-vy(ind1)) - (vx(ind4)-vx(ind1))*(vy(ind3)-vy(ind1)) ) ) &
        / 6.d0 )

  end function tvol


  subroutine compute_metric_tensor (b, ib, level)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_block_type), intent(inout) :: b
    integer(kind=I4),     intent(in)    :: ib, level   !< Block number and multigrid level, for the messages
    integer(kind=I4) :: i, j, k, im, jm, km, h
    real(kind=R8)    :: det, A(3,3), cofactor(3,3)
    logical          :: degenerate
    character(len=32) :: lev

    im = b%dim(1) ; jm = b%dim(2) ; km = b%dim(3)

    lev = ''
    if (level > 1) write(lev,'(A,I0,A)') ' (multigrid level ', level, ')'
    degenerate = .false.
    !$omp parallel do collapse(3), private(i,j,k,h,A,det,cofactor), shared(b,im,jm,km,ib,lev,degenerate)
    do k = 0, km+1
    do j = 0, jm+1
    do i = 0, im+1

      A(1,:) = b%node(i  ,j  ,k  )%c - b%node(i-1,j  ,k  )%c + &
               b%node(i  ,j-1,k  )%c - b%node(i-1,j-1,k  )%c + &
               b%node(i  ,j  ,k-1)%c - b%node(i-1,j  ,k-1)%c + &
               b%node(i  ,j-1,k-1)%c - b%node(i-1,j-1,k-1)%c

      A(2,:) = b%node(i  ,j  ,k  )%c - b%node(i  ,j-1,k  )%c + &
               b%node(i-1,j  ,k  )%c - b%node(i-1,j-1,k  )%c + &
               b%node(i  ,j  ,k-1)%c - b%node(i  ,j-1,k-1)%c + &
               b%node(i-1,j  ,k-1)%c - b%node(i-1,j-1,k-1)%c

      A(3,:) = b%node(i  ,j  ,k  )%c - b%node(i  ,j  ,k-1)%c + &
               b%node(i-1,j  ,k  )%c - b%node(i-1,j  ,k-1)%c + &
               b%node(i  ,j-1,k  )%c - b%node(i  ,j-1,k-1)%c + &
               b%node(i-1,j-1,k  )%c - b%node(i-1,j-1,k-1)%c

      A = 0.25d0 * A

      det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

      if ( abs(det) == 0d0 ) then
        b%M(i,j,k)%c = 0d0
        do h = 1, 3
          b%M(i,j,k)%c(h,h) = 1d0
        end do
      else
        cofactor(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        cofactor(1,2) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
        cofactor(1,3) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        cofactor(2,1) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
        cofactor(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        cofactor(2,3) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
        cofactor(3,1) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        cofactor(3,2) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
        cofactor(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
        b%M(i,j,k)%c = cofactor / det
      end if

      do h = 1, 3
        b%dl(i,j,k)%c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
      end do

      !> An interior cell has a volume: |det| above 1e-12 of the product of the lengths of its
      !> mean edges (Hadamard's bound). A ghost cell is extrapolated and only its length across
      !> the face is read, so its metric keeps the identity when det is zero.
      if (i >= 1 .and. i <= im .and. j >= 1 .and. j <= jm .and. k >= 1 .and. k <= km) then
        if (.not. (abs(det) > 1d-12*b%dl(i,j,k)%c(1)*b%dl(i,j,k)%c(2)*b%dl(i,j,k)%c(3))) then
          !$omp critical (ICE_metrics_report)
          if (.not. degenerate) &
            write(*,'(A,I0,A,3(1X,I0),A,ES12.4,A,ES12.4)') ' [ERROR] [ICE::compute_metric_tensor] block ', ib, &
              ', cell', i, j, k, trim(lev)//': degenerate; det = ', det, &
              ', |A1| |A2| |A3| = ', b%dl(i,j,k)%c(1)*b%dl(i,j,k)%c(2)*b%dl(i,j,k)%c(3)
          degenerate = .true.
          !$omp end critical (ICE_metrics_report)
        end if
      end if

    end do
    end do
    end do
    if (degenerate) then
      write(*,'(A)') ' [ERROR] [ICE::compute_metric_tensor] degenerate cell in the mesh'
      error stop 1
    endif

  end subroutine compute_metric_tensor


end module ICE_Mod_Metrics
