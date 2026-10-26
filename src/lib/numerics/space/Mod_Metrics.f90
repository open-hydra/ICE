module ICE_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: setup_metrics, delthe, meshType, ncollapsed_corners

  integer(kind=I4) :: meshType ! 1 -> 1D, 2 -> 2D, 3 -> 3D
  real(kind=R8)    :: delthe   ! grid axisymmetric angle
  integer(kind=I4) :: ncollapsed_corners = 0 ! cell corners of the last setup_metrics on a collapsed edge, not checked
                                             ! for folding (check_folded_cells)

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

    call check_mesh_type (grid)

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


  subroutine check_mesh_type (grid)
    use ICE_Advanced_Types_m
    use ICE_Global_m, only: ndir
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    real(kind=R8) :: theta(2)

    if (grid%blk(1)%dim(3) > 1) then
      meshType = 3
    elseif (grid%blk(1)%dim(3) == 1 .and. grid%blk(1)%dim(2) == 1) then
      meshType = 1
    else
      meshType = 2
      associate( node => grid%blk(1)%node, jm => grid%blk(1)%dim(2) )
        theta(1) = atan2( node(0,jm,1)%c(3), node(0,jm,1)%c(2) ) - atan2( node(0,jm,0)%c(3), node(0,jm,0)%c(2) )
        theta(2) = atan2( node(0,1,1)%c(3),  node(0,1,1)%c(2) )  - atan2( node(0,1,0)%c(3),  node(0,1,0)%c(2) )
      end associate
      if ( abs(theta(1) - theta(2)) < 1.d-5 .and. theta(1) /= 0.d0 ) then
        delthe = theta(1)
      else
        delthe = 0.d0
      endif
    endif

    ndir = meshType

  end subroutine check_mesh_type


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

    ! k-faces
    if ( meshType == 2 .and. delthe == 0d0 .or. meshType == 1 ) then
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
