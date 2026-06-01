module ICE_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: setup_metrics, delthe, meshType

  integer(kind=I4) :: meshType ! 1 -> 1D, 2 -> 2D, 3 -> 3D
  real(kind=R8)    :: delthe   ! grid axisymmetric angle

contains

  subroutine setup_metrics (grid, IOfield)
    use Lib_ORION_data
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    type(orion_data),      intent(inout) :: IOfield
    integer(kind=I4) :: b

    call import_nodes (grid, IOfield)

    call check_mesh_type (grid)

    do b = 1, grid%nb
      call extrapolate_nodes (grid%blk(b)%dim(1), &
                              grid%blk(b)%dim(2), &
                              grid%blk(b)%dim(3), &
                              grid%blk(b)%node)

      call compute_metrics (grid%blk(b))
    end do

    do b = 1, grid%nb
      call compute_metric_tensor (grid%blk(b))
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
    real(kind=R8) :: theta1, theta2, theta(2)

    if (grid%blk(1)%dim(3) > 1) then
      meshType = 3
    elseif (grid%blk(1)%dim(3) == 1 .and. grid%blk(1)%dim(2) == 1) then
      meshType = 1
    else
      meshType = 2
      associate( node => grid%blk(1)%node, jm => grid%blk(1)%dim(2) )
        theta1 = atan2( node(0,1,0)%c(3), node(0,1,0)%c(2) )
        theta2 = atan2( node(0,jm,1)%c(3), node(0,jm,1)%c(2) )
      end associate
      theta(2) = theta2 - theta1
      if ( (theta(1) - theta(2)) < 1.d-5 ) then
        delthe = theta(1)
      else
        delthe = 0.d0
      endif
    endif

    ndir = meshType

  end subroutine check_mesh_type


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
        write(*,*) 'Negative volume in i,j,k', i, j, k
        stop
      endif
      b%vol(i,j,k) = vol
    end do ; end do ; end do
    !$omp end parallel

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


  subroutine compute_metric_tensor (b)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_block_type), intent(inout) :: b
    integer(kind=I4) :: i, j, k, im, jm, km, h
    real(kind=R8)    :: det, A(3,3), cofactor(3,3)

    im = b%dim(1) ; jm = b%dim(2) ; km = b%dim(3)

    !$omp parallel do collapse(3), private(i,j,k,h,A,det,cofactor), shared(b,im,jm,km)
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

    end do
    end do
    end do

  end subroutine compute_metric_tensor


end module ICE_Mod_Metrics
