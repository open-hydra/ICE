module ICE_Lib_Multigrid
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Check_Multigrid, Coarse_Grid, Coarse_IOfield
  public :: fine2coarse_prim, coarse2fine_prim

contains


  subroutine Check_Multigrid(domain)
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_multigrid
    use ir_precision, only: str
    implicit none
    type(ICE_domain_type), intent(in) :: domain
    integer :: b, d, rap, check

    rap = 2**(obj_multigrid%MGL - 1)
    do b = 1, domain%nb
      check = 0
      do d = 1, 3
        if (mod(domain%blk(b)%dim(d), rap) == 0) check = check + 1
        if (d == 3) then
          if (domain%blk(b)%dim(d) == 1) check = check + 1
        end if
      end do
      if (check < 3) then
        write(*,'(A)') ' [ERROR] [ICE::Check_Multigrid] block '//trim(str(.true.,b))//' not divisible by 2^(MGL-1)'
        error stop 1
      end if
    end do

  end subroutine Check_Multigrid


  ! Volume-weighted restriction: fine prim -> coarse prim, for one particle group.
  ! Works in prim space for non-conserved extras, converts to cons for conserved vars.
  subroutine fine2coarse_prim(p, fPrim, cPrim, fVol, cVol, fDim, cDim)
    use ICE_Global_m,    only: ncond
    use ICE_Lib_Model,   only: prim_2_cons, cons_2_prim, assign_prim_2_cons, assign_cons_2_prim
    implicit none
    integer(I4), intent(in) :: p
    integer(I4), intent(in) :: fDim(3), cDim(3)
    real(R8), intent(in)  :: fPrim(1:ncond(p), 0:fDim(1)+1, 0:fDim(2)+1, 0:fDim(3)+1)
    real(R8), intent(out) :: cPrim(1:ncond(p), 0:cDim(1)+1, 0:cDim(2)+1, 0:cDim(3)+1)
    real(R8), intent(in)  :: fVol(1:fDim(1), 1:fDim(2), 1:fDim(3))
    real(R8), intent(in)  :: cVol(1:cDim(1), 1:cDim(2), 1:cDim(3))
    ! Local
    integer  :: i, j, k, i2, j2, k2, i2d, j2d, k2d
    real(R8) :: fCons(1:ncond(p), 1:fDim(1), 1:fDim(2), 1:fDim(3))
    real(R8) :: cCons(1:ncond(p), 1:cDim(1), 1:cDim(2), 1:cDim(3))

    call assign_prim_2_cons(p)
    call assign_cons_2_prim(p)

    !$omp do collapse(2)
    do k = 1, fDim(3)
    do j = 1, fDim(2)
    do i = 1, fDim(1)
      fCons(:,i,j,k) = prim_2_cons(fPrim(:,i,j,k))
    end do ; end do ; end do

    !$omp do collapse(2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      i2  = 2*i  ; j2  = 2*j  ; k2  = 2*k
      i2d = i2-1 ; j2d = j2-1 ; k2d = k2-1

      if (fDim(3) == 1) then
        cCons(:,i,j,k) = fCons(:,i2d,j2d,1)*fVol(i2d,j2d,1) &
                       + fCons(:,i2, j2d,1)*fVol(i2, j2d,1) &
                       + fCons(:,i2d,j2, 1)*fVol(i2d,j2, 1) &
                       + fCons(:,i2, j2, 1)*fVol(i2, j2, 1)
      else
        cCons(:,i,j,k) = fCons(:,i2d,j2d,k2d)*fVol(i2d,j2d,k2d) &
                       + fCons(:,i2, j2d,k2d)*fVol(i2, j2d,k2d) &
                       + fCons(:,i2d,j2, k2d)*fVol(i2d,j2, k2d) &
                       + fCons(:,i2, j2, k2d)*fVol(i2, j2, k2d) &
                       + fCons(:,i2d,j2d,k2 )*fVol(i2d,j2d,k2 ) &
                       + fCons(:,i2, j2d,k2 )*fVol(i2, j2d,k2 ) &
                       + fCons(:,i2d,j2, k2 )*fVol(i2d,j2, k2 ) &
                       + fCons(:,i2, j2, k2 )*fVol(i2, j2, k2 )
      end if
    end do ; end do ; end do

    !$omp do collapse(2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      cCons(:,i,j,k) = cCons(:,i,j,k) / cVol(i,j,k)
      cPrim(:,i,j,k) = cons_2_prim(cCons(:,i,j,k))
    end do ; end do ; end do

  end subroutine fine2coarse_prim


  ! Trilinear prolongation: coarse prim -> fine prim, for one particle group.
  subroutine coarse2fine_prim(p, fPrim, cPrim, fDim, cDim)
    use ICE_Global_m, only: ncond
    implicit none
    integer(I4), intent(in) :: p
    integer(I4), intent(in) :: fDim(3), cDim(3)
    real(R8), intent(out) :: fPrim(1:ncond(p), 0:fDim(1)+1, 0:fDim(2)+1, 0:fDim(3)+1)
    real(R8), intent(in)  :: cPrim(1:ncond(p), 0:cDim(1)+1, 0:cDim(2)+1, 0:cDim(3)+1)
    ! Local
    integer :: i, j, k, ii, jj, kk, counter
    integer :: i2, j2, k2, i2d, j2d, k2d, im, jm, km, ip, jp, kp
    integer :: mask(3), id(6)
    real(R8) :: a1, a2, a3, a4, coeffs(8), interp(ncond(p))

    a1 = 27.d0/64.d0 ; a2 = 9.d0/64.d0 ; a3 = 3.d0/64.d0 ; a4 = 1.d0/64.d0
    coeffs(1:8) = [ a1, a2, a2, a2, a3, a3, a3, a4 ]

    if (fDim(3) == 1) then
      a1 = 9.d0/16.d0 ; a2 = 3.d0/16.d0 ; a3 = 1.d0/16.d0 ; a4 = 0.d0
      coeffs(1:8) = [ a1, a2, a2, a4, a3, a4, a4, a4 ]
    end if

    !$omp do collapse(2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2  = 2*i  ; j2  = 2*j  ; k2  = 2*k
      i2d = i2-1 ; j2d = j2-1 ; k2d = k2-1

      im = max(1, i-1) ; ip = min(cDim(1), i+1)
      jm = max(1, j-1) ; jp = min(cDim(2), j+1)
      km = max(1, k-1) ; kp = min(cDim(3), k+1)

      if (fDim(3) == 1) then
        k2  = 1 ; k2d = 1
      end if

      id(1:6) = [ im, jm, km, ip, jp, kp ]

      counter = 1
      do kk = k2d, k2
      do jj = j2d, j2
      do ii = i2d, i2

        if (fDim(3) == 1) then
          if (counter==1) mask = [1,2,3]
          if (counter==2) mask = [4,2,3]
          if (counter==3) mask = [1,5,3]
          if (counter==4) mask = [4,5,3]
        else
          if (counter==1) mask = [1,2,3]
          if (counter==2) mask = [4,2,3]
          if (counter==3) mask = [1,5,3]
          if (counter==4) mask = [4,5,3]
          if (counter==5) mask = [1,2,6]
          if (counter==6) mask = [4,2,6]
          if (counter==7) mask = [1,5,6]
          if (counter==8) mask = [4,5,6]
        end if

        interp = coeffs(1)*cPrim(:,i,j,k)                              &
               + coeffs(2)*cPrim(:,id(mask(1)),j,k)                    &
               + coeffs(3)*cPrim(:,i,id(mask(2)),k)                    &
               + coeffs(4)*cPrim(:,i,j,id(mask(3)))                    &
               + coeffs(5)*cPrim(:,id(mask(1)),id(mask(2)),k)          &
               + coeffs(6)*cPrim(:,id(mask(1)),j,id(mask(3)))          &
               + coeffs(7)*cPrim(:,i,id(mask(2)),id(mask(3)))          &
               + coeffs(8)*cPrim(:,id(mask(1)),id(mask(2)),id(mask(3)))

        fPrim(:,ii,jj,kk) = interp
        counter = counter + 1

      end do ; end do ; end do

    end do ; end do ; end do

  end subroutine coarse2fine_prim


  subroutine Coarse_Grid(Fine, Coarse)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(in)    :: Fine
    type(ICE_domain_type), intent(inout) :: Coarse
    integer :: b, i, j, k, i2, j2, k2

    do b = 1, Coarse%nb
      !$omp do collapse(2)
      do k = 0, Fine%blk(b)%dim(3), 2 - mod(Fine%blk(b)%dim(3), 2)
      do j = 0, Fine%blk(b)%dim(2), 2
      do i = 0, Fine%blk(b)%dim(1), 2
        i2 = i/2 ; j2 = j/2 ; k2 = k/2
        if (Fine%blk(b)%dim(3) == 1) k2 = k
        Coarse%blk(b)%node(i2,j2,k2)%c = Fine%blk(b)%node(i,j,k)%c
      end do ; end do ; end do
    end do

  end subroutine Coarse_Grid


  subroutine Coarse_IOfield(fIOfield, cIOfield)
    use Lib_ORION_data
    implicit none
    type(ORION_data), intent(in)    :: fIOfield
    type(ORION_data), intent(inout) :: cIOfield
    integer :: b, nb, Ni, Nj, Nk, Onvar

    Onvar = size(fIOfield%block(1)%vars, 1)
    nb    = size(fIOfield%block)
    allocate(cIOfield%block(nb))

    do b = 1, nb
      Ni = fIOfield%block(b)%Ni / 2
      Nj = fIOfield%block(b)%Nj / 2
      Nk = fIOfield%block(b)%Nk / 2
      Nk = max(1, Nk)
      cIOfield%block(b)%Ni = Ni
      cIOfield%block(b)%Nj = Nj
      cIOfield%block(b)%Nk = Nk
      allocate( cIOfield%block(b)%mesh(1:3, 0:Ni, 0:Nj, 0:Nk) )
      allocate( cIOfield%block(b)%vars(1:Onvar, Ni, Nj, Nk) )
    end do

  end subroutine Coarse_IOfield


end module ICE_Lib_Multigrid
