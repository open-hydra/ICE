module ICE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Lib_Reconstruction, only: state_reconstruction

  implicit none
  private
  public :: compute_flux, state_reconstruction

contains

  subroutine compute_flux (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m,     only: obj_space_scheme
    use ICE_Lib_Shock_Detector, only: SD_rho
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, pass

    !$OMP PARALLEL
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !> Shock detector: beta=1 (smooth) or beta=0 (shock)
      if (obj_space_scheme%SD) then
        !$OMP DO COLLAPSE (3)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)
          grid%blk(b)%cond_phase(p)%beta(i,j,k) = &
            SD_rho(grid%blk(b)%cond_phase(p)%prim(1, i-1:i+1, j-1:j+1, k-1:k+1))
        enddo ; enddo ; enddo
        !$OMP END DO
      else
        !$OMP DO COLLAPSE (3)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)
          grid%blk(b)%cond_phase(p)%beta(i,j,k) = 1d0
        enddo ; enddo ; enddo
        !$OMP END DO
      end if

      !> Face i adds to cells i and i+1, so two faces sharing a cell must not run
      !> concurrently: sweep odd faces, then even faces (the END DO barrier separates them).
      do pass = 1, 2

      !$OMP DO COLLAPSE (3)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = pass, grid%blk(b)%dim(1)-1, 2

        call compute_flux_ (grid%blk(b)%cond_phase(p)%prim(1:ncond(p),i-1:i+2,j,k),    &
                            [grid%blk(b)%dl(i-1,j,k)%c(1), grid%blk(b)%dl(i,j,k)%c(1), &
                             grid%blk(b)%dl(i+1,j,k)%c(1), grid%blk(b)%dl(i+2,j,k)%c(1)], &
                            grid%blk(b)%dir(1)%f(i,j,k)%N,                     &
                            grid%blk(b)%dir(1)%f(i,j,k)%A,                     &
                            grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i:i+1,j,k),  &
                            ncond(p),                                                    &
                            grid%blk(b)%cond_phase(p)%beta(i,j,k)                      )

      end do ; end do ; end do
      !$OMP END DO

      enddo

      do pass = 1, 2

      !$OMP DO COLLAPSE (3)
      do k = 1, grid%blk(b)%dim(3)
      do i = 1, grid%blk(b)%dim(1)
      do j = pass, grid%blk(b)%dim(2)-1, 2

        call compute_flux_ (grid%blk(b)%cond_phase(p)%prim(1:ncond(p),i,j-1:j+2,k),    &
                            [grid%blk(b)%dl(i,j-1,k)%c(2), grid%blk(b)%dl(i,j,k)%c(2), &
                             grid%blk(b)%dl(i,j+1,k)%c(2), grid%blk(b)%dl(i,j+2,k)%c(2)], &
                            grid%blk(b)%dir(2)%f(i,j,k)%N,                     &
                            grid%blk(b)%dir(2)%f(i,j,k)%A,                     &
                            grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j:j+1,k),  &
                            ncond(p),                                                    &
                            grid%blk(b)%cond_phase(p)%beta(i,j,k)                      )

      enddo ; enddo ; enddo
      !$OMP END DO

      enddo

      do pass = 1, 2

      !$OMP DO COLLAPSE (3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)
      do k = pass, grid%blk(b)%dim(3)-1, 2

        call compute_flux_ (grid%blk(b)%cond_phase(p)%prim(1:ncond(p),i,j,k-1:k+2),    &
                            [grid%blk(b)%dl(i,j,k-1)%c(3), grid%blk(b)%dl(i,j,k)%c(3), &
                             grid%blk(b)%dl(i,j,k+1)%c(3), grid%blk(b)%dl(i,j,k+2)%c(3)], &
                            grid%blk(b)%dir(3)%f(i,j,k)%N,                     &
                            grid%blk(b)%dir(3)%f(i,j,k)%A,                     &
                            grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k:k+1),  &
                            ncond(p),                                                    &
                            grid%blk(b)%cond_phase(p)%beta(i,j,k)                      )

      enddo ; enddo ; enddo
      !$OMP END DO

      enddo

    enddo
    !$OMP END PARALLEL

  end subroutine compute_flux
            

  subroutine compute_flux_ (prim, length, normal, area, residual, nvar, beta)
    use ICE_Lib_Riemann
    implicit none
    real(kind=R8), dimension(1:nvar,-1:2), intent(in)    :: prim
    real(kind=R8), dimension(-1:2),        intent(in)    :: length
    real(kind=R8), dimension(3),           intent(in)    :: normal
    real(kind=R8),                         intent(in)    :: area
    real(kind=R8), dimension(1:nvar,0:1),  intent(inout) :: residual
    integer(kind=I4),                      intent(in)    :: nvar
    real(kind=R8),                         intent(in)    :: beta

    real(kind=R8) :: dl0, dl1, dl2, dll, dlr
    real(kind=R8), dimension(size(prim(:,1))) :: prim_1, prim_4, flux

    dl0 = 0.5_R8 * (length(-1)+length(0))
    dl1 = 0.5_R8 * (length(0)+length(1))
    dl2 = 0.5_R8 * (length(1)+length(2))
    dll = 0.5d0 * length(0)
    dlr = 0.5d0 * length(1)

    call state_reconstruction (prim(:,-1), prim(:,0), prim(:,1), prim(:,2), &
                               dl0, dl1, dl2, dll, dlr, prim_1, prim_4, beta)

    flux = riemann (prim_1, prim_4, normal)
    flux = flux * area

    residual (:,0) = residual (:,0) - flux
    residual (:,1) = residual (:,1) + flux

  end subroutine compute_flux_



end module ICE_Mod_Fluxes