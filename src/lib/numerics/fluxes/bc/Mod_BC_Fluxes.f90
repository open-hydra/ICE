module ICE_Mod_BC_Fluxes
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Lib_Ghost, only: compute_ghost, fill_second_ghost

  implicit none
  private
  public :: compute_ghost, fill_second_ghost, compute_bound

contains

  subroutine compute_bound (grid)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Lib_Riemann
    use ICE_Mod_Fluxes, only: state_reconstruction
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: n, b, f, p, i, j, k
    integer(kind=I4) :: ig, jg, kg, ig2, jg2, kg2, ip, jp, kp
    integer(kind=I4) :: dir
    real(kind=R8)    :: normal(3), area
    real(kind=R8)    :: dl0, dl1, dl2, dll, dlr, dl_g1, dl_m, dl_4th
    real(kind=R8)    :: beta_val
    real(kind=R8)    :: priml(12), primr(12)

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(grid, ngroups, ncond, riemann), &
    !$OMP PRIVATE(n, b, f, p, i, j, k, ig, jg, kg, ig2, jg2, kg2, ip, jp, kp, &
    !$OMP         dir, normal, area, dl0, dl1, dl2, dll, dlr, dl_g1, dl_m, dl_4th, &
    !$OMP         beta_val, priml, primr)
    do b = 1, grid%nb
      do p = 1, ngroups

        !$OMP DO COLLAPSE (3)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)

          grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) = 0._R8

        enddo ; enddo ; enddo
        !$OMP END DO
      enddo
    enddo

    !$OMP DO SCHEDULE (DYNAMIC)
    do n = 1, size(grid%bc)

      if (grid%bc(n)%type == 0 .or. grid%bc(n)%type == 2) cycle

      b = grid%bc(n)%b
      i = grid%bc(n)%i ; j = grid%bc(n)%j ; k = grid%bc(n)%k
      p = grid%bc(n)%p ; f = grid%bc(n)%f

      !> Ghost and interior neighbor coordinates
      ig  = i -   guide(f,1) ; jg  = j -   guide(f,2) ; kg  = k -   guide(f,3)
      ig2 = i - 2*guide(f,1) ; jg2 = j - 2*guide(f,2) ; kg2 = k - 2*guide(f,3)
      ip  = i +   guide(f,1) ; jp  = j +   guide(f,2) ; kp  = k +   guide(f,3)

      dir      = (f+1)/2
      beta_val = grid%blk(b)%cond_phase(p)%beta(i,j,k)

      select case (f)

        case (1,3,5)
        !> Odd faces: stencil (g2, g1, m, m+1) → priml=ghost side, primr=interior side
        normal = grid%blk(b)%dir(dir)%f(ig,jg,kg)%N
        area   = grid%blk(b)%dir(dir)%f(ig,jg,kg)%A

        dl_g1 = grid%blk(b)%dl(ig,jg,kg)%c(dir)     ! dl_g2 ≈ dl_g1 (out of bounds)
        dl_m  = grid%blk(b)%dl(i,j,k)%c(dir)
        dl_4th= grid%blk(b)%dl(ip,jp,kp)%c(dir)
        dl0 = dl_g1 ;  dl1 = 0.5_R8*(dl_g1+dl_m) ;  dl2 = 0.5_R8*(dl_m+dl_4th)
        dll = 0.5_R8*dl_g1 ;  dlr = 0.5_R8*dl_m

        call state_reconstruction(grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig2,jg2,kg2), &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig,jg,kg),    &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),i,j,k),       &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ip,jp,kp),    &
                                  dl0, dl1, dl2, dll, dlr,                                 &
                                  priml(1:ncond(p)), primr(1:ncond(p)), beta_val)

        grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) = &
          grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) + &
          riemann(priml(1:ncond(p)), primr(1:ncond(p)), normal) * area

        case (2,4,6)
        !> Even faces: stencil (m-1, m, g1, g2) → priml=interior side, primr=ghost side
        normal = grid%blk(b)%dir(dir)%f(i,j,k)%N
        area   = grid%blk(b)%dir(dir)%f(i,j,k)%A

        dl_4th= grid%blk(b)%dl(ip,jp,kp)%c(dir)   ! ip = m-1 for even faces
        dl_m  = grid%blk(b)%dl(i,j,k)%c(dir)
        dl_g1 = grid%blk(b)%dl(ig,jg,kg)%c(dir)   ! dl_g2 ≈ dl_g1 (out of bounds)
        dl0 = 0.5_R8*(dl_4th+dl_m) ;  dl1 = 0.5_R8*(dl_m+dl_g1) ;  dl2 = dl_g1
        dll = 0.5_R8*dl_m ;  dlr = 0.5_R8*dl_g1

        call state_reconstruction(grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ip,jp,kp),    &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),i,j,k),       &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig,jg,kg),    &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig2,jg2,kg2), &
                                  dl0, dl1, dl2, dll, dlr,                                 &
                                  priml(1:ncond(p)), primr(1:ncond(p)), beta_val)

        grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) = &
          grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) - &
          riemann(priml(1:ncond(p)), primr(1:ncond(p)), normal) * area

      end select

    enddo
    !$OMP END PARALLEL

  end subroutine compute_bound


end module ICE_Mod_BC_Fluxes
