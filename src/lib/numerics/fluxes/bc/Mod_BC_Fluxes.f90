module ICE_Mod_BC_Fluxes
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Lib_Ghost, only: compute_ghost, fill_second_ghost

  implicit none
  private
  public :: compute_ghost, fill_second_ghost, compute_bound

contains

  subroutine compute_bound (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Lib_Riemann
    use ICE_Mod_Fluxes, only: state_reconstruction, bad_recon
    use ICE_Config_Types_m, only: obj_time_scheme, obj_condensed
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4),      intent(in)    :: p
    integer(kind=I4) :: n, nn, b, f, i, j, k
    integer(kind=I4) :: ig, jg, kg, ig2, jg2, kg2, ip, jp, kp
    integer(kind=I4) :: dir
    real(kind=R8)    :: normal(3), area
    real(kind=R8)    :: dl0, dl1, dl2, dll, dlr, dl_g1, dl_m, dl_4th
    real(kind=R8)    :: beta_val
    real(kind=R8)    :: priml(ncond_max), primr(ncond_max), flux(ncond_max)
    integer(kind=I4) :: v

    bad_recon = .false.
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(grid, ngroups, ncond, riemann), &
    !$OMP PRIVATE(n, nn, b, f, i, j, k, ig, jg, kg, ig2, jg2, kg2, ip, jp, kp, &
    !$OMP         dir, normal, area, dl0, dl1, dl2, dll, dlr, dl_g1, dl_m, dl_4th, &
    !$OMP         beta_val, priml, primr, flux, v)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$OMP DO COLLAPSE (3)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) = 0._R8

      enddo ; enddo ; enddo
      !$OMP END DO
    enddo

    !> This rank's own boundary entries only; see build_local_bc_index.
    !$OMP DO SCHEDULE (DYNAMIC, 64)
    do nn = 1, grid%n_local_bc
      n = grid%local_bc_idx(nn)

      if (grid%bc(n)%type == 0) cycle

      b = grid%bc(n)%b
      i = grid%bc(n)%i ; j = grid%bc(n)%j ; k = grid%bc(n)%k
      if (grid%bc(n)%p /= p) cycle
      f = grid%bc(n)%f

      !> Nothing crosses a symmetry face of a pressureless cloud
      if (grid%bc(n)%type == 200 .and. trim(obj_time_scheme%model(p)) == 'MK') cycle

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
                                  priml(1:ncond(p)), primr(1:ncond(p)), beta_val,        &
                                  obj_condensed(mat_of(p)), solid_of(p))

        flux(1:ncond(p)) = riemann(priml(1:ncond(p)), primr(1:ncond(p)), normal, obj_condensed(mat_of(p))) * area

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
                                  priml(1:ncond(p)), primr(1:ncond(p)), beta_val,        &
                                  obj_condensed(mat_of(p)), solid_of(p))

        flux(1:ncond(p)) = - riemann(priml(1:ncond(p)), primr(1:ncond(p)), normal, obj_condensed(mat_of(p))) * area

      end select

      !> Edge and corner cells carry one entry per boundary face, handled by different
      !> threads: accumulate atomically.
      do v = 1, ncond(p)
        !$OMP ATOMIC
        grid%blk(b)%cond_phase(p)%residual(v,i,j,k) = grid%blk(b)%cond_phase(p)%residual(v,i,j,k) + flux(v)
      enddo

    enddo
    !$OMP END PARALLEL
    if (bad_recon) then
      write(*,'(A)') ' [ERROR] [ICE::compute_bound] unphysical state at first order on a boundary face'
      error stop 1
    endif

  end subroutine compute_bound


end module ICE_Mod_BC_Fluxes
