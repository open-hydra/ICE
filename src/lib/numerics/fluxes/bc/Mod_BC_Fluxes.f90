module ICE_Mod_BC_Fluxes
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Lib_Ghost, only: compute_ghost, fill_second_ghost
  use ICE_Mod_GhostExchange, only: exchange_conn_flux
  use ICE_Mod_ThreadGroups, only: n_tgroups, thread_group, split_range
  !$ use omp_lib, only: omp_get_thread_num

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
    integer(kind=I4) :: c, tid, g, m, s, a, z, jj

    !> The residual was zeroed by zero_residual at the start of the stage and
    !> nothing has written it since, so the boundary fluxes are its first terms.
    bad_recon = .false.

    !> The flux of every connection face (101/201) of family p, once, by the face's owner record,
    !> before any cell takes it: the owner's cell then adds it and the partner's cell its negative
    !> (bound_cell), so the two blocks exchange the same flux. The end of the region is the barrier
    !> before the cells: a partner's cell may belong to another thread group than its owner's.
    !> With MPI, the owner's rank sends the flux to the partner's (exchange_conn_flux).
    if (grid%conn(p)%nown > 0) then
      if (n_tgroups > 1) then
        !$OMP PARALLEL DEFAULT(NONE) SHARED(grid, p) PRIVATE(tid, g, m, s, a, z, jj)
        tid = 0
        !$ tid = omp_get_thread_num()
        call thread_group(tid, g, m, s)
        call split_range(grid%conn(p)%tg_first(g), grid%conn(p)%tg_first(g+1) - 1, m, s, a, z)
        do jj = a, z
          call owner_flux(grid, p, grid%conn(p)%own(jj))
        end do
        !$OMP END PARALLEL
      else
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(grid, p) PRIVATE(jj) SCHEDULE(DYNAMIC, 16)
        do jj = 1, grid%conn(p)%nown
          call owner_flux(grid, p, grid%conn(p)%own(jj))
        end do
        !$OMP END PARALLEL DO
      end if
    end if
    call exchange_conn_flux(grid, p)

    !> Thread groups: each thread its static share of its own group's boundary cells
    if (n_tgroups > 1) then
      !$OMP PARALLEL DEFAULT(NONE) SHARED(grid, p) PRIVATE(tid, g, m, s, a, z, jj)
      tid = 0
      !$ tid = omp_get_thread_num()
      call thread_group(tid, g, m, s)
      call split_range(grid%bcells(p)%tg_first(g), grid%bcells(p)%tg_first(g+1) - 1, m, s, a, z)
      do jj = a, z
        call bound_cell(grid, p, grid%bcells(p)%tg_cell(jj))
      end do
      !$OMP END PARALLEL
    else
    !$OMP PARALLEL DEFAULT(NONE) SHARED(grid, p) PRIVATE(c)

    !> One iteration per boundary cell of this rank (see build_local_bc_index),
    !> its records in face order: an edge or corner cell receives its
    !> two or three fluxes from one thread in one fixed order, whatever the thread
    !> count, and no accumulation is shared.
    !$OMP DO SCHEDULE (DYNAMIC, 16)
    do c = 1, grid%bcells(p)%n
      call bound_cell(grid, p, c)
    enddo
    !$OMP END PARALLEL
    end if
    if (bad_recon) then
      write(*,'(A)') ' [ERROR] [ICE::compute_bound] unphysical state at first order on a boundary face'
      error stop 1
    endif

  end subroutine compute_bound


  !> The boundary fluxes of boundary cell c of family p, its records in face order: the
  !> body of compute_bound's loop.
  subroutine bound_cell(grid, p, c)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_time_scheme
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4),      intent(in)    :: p
    integer(kind=I4),      intent(in)    :: c
    integer(kind=I4) :: n, r, b, i, j, k
    real(kind=R8)    :: flux(ncond_max)

    do r = grid%bcells(p)%first(c), grid%bcells(p)%first(c+1) - 1
      n = grid%bcells(p)%rec(r)

      b = grid%bc(n)%b
      i = grid%bc(n)%i ; j = grid%bc(n)%j ; k = grid%bc(n)%k

      !> Nothing crosses a symmetry face of a pressureless cloud
      if (grid%bc(n)%type == 200 .and. trim(obj_time_scheme%model(p)) == 'MK') cycle

      select case (grid%bc(n)%type)
        case (101, 201)
          !> A connection face's flux is its owner record's, computed before the cells (compute_bound).
          !>  Every conserved component is in the global frame, so what enters the owner's cell leaves the
          !>  partner's whatever the orientation of the two blocks' indices
          if (grid%bc(n)%owner) then
            flux(1:ncond(p)) =  grid%conn(p)%F(1:ncond(p), grid%bc(n)%cf)
          else
            flux(1:ncond(p)) = -grid%conn(p)%F(1:ncond(p), grid%bc(n)%cf)
          end if
        case default
          call record_flux(grid, p, n, flux(1:ncond(p)))
      end select

      grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) = &
        grid%blk(b)%cond_phase(p)%residual(1:ncond(p),i,j,k) + flux(1:ncond(p))

    enddo

  end subroutine bound_cell


  !> The flux of connection face grid%bc(n)%cf of family p, from its owner record n (record_flux,
  !> which gives the face the shock-detector weight of the owner's cell: under MPI the partner's
  !> cell may lie on another rank, whose weights this rank does not compute).
  subroutine owner_flux(grid, p, n)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4),      intent(in)    :: p, n
    real(kind=R8)    :: flux(ncond_max)

    call record_flux(grid, p, n, flux(1:ncond(p)))
    grid%conn(p)%F(1:ncond(p), grid%bc(n)%cf) = flux(1:ncond(p))

  end subroutine owner_flux


  !> The flux through the face of boundary record n of family p into its cell. A first ghost
  !> that holds an injected state (inj) has length 0: its side of the face is that state, and
  !> the cell's side takes the one-sided slope from the face; the length between the two
  !> ghosts stays, as both hold the state and their slope is 0. In a block one cell across,
  !> ip is the opposite face's first ghost, which inj may flag as well; inj is never set on
  !> an interior cell.
  subroutine record_flux(grid, p, n, flux)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Lib_Riemann
    use ICE_Mod_Fluxes, only: state_reconstruction, bad_recon
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    type(ICE_domain_type), intent(in)  :: grid
    integer(kind=I4),      intent(in)  :: p, n
    real(kind=R8),         intent(out) :: flux(:)
    integer(kind=I4) :: b, f, i, j, k
    integer(kind=I4) :: ig, jg, kg, ig2, jg2, kg2, ip, jp, kp
    integer(kind=I4) :: dir
    real(kind=R8)    :: normal(3), area
    real(kind=R8)    :: dl0, dl1, dl2, dll, dlr, dl_g1, dl_m, dl_4th
    real(kind=R8)    :: beta_val
    real(kind=R8)    :: priml(ncond_max), primr(ncond_max)

      b = grid%bc(n)%b
      i = grid%bc(n)%i ; j = grid%bc(n)%j ; k = grid%bc(n)%k
      f = grid%bc(n)%f

      !> Ghost and interior neighbor coordinates
      ig  = i -   guide(f,1) ; jg  = j -   guide(f,2) ; kg  = k -   guide(f,3)
      ig2 = i - 2*guide(f,1) ; jg2 = j - 2*guide(f,2) ; kg2 = k - 2*guide(f,3)
      ip  = i +   guide(f,1) ; jp  = j +   guide(f,2) ; kp  = k +   guide(f,3)

      dir      = (f+1)/2
      !> The face's shock-detector weight is its cell's, the only one of its two cells the block holds
      !> (an interior face takes the smaller of its two cells' weights, face_flux)
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
        if (grid%blk(b)%cond_phase(p)%inj(ig,jg,kg)) then
          dl1 = 0.5_R8*dl_m ;  dll = 0._R8
        end if
        if (grid%blk(b)%cond_phase(p)%inj(ip,jp,kp)) dl2 = 0.5_R8*dl_m

        call state_reconstruction(grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig2,jg2,kg2), &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig,jg,kg),    &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),i,j,k),       &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ip,jp,kp),    &
                                  dl0, dl1, dl2, dll, dlr,                                 &
                                  priml(1:ncond(p)), primr(1:ncond(p)), beta_val,        &
                                  obj_condensed(mat_of(p)), solid_of(p), [b, i, j, k, f])

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
        if (grid%blk(b)%cond_phase(p)%inj(ig,jg,kg)) then
          dl1 = 0.5_R8*dl_m ;  dlr = 0._R8
        end if
        if (grid%blk(b)%cond_phase(p)%inj(ip,jp,kp)) dl0 = 0.5_R8*dl_m

        call state_reconstruction(grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ip,jp,kp),    &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),i,j,k),       &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig,jg,kg),    &
                                  grid%blk(b)%cond_phase(p)%prim(1:ncond(p),ig2,jg2,kg2), &
                                  dl0, dl1, dl2, dll, dlr,                                 &
                                  priml(1:ncond(p)), primr(1:ncond(p)), beta_val,        &
                                  obj_condensed(mat_of(p)), solid_of(p), [b, i, j, k, f])

        flux(1:ncond(p)) = - riemann(priml(1:ncond(p)), primr(1:ncond(p)), normal, obj_condensed(mat_of(p))) * area

      end select

  end subroutine record_flux


end module ICE_Mod_BC_Fluxes
