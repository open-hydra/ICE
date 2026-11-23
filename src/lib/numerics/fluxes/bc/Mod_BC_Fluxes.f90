module ICE_Mod_BC_Fluxes
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Lib_Ghost, only: compute_ghost, fill_second_ghost
  use ICE_Mod_GhostExchange, only: exchange_conn_flux, exchange_conn_beta
  use ICE_Mod_ThreadGroups, only: n_tgroups, thread_group, split_range
  !$ use omp_lib, only: omp_get_thread_num

  implicit none
  private
  public :: compute_ghost, fill_second_ghost, compute_bound, divert_save, divert_burnout

  !> divert_save's residuals, per block (allocated on the first call)
  type :: res_keep_t
    real(kind=R8), allocatable :: a(:,:,:,:)
  end type res_keep_t
  type(res_keep_t), allocatable :: res_keep(:)

  !> The steady burn-out's band, in multiples of the removal size (d-burnout, else m_burnout's drop): all of a
  !> cell's outflow burns at its exit at or below burnout_band(1) times it, none at or above burnout_band(2) times
  !> it. The one copy of its value: the setup prints it (Wrap_Setup), the gate p31 reads it from there
  real(kind=R8), parameter, public :: burnout_band(2) = [0.90_R8, 1._R8]
  real(kind=R8), parameter :: band_lo = burnout_band(1), band_hi = burnout_band(2)

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
    !> With MPI, the partner's rank first sends its cell's detector weight to the owner's
    !> (exchange_conn_beta), and the owner's rank then sends the flux to the partner's
    !> (exchange_conn_flux).
    call exchange_conn_beta(grid, p)
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


  !> The flux of connection face grid%bc(n)%cf of family p, from its owner record n. A connection
  !> face lies inside the domain, so it takes the shock-detector weight of an interior face, the
  !> smaller of its two cells' (face_flux); the partner's cell is in another block, and its weight
  !> comes from conn(p)%bp when that block is on another rank (exchange_conn_beta).
  subroutine owner_flux(grid, p, n)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4),      intent(in)    :: p, n
    integer(kind=I4) :: q
    real(kind=R8)    :: flux(ncond_max), beta_q

    q = grid%bc(n)%partner
    if (is_local_block(grid%bc(q)%b)) then
      beta_q = grid%blk(grid%bc(q)%b)%cond_phase(p)%beta(grid%bc(q)%i, grid%bc(q)%j, grid%bc(q)%k)
    else
      beta_q = grid%conn(p)%bp(grid%bc(n)%cf)
    end if
    call record_flux(grid, p, n, flux(1:ncond(p)), &
                     min(grid%blk(grid%bc(n)%b)%cond_phase(p)%beta(grid%bc(n)%i, grid%bc(n)%j, grid%bc(n)%k), beta_q))
    grid%conn(p)%F(1:ncond(p), grid%bc(n)%cf) = flux(1:ncond(p))

  end subroutine owner_flux


  !> The flux through the face of boundary record n of family p into its cell, with the shock-
  !> detector weight beta_face, or the cell's own. A first ghost
  !> that holds an injected state (inj) has length 0: its side of the face is that state, and
  !> the cell's side takes the one-sided slope from the face; the length between the two
  !> ghosts stays, as both hold the state and their slope is 0. In a block one cell across,
  !> ip is the opposite face's first ghost, which inj may flag as well; inj is never set on
  !> an interior cell.
  subroutine record_flux(grid, p, n, flux, beta_face)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Lib_Riemann
    use ICE_Mod_Fluxes, only: state_reconstruction, bad_recon
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    type(ICE_domain_type), intent(in)  :: grid
    integer(kind=I4),      intent(in)  :: p, n
    real(kind=R8),         intent(out) :: flux(:)
    real(kind=R8), optional, intent(in) :: beta_face
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
      !> A physical boundary face or a chimera face (102) takes its cell's shock-detector weight, the
      !> only cell it has: a chimera ghost is a blend of donor cells, not a cell (a connection face
      !> is given the interior faces' weight by owner_flux)
      if (present(beta_face)) then
        beta_val = beta_face
      else
        beta_val = grid%blk(b)%cond_phase(p)%beta(i,j,k)
      end if

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

  !> Steady runs, a burning family without a residue (d-stop 0): the steady burn-out diverts this family (its stages
  !> keep the residual after compute_bound for divert_burnout)
  logical function diverts(p)
    use ICE_Global_m,       only: mat_of
    use ICE_Config_Types_m, only: obj_condensed, obj_time_scheme
    implicit none
    integer(kind=I4), intent(in) :: p
    diverts = .not. obj_time_scheme%time_accurate .and. obj_condensed(mat_of(p))%combSelect > 0
    if (diverts) diverts = .not. obj_condensed(mat_of(p))%dRes > 0._R8
  end function diverts

  !> The residual of a family the steady burn-out diverts, after compute_bound and before compute_flux: the boundary
  !> fluxes, from which divert_burnout rebuilds the residual of the cells it changes. Outside parallel regions.
  subroutine divert_save (grid, p)
    use ICE_Global_m,         only: ncond
    use ICE_Advanced_Types_m, only: ICE_domain_type
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer(kind=I4),      intent(in) :: p
    integer(kind=I4) :: b

    if (.not. diverts(p)) return
    if (.not. allocated(res_keep)) allocate(res_keep(grid%nb))
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      res_keep(b)%a = grid%blk(b)%cond_phase(p)%residual(1:ncond(p), 1:grid%blk(b)%dim(1), 1:grid%blk(b)%dim(2), &
                                                          1:grid%blk(b)%dim(3))
    enddo
  end subroutine divert_save

  !> Steady runs, a burning family without a residue (d-stop 0): the steady form of burnout_switch, smoothed. A cell
  !> that holds drops (above ICE's empty state, as the switch tests) which burn (direct at their burning temperature,
  !> Beckstead from T-ign) and lie within the band of the removal size (d-burnout, else a drop mass of m_burnout) is
  !> flagged with a fraction f (burnout_fraction: 1 at or below band_lo times it, 0 at or above band_hi times it):
  !> f of what flows out of it through an interior face or a physical boundary face does not reach the cell beyond
  !> (or leave the domain), it burns at the cell's exit, a sink of the cell in source and in burn, so the coupler
  !> hands it to the gas as every burn, and the mass it generates is born at the cell's velocity and temperature; 1 - f
  !> of the face's flux reaches the cell beyond. The diverted outflow is the flux of the cell's own state on both sides
  !> of the face (its own velocity: the MK flux upwinds at the mean of the two cells' velocities, and the cell beyond is
  !> emptied by this), so at f = 1 the cell's content, and so its drop size, does not depend on the cell beyond, and a
  !> cell fed only by cells at f = 1 holds nothing at the fixed point, as the switch leaves it. A flag that is only 0
  !> or 1 cycled with MUSCL (a cell at the threshold flipped, and the faces around it followed); f moves the diverted
  !> flux smoothly with the drops' size. A connection (101/201) or chimera (102) face carries its flux as it is: the
  !> cell across it tests itself. After compute_flux, before compute_residual: the residual of each flagged cell and of
  !> each cell next to one is rebuilt from the one divert_save kept, the interior faces in compute_flux's order
  !> (tile_sweep), a flagged cell's faces at first order in the measure of the larger f of its two cells, a face whose
  !> flux leaves a flagged cell carried in the measure 1 - f by both its cells (a cell none of whose faces touches a
  !> flagged cell gets compute_flux's residual bit for bit; a flux is never added and taken off again, which would
  !> leave its round-off in an empty cell; at f = 1 every face is a full diversion's, bit for bit).
  subroutine divert_burnout (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Parameters_m,   only: pi
    use ICE_Config_Types_m, only: obj_condensed, obj_time_scheme, COMB_DIRECT
    use ICE_Lib_Model,      only: cons_2_prim
    use ICE_Lib_Properties, only: mat_rho
    use ICE_Mod_Fluxes,     only: face_flux
    use ICE_Mod_Freeze,     only: nfrozen, frozen
    use ICE_Lib_Riemann,    only: riemann
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4),      intent(in)    :: p
    integer(kind=I4) :: b, i, j, k, n, nc, ni, nj, nk, c, r, m, nf
    real(kind=R8)    :: F(ncond(p)), Ton, mb, zero(ncond(p)), empty(ncond(p)), empty1
    real(kind=R8), allocatable :: sink(:,:,:,:), frac(:,:,:)
    logical,       allocatable :: flag(:,:,:), redo(:,:,:)

    if (.not. diverts(p)) return
    n  = nbase(p)
    nc = ncond(p)
    zero   = 0._R8
    empty  = cons_2_prim(zero, obj_condensed(mat_of(p)))
    empty1 = empty(1)
    nf = 0
    if (allocated(nfrozen)) nf = nfrozen(p)
    Ton = obj_condensed(mat_of(p))%Tign
    if (obj_condensed(mat_of(p))%combSelect == COMB_DIRECT) Ton = obj_condensed(mat_of(p))%Tboil*(1._R8 - 1.e-10_R8)

    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
      ! padded by a cell on each side, so that a face's two cells need no bounds test
      allocate(flag(0:ni+1,0:nj+1,0:nk+1), frac(0:ni+1,0:nj+1,0:nk+1))
      flag = .false.
      frac = 0._R8
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        associate (q => grid%blk(b)%cond_phase(p)%prim(:,i,j,k))
          flag(i,j,k) = q(1) > empty1 .and. q(n) > 0._R8
          if (flag(i,j,k)) flag(i,j,k) = q(n-1) >= Ton
          if (flag(i,j,k)) then
            mb = m_burnout
            if (obj_condensed(mat_of(p))%dBurn > 0._R8) &
              mb = pi/6._R8*mat_rho(obj_condensed(mat_of(p)), q(n-1))*obj_condensed(mat_of(p))%dBurn**3
            frac(i,j,k) = burnout_fraction(q(1), mb*q(n))
            flag(i,j,k) = frac(i,j,k) > 0._R8
          end if
        end associate
      enddo ; enddo ; enddo
      if (.not. any(flag)) then
        deallocate(flag, frac)
        cycle
      end if
      allocate(sink(nc,ni,nj,nk), redo(ni,nj,nk))
      sink = 0._R8
      redo = .false.
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        if (.not. flag(i,j,k)) cycle
        redo(max(i-1,1):min(i+1,ni), j, k) = .true.
        redo(i, max(j-1,1):min(j+1,nj), k) = .true.
        redo(i, j, max(k-1,1):min(k+1,nk)) = .true.
        ! the cell's own outflow through each interior face, the face above (s = 1) and the face below (s = -1)
        if (i < ni) call own_out(1, i,   j, k,  1, i, j, k)
        if (i > 1)  call own_out(1, i-1, j, k, -1, i, j, k)
        if (j < nj) call own_out(2, i, j,   k,  1, i, j, k)
        if (j > 1)  call own_out(2, i, j-1, k, -1, i, j, k)
        if (nk > 1) then
          if (k < nk) call own_out(3, i, j, k,    1, i, j, k)
          if (k > 1)  call own_out(3, i, j, k-1, -1, i, j, k)
        end if
      enddo ; enddo ; enddo
      ! the residuals, rebuilt as tile_sweep builds them: each direction in the odd-then-even face order
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        if (.not. redo(i,j,k)) cycle
        associate (rs => grid%blk(b)%cond_phase(p)%residual(1:nc,i,j,k))
          rs = res_keep(b)%a(:,i,j,k)
          if (mod(i,2) == 1) then
            if (i < ni) rs = rs - part(1, i,   j, k)
            if (i > 1)  rs = rs + part(1, i-1, j, k)
          else
            rs = rs + part(1, i-1, j, k)
            if (i < ni) rs = rs - part(1, i, j, k)
          end if
          if (mod(j,2) == 1) then
            if (j < nj) rs = rs - part(2, i, j,   k)
            if (j > 1)  rs = rs + part(2, i, j-1, k)
          else
            rs = rs + part(2, i, j-1, k)
            if (j < nj) rs = rs - part(2, i, j, k)
          end if
          if (nk > 1) then
            if (mod(k,2) == 1) then
              if (k < nk) rs = rs - part(3, i, j, k)
              if (k > 1)  rs = rs + part(3, i, j, k-1)
            else
              rs = rs + part(3, i, j, k-1)
              if (k < nk) rs = rs - part(3, i, j, k)
            end if
          end if
        end associate
      enddo ; enddo ; enddo
      ! the physical boundary faces of the flagged cells of this block: what leaves through them
      do c = 1, grid%bcells(p)%n
        r = grid%bcells(p)%rec(grid%bcells(p)%first(c))
        if (grid%bc(r)%b /= b) cycle
        if (.not. flag(grid%bc(r)%i, grid%bc(r)%j, grid%bc(r)%k)) cycle
        do r = grid%bcells(p)%first(c), grid%bcells(p)%first(c+1) - 1
          m = grid%bcells(p)%rec(r)
          select case (grid%bc(m)%type)
            case (101, 201, 102)
              cycle
            case (200)
              if (trim(obj_time_scheme%model(p)) == 'MK') cycle
          end select
          call record_flux(grid, p, m, F)
          if (F(1) < 0._R8) then
            i = grid%bc(m)%i ; j = grid%bc(m)%j ; k = grid%bc(m)%k
            F = frac(i,j,k)*F
            grid%blk(b)%cond_phase(p)%residual(1:nc,i,j,k) = grid%blk(b)%cond_phase(p)%residual(1:nc,i,j,k) - F
            sink(:,i,j,k) = sink(:,i,j,k) + F
          end if
        enddo
      enddo
      ! the sinks, per unit volume and time, in source and (the mass) in burn; the frozen slots exchange nothing
      if (nf > 0) sink(frozen(1:nf,p),:,:,:) = 0._R8
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        if (.not. flag(i,j,k)) cycle
        grid%blk(b)%cond_phase(p)%source(1:nc,i,j,k) = grid%blk(b)%cond_phase(p)%source(1:nc,i,j,k) &
                                                       + sink(:,i,j,k)/grid%blk(b)%vol(i,j,k)
        grid%blk(b)%cond_phase(p)%burn(i,j,k) = grid%blk(b)%cond_phase(p)%burn(i,j,k) &
                                                - sink(1,i,j,k)/grid%blk(b)%vol(i,j,k)
      enddo ; enddo ; enddo
      deallocate(flag, frac, sink, redo)
    enddo

  contains

    !> The flux of the face of direction dd between cell (il,jl,kl) and its upper neighbour, as compute_flux evaluates
    !> it, or at first order in the measure w, the larger f of its two cells (the reconstruction's stencil would reach
    !> the cell beyond, which the diversion empties); 1 - f of it when it leaves a diverted cell (the lower cell when
    !> it is positive, the upper one when negative)
    function part(dd, il, jl, kl) result(Ff)
      integer(kind=I4), intent(in) :: dd, il, jl, kl
      real(kind=R8) :: Ff(nc), Fm(nc), w
      integer(kind=I4) :: iu, ju, ku
      iu = il ; ju = jl ; ku = kl
      select case (dd)
        case (1) ; iu = il + 1
        case (2) ; ju = jl + 1
        case default ; ku = kl + 1
      end select
      w = max(frac(il,jl,kl), frac(iu,ju,ku))
      if (w > 0._R8) then
        ! a face of a flagged cell: first order, its two cells' states (a reconstruction would reach the cell beyond),
        ! blended with compute_flux's face below w = 1
        Ff = riemann(grid%blk(b)%cond_phase(p)%prim(1:nc,il,jl,kl), grid%blk(b)%cond_phase(p)%prim(1:nc,iu,ju,ku), &
                     grid%blk(b)%dir(dd)%f(il,jl,kl)%N, obj_condensed(mat_of(p)))*grid%blk(b)%dir(dd)%f(il,jl,kl)%A
        if (w < 1._R8) then
          call face_flux(grid%blk(b), b, p, nc, dd, il, jl, kl, Fm)
          Ff = w*Ff + (1._R8 - w)*Fm
        end if
      else
        call face_flux(grid%blk(b), b, p, nc, dd, il, jl, kl, Ff)
      end if
      if (Ff(1) > 0._R8 .and. frac(il,jl,kl) > 0._R8) then
        Ff = keep(frac(il,jl,kl), Ff)
      else if (Ff(1) < 0._R8 .and. frac(iu,ju,ku) > 0._R8) then
        Ff = keep(frac(iu,ju,ku), Ff)
      end if
    end function part

    !> 1 - f of a face's flux, exactly zero at f = 1 (a full diversion's zero, not a signed one)
    pure function keep(f, Fin) result(Fk)
      real(kind=R8), intent(in) :: f, Fin(:)
      real(kind=R8) :: Fk(size(Fin))
      Fk = 0._R8
      if (f < 1._R8) Fk = (1._R8 - f)*Fin
    end function keep

    !> The flagged cell (ic,jc,kc), on the lower side (s = 1) or the upper side (s = -1) of the face of direction dd
    !> between cell (il,jl,kl) and its upper neighbour: f of the flux of its own state on both sides, to its sink when
    !> it points out of the cell
    subroutine own_out(dd, il, jl, kl, s, ic, jc, kc)
      integer(kind=I4), intent(in) :: dd, il, jl, kl, s, ic, jc, kc
      real(kind=R8) :: Fo(nc)
      Fo = riemann(grid%blk(b)%cond_phase(p)%prim(1:nc,ic,jc,kc), grid%blk(b)%cond_phase(p)%prim(1:nc,ic,jc,kc), &
                   grid%blk(b)%dir(dd)%f(il,jl,kl)%N, obj_condensed(mat_of(p)))*grid%blk(b)%dir(dd)%f(il,jl,kl)%A
      if (s*Fo(1) > 0._R8) sink(:,ic,jc,kc) = sink(:,ic,jc,kc) - frac(ic,jc,kc)*s*Fo
    end subroutine own_out

  end subroutine divert_burnout


  !> The fraction of a diverted cell's outflow that burns at its exit, from the drops' bulk density m against mb (the
  !> drops' number times d-burnout's drop mass, else m_burnout): 1 at or below band_lo**3*mb, 0 at or above
  !> band_hi**3*mb, a cubic step in the diameter ratio between (C1)
  pure function burnout_fraction(m, mb) result(f)
    real(kind=R8), intent(in) :: m, mb
    real(kind=R8) :: f, t
    f = 0._R8
    if (m <= band_lo**3*mb) then
      f = 1._R8
    else if (m < band_hi**3*mb) then
      t = (band_hi - (m/mb)**(1._R8/3._R8))/(band_hi - band_lo)
      f = t*t*(3._R8 - 2._R8*t)
    end if
  end function burnout_fraction


end module ICE_Mod_BC_Fluxes
