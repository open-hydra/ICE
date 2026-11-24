module ICE_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block

  implicit none
  private
  public :: compute_dt, set_dt_global, bound_empty_dt

contains

  subroutine compute_dt (p, cfl, cfl_rampa_iter, dt_max, tau_factor, grid)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_sim_param
    use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, thread_group, cell_slice
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    integer(kind=I4), intent(in)  :: p
    real(kind=R8), intent(in)     :: cfl
    integer(kind=I4), intent(in)  :: cfl_rampa_iter
    real(kind=R8), intent(in)     :: dt_max
    real(kind=R8), intent(in)     :: tau_factor
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, i, j, k
    integer(kind=I4) :: n, tid, g, m, s, i1, j1, k1, i2, j2, k2
    real(kind=R8)    :: dtmin, tau
    logical          :: tau_limit

    tau_limit = (obj_sim_param%owcoupled .or. obj_sim_param%twcoupled) .and. tau_factor > 0._R8

    !> Per-thread minimum, merged once below: updating grid%dtglobal inside the
    !> worksharing loop would race between threads, and so would reading it here
    !> while another thread merges (the blocks carry no barrier between them).
    dtmin = huge(dtmin)

    !> Thread groups: each thread its slice of its group's blocks, the same statements
    if (n_tgroups > 1) then
      tid = 0
      !$ tid = omp_get_thread_num()
      call thread_group(tid, g, m, s)
      do n = tg_first(g), tg_first(g+1) - 1
        b = tg_blocks(n)
        call cell_slice(grid%blk(b)%dim, m, s, i1, j1, k1, i2, j2, k2)
        do k = k1, k2
        do j = merge(j1, 1, k == k1), merge(j2, grid%blk(b)%dim(2), k == k2)
        do i = merge(i1, 1, k == k1 .and. j == j1), merge(i2, grid%blk(b)%dim(1), k == k2 .and. j == j2)
            grid%blk(b)%cond_phase(p)%dt(i,j,k) = 1e+5
            call compute_dt_ (grid%blk(b)%cond_phase(p)%prim(:,i,j,k), &
                              grid%blk(b)%M(i,j,k)%c,                  &
                              grid%blk(b)%cond_phase(p)%dt(i,j,k))
            grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%blk(b)%cond_phase(p)%dt(i,j,k) * cfl
            if (grid%iter <= cfl_rampa_iter) then
              grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%blk(b)%cond_phase(p)%dt(i,j,k) * grid%iter / cfl_rampa_iter
            endif
            grid%blk(b)%cond_phase(p)%dt(i,j,k) = min(grid%blk(b)%cond_phase(p)%dt(i,j,k), dt_max)
            if (tau_limit) then
              tau = grid%blk(b)%cond_phase(p)%tau(i,j,k)
              if (grid%blk(b)%cond_phase(p)%prim(1,i,j,k) > rho_empty .and. ieee_is_finite(tau) .and. tau > 0._R8) &
                grid%blk(b)%cond_phase(p)%dt(i,j,k) = min(grid%blk(b)%cond_phase(p)%dt(i,j,k), tau_factor*tau)
            endif
            dtmin = min (dtmin, grid%blk(b)%cond_phase(p)%dt(i,j,k))
        enddo ; enddo ; enddo
      enddo
    else
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle

      !$omp do collapse(3) schedule(runtime)
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)

        grid%blk(b)%cond_phase(p)%dt(i,j,k) = 1e+5

        call compute_dt_ (grid%blk(b)%cond_phase(p)%prim(:,i,j,k), &
                          grid%blk(b)%M(i,j,k)%c,                  &
                          grid%blk(b)%cond_phase(p)%dt(i,j,k))

         
        grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%blk(b)%cond_phase(p)%dt(i,j,k) * cfl

        if (grid%iter <= cfl_rampa_iter) then
          grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%blk(b)%cond_phase(p)%dt(i,j,k) * grid%iter / cfl_rampa_iter
        endif

        !> dt-max is a ceiling on the step actually taken, so it comes after the CFL
        !  factor and after the ramp. Its job is to bound the explicit source terms,
        !  whose relaxation times the CFL condition knows nothing about.
        grid%blk(b)%cond_phase(p)%dt(i,j,k) = min(grid%blk(b)%cond_phase(p)%dt(i,j,k), dt_max)

        !> Coupled runs: the explicit relaxation sources bound the step by the relaxation time
        !  of the cells that carry the phase; an empty cell has no meaningful tau
        if (tau_limit) then
          tau = grid%blk(b)%cond_phase(p)%tau(i,j,k)
          if (grid%blk(b)%cond_phase(p)%prim(1,i,j,k) > rho_empty .and. ieee_is_finite(tau) .and. tau > 0._R8) &
            grid%blk(b)%cond_phase(p)%dt(i,j,k) = min(grid%blk(b)%cond_phase(p)%dt(i,j,k), tau_factor*tau)
        endif

        dtmin = min (dtmin, grid%blk(b)%cond_phase(p)%dt(i,j,k))

      enddo ; enddo ; enddo
      !$omp end do nowait
    enddo
    end if

    !$omp critical (ice_dtglobal)
    grid%dtglobal = min (grid%dtglobal, dtmin)
    !$omp end critical (ice_dtglobal)
    !$omp barrier

  end subroutine compute_dt    
  
  



  !> The step of the cells a family does not occupy (bulk density at or below rho_empty) where a source fills them (the
  !> births of a generation, srate; another solver's generated mass, ext; in any cell) or where the family burns
  !> (empties: its drops burn down to ICE's empty state and its flux fills the cells beyond). Their own state bounds
  !> nothing (velocity ~0, no relaxation time), but the source fills them in the first stage of this step, at its own
  !> velocity, or the cells it fills fill them, and the later stages and the reconstruction move that content with the
  !> cell's step.
  !> Steady (each cell its own step): an empty cell steps no further than the family's smallest step over the cells it
  !> occupies and over the empty cells the source fills (at the source's velocity: its momentum over its mass), nor
  !> than its own width allows at the family's largest speed, each Cartesian component's largest magnitude over those
  !> cells (sum_d sum_c umax_c |M(d,c)| bounds sum_d |u.M(d,:)| for each of their velocities: on a graded mesh an empty
  !> cell can be narrower than every occupied one).
  !> Time-accurate (one step, compute_dt's minimum over the cells, the occupied ones with it): the step no further than
  !> the source's step in the empty cells it fills; the minimum over the ranks follows, with the step's own.
  !> A family no source fills and that does not burn keeps compute_dt's steps, bit for bit. A steady run's bound takes
  !> this rank's cells: a source on another rank fills this rank's cells through a connection face unbounded.
  !> Outside parallel regions.
  subroutine bound_empty_dt (p, cfl, grid, empties)
    use ICE_Global_m,         only: rho_empty, ndir
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Config_Types_m,   only: obj_time_scheme
    implicit none
    integer(kind=I4),      intent(in)    :: p
    real(kind=R8),         intent(in)    :: cfl
    type(ICE_domain_type), intent(inout) :: grid
    logical, optional,     intent(in)    :: empties
    integer(kind=I4) :: b, i, j, k, d
    real(kind=R8)    :: floor_dt, src_dt, dts, s(4), umax(3), rate
    real(kind=R8), allocatable :: prim_s(:)
    logical          :: filled, any_src

    filled = .false.
    if (present(empties)) filled = empties
    any_src = .false.
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      any_src = any_src .or. allocated(grid%blk(b)%cond_phase(p)%srate) .or. allocated(grid%blk(b)%cond_phase(p)%ext)
    enddo
    if (.not. (filled .or. any_src)) return
    floor_dt = huge(floor_dt)
    src_dt   = huge(src_dt)
    umax     = 0._R8
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      associate (cp => grid%blk(b)%cond_phase(p))
      if (.not. allocated(prim_s)) allocate(prim_s(size(cp%prim, 1)))
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)
        s = 0._R8
        if (allocated(cp%srate)) s = s + cp%srate(1:4,i,j,k)
        if (allocated(cp%ext))   s = s + cp%ext(1:4,i,j,k)
        if (s(1) > 0._R8) filled = .true.
        if (cp%prim(1,i,j,k) > rho_empty) then
          floor_dt = min(floor_dt, cp%dt(i,j,k))
          umax     = max(umax, abs(cp%prim(2:4,i,j,k)))
          cycle
        endif
        if (.not. s(1) > 0._R8) cycle
        prim_s = cp%prim(:,i,j,k)
        prim_s(2:4) = s(2:4)/s(1)
        umax = max(umax, abs(prim_s(2:4)))
        dts = 1e+5
        call compute_dt_(prim_s, grid%blk(b)%M(i,j,k)%c, dts)
        src_dt = min(src_dt, dts*cfl)
      enddo ; enddo ; enddo
      end associate
    enddo
    if (.not. filled) return
    if (obj_time_scheme%time_accurate) then
      grid%dtglobal = min(grid%dtglobal, src_dt)
      return
    endif
    floor_dt = min(floor_dt, src_dt)
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      associate (cp => grid%blk(b)%cond_phase(p))
      do k = 1, grid%blk(b)%dim(3)
      do j = 1, grid%blk(b)%dim(2)
      do i = 1, grid%blk(b)%dim(1)
        if (cp%prim(1,i,j,k) > rho_empty) cycle
        rate = 0._R8
        do d = 1, ndir
          rate = rate + sum(umax*abs(grid%blk(b)%M(i,j,k)%c(d,1:3)))
        enddo
        dts = floor_dt
        if (rate > 0._R8) dts = min(dts, cfl/rate)
        cp%dt(i,j,k) = min(cp%dt(i,j,k), dts)
      enddo ; enddo ; enddo
      end associate
    enddo

  end subroutine bound_empty_dt


  !> The step of a cell, one rule over the directions the mesh resolves:
  !>   1/dt = sum_d ( |u . grad xi_d| + a_d |grad xi_d| ),  d = 1, ..., ndir,
  !> with grad xi_d the row d of the metric tensor, normal to the d faces with the modulus of the
  !> inverse width of the cell across them (the edge of a skewed cell is longer), and a_d the signal
  !> speed of the closure along it. A scheme that is not split by direction needs the sum of the
  !> direction limits, not the smallest of them. The directions of a single cell add no limit: their
  !> faces give a mode only to the components odd across mirror faces, which setup_freeze does not
  !> integrate (it warns where the faces are not mirror planes). With the linearised MUSCL scheme
  !> (the van Leer, van Albada and MC slopes in smooth flow) the rule is stable for cfl <= 1 with RK2
  !> (Heun) and cfl <= 1.17 with RK3. The caller applies the cfl factor, the ramp and the ceilings; a
  !> cell with nothing to carry (rate 0) keeps the dt it was given.
  subroutine compute_dt_ (prim, tensor, dt)
    use ICE_Lib_Model
    use ICE_Global_m, only: ndir
    implicit none
    real(kind=R8), dimension(:),   intent(in)    :: prim
    real(kind=R8), dimension(:,:), intent(in)    :: tensor
    real(kind=R8),                 intent(inout) :: dt

    integer(kind=I4) :: d
    real(kind=R8)    :: versor(3), gn, rate

    rate = 0._R8
    do d = 1, ndir
      gn     = norm2 ( tensor(d,:) )
      versor = tensor(d,:) / gn
      rate   = rate + ( abs( dot_product (prim(2:4), versor) ) + wavespeed_make(prim, versor) ) * gn
    enddo
    if (rate > 0._R8) dt = min (dt, 1._R8/rate)

  end subroutine compute_dt_

  
  subroutine set_dt_global (grid)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, thread_group, cell_slice
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: b, p, i, j, k
    integer(kind=I4) :: n, tid, g, m, s, i1, j1, k1, i2, j2, k2

    !> Thread groups: each thread its slice of its group's blocks (orphaned: no worksharing)
    if (n_tgroups > 1) then
      tid = 0
      !$ tid = omp_get_thread_num()
      call thread_group(tid, g, m, s)
      do n = tg_first(g), tg_first(g+1) - 1
        b = tg_blocks(n)
        call cell_slice(grid%blk(b)%dim, m, s, i1, j1, k1, i2, j2, k2)
        do p = 1, ngroups
        do k = k1, k2
        do j = merge(j1, 1, k == k1), merge(j2, grid%blk(b)%dim(2), k == k2)
        do i = merge(i1, 1, k == k1 .and. j == j1), merge(i2, grid%blk(b)%dim(1), k == k2 .and. j == j2)
            grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%dtglobal
        enddo ; enddo ; enddo
        enddo
      enddo
      return
    end if

    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      do p = 1, ngroups

        !$omp do collapse(3) schedule(runtime)
        do k = 1, grid%blk(b)%dim(3)
        do j = 1, grid%blk(b)%dim(2)
        do i = 1, grid%blk(b)%dim(1)
        
          grid%blk(b)%cond_phase(p)%dt(i,j,k) = grid%dtglobal
        
        end do ; end do ; end do
        !$omp end do nowait
      enddo
    enddo  
  
  end subroutine set_dt_global
  

end module ICE_Mod_dt