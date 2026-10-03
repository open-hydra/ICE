module ICE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Lib_Reconstruction, only: state_reconstruction, bad_recon
  use ICE_Mod_ThreadGroups, only: n_tgroups, tg_first, tg_blocks, tg_thr_lo, tg_thr_hi, tg_rounds, thread_group

  implicit none
  private
  public :: compute_flux, state_reconstruction, bad_recon

  !> The tile of one thread: the cells (i = 1..ni, j = j1..j2, k = k1..k2) of a block.
  !> Deep tiling cuts ranges of whole k-planes and keeps the k-faces below and above
  !> the plane being swept in two private plane buffers per thread; shallow tiling
  !> cuts pieces of a plane and computes every k-face of the block once into the
  !> shared kface array. Deep is taken only while the threads' plane buffers together
  !> fit plane_buffer_budget (on a 192 x 192 plane with seven variables a buffer is
  !> 2.1 MB: eighty threads of two each are 340 MB). The shallow tiling is not faster
  !> by itself -- at 80 threads on the 192^3 block the sweep took 218 ms against 221 --
  !> but its plane pieces give the dynamic schedule tiles to balance: with four per
  !> thread it took 181 ms (monolith, job 303578, steady state). Every thread builds
  !> the same list from the block's dimensions; the list holds tiles_per_thread tiles
  !> per thread.
  type :: tile_type
    integer :: k1 = 0, k2 = 0, j1 = 0, j2 = 0
  end type tile_type

  !> A k-face between two tiles is read by both: by the tile owning the cell
  !> below it and by the one owning the cell above. It is computed once, by the
  !> owner of the cell below, into a shared array -- one plane per tile when the
  !> tiles are ranges of planes (`seam`), every k-face of the block when they are
  !> pieces of a plane (`kface`). Grown to the largest block that needed them;
  !> freed with the program.
  real(R8), allocatable, private :: kface(:,:,:,:), seam(:,:,:,:)
  !> The same per thread group (last index), when the blocks are owned by groups.
  real(R8), allocatable, private :: kface_g(:,:,:,:,:), seam_g(:,:,:,:,:)

  !> Bytes the threads' private plane buffers may take together before the shallow
  !> tiling is used instead (the L3 of one node's socket, roughly).
  integer(kind=8), parameter :: plane_buffer_budget = 16_8 * 1024_8 * 1024_8

contains

  !> The interior fluxes of family p, gathered by cell. Every face flux is computed
  !> once by the thread that owns the cell below it (the lower i, j or k) and kept in
  !> that thread's line, row or plane buffer; each cell then receives the two faces
  !> of each direction in the order the odd-then-even face passes of the scatter form
  !> gave them -- an odd cell its own face first, then the one below; an even cell the
  !> one below first -- so the residual is bit for bit what those passes produced,
  !> without their six barriers and the strided copies of the j and k stencils.
  subroutine compute_flux (grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m,     only: obj_space_scheme, obj_condensed
    use ICE_Lib_Shock_Detector, only: SD_rho
    !$ use omp_lib,             only: omp_get_max_threads
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, t, nt, nthreads, ni, nj, nk, nc, ni_max, nj_max
    logical          :: deep
    !> tiles is each thread's inside the region; the serial pass before it uses its
    !> own tiles0: ifx 2023.1 reports "already allocated" at the first allocation of a
    !> PRIVATE copy whose original was allocated and freed before the region.
    type(tile_type), allocatable :: tiles(:), tiles0(:)
    !> Per-thread face buffers, sized once per region to the largest local block:
    !> the i-faces of a line, the j-faces below (jlo) and above (jhi) a row, the
    !> k-faces below (klo) and above (khi) a plane (deep blocks only). Moving up a
    !> row or a plane swaps the index: the faces above become the faces below.
    real(R8), allocatable :: fi(:,:), gj(:,:,:), hk(:,:,:,:)
    integer(kind=I4) :: jlo, jhi, klo, khi, sw

    if (n_tgroups > 1) then
      call compute_flux_grouped(grid, p)
      return
    end if

    nc = ncond(p)
    bad_recon = .false.
    !> The tiling is decided by the thread count the region will have at most, so
    !> that the shared seam or k-face planes every local block needs can be grown
    !> here, once, instead of by a SINGLE (and its barrier) per block.
    nthreads = 1
    !$ nthreads = omp_get_max_threads()
    ni_max = 0 ; nj_max = 0
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
      ni_max = max(ni_max, ni) ; nj_max = max(nj_max, nj)
      call make_tiles(ni, nj, nk, nc, nthreads, tiles0, deep)
      if (deep) then
        call grow(seam, nc, ni, nj, size(tiles0))
      else if (nk > 1) then
        call grow(kface, nc, ni, nj, nk)
      end if
    end do
    if (allocated(tiles0)) deallocate(tiles0)
    call grow(seam, 1, 1, 1, 1) ; call grow(kface, 1, 1, 1, 1)   ! both passed to the tile procedures
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(b, i, j, k, t, nt, ni, nj, nk, deep, tiles, &
    !$OMP                                  fi, gj, hk, jlo, jhi, klo, khi, sw)
    allocate(fi(nc, 0:ni_max), gj(nc, ni_max, 0:1), hk(nc, ni_max, nj_max, 0:1))
    jlo = 0 ; jhi = 1 ; klo = 0 ; khi = 1
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
      call make_tiles(ni, nj, nk, nc, nthreads, tiles, deep)
      nt = size(tiles)

      !> The shock-detector weight of every cell, before any face reads it: a face
      !> takes the weight of the cell below it, which may belong to another tile.
      !$OMP DO SCHEDULE (DYNAMIC, 1)
      do t = 1, nt
        call tile_beta(grid%blk(b), p, ni, tiles(t))
      enddo
      !$OMP END DO

      !> The k-faces shared between tiles, once, by the owner of the cell below
      !> each: the top faces of every tile's last plane (deep block), or every
      !> k-face of the block (shallow block).
      if (deep .or. nk > 1) then
        !$OMP DO SCHEDULE (DYNAMIC, 1)
        do t = 1, nt
          call tile_kfaces(grid%blk(b), p, nc, ni, nj, nk, tiles(t), t, deep, seam, kface)
        enddo
        !$OMP END DO
      end if

      !$OMP DO SCHEDULE (DYNAMIC, 1)
      do t = 1, nt
        call tile_sweep(grid%blk(b), p, nc, ni, nj, nk, tiles, t, deep, seam, kface, &
                        fi, gj, hk, jlo, jhi, klo, khi)
      enddo
      !$OMP END DO

      deallocate(tiles)
    enddo
    deallocate(fi, gj, hk)
    !$OMP END PARALLEL
    if (bad_recon) then
      write(*,'(A)') ' [ERROR] [ICE::compute_flux] unphysical state at first order on an interior face'
      error stop 1
    endif

  end subroutine compute_flux


  !> compute_flux with thread groups: the blocks in rounds -- round r sweeps the r-th
  !> block of every group at once, each group with its own threads, its own tiles
  !> (cut for the group's size) and its own shared planes -- so the barriers between
  !> the three phases stay barriers of the whole team, one set per round instead of
  !> one per block. Inside a group the tiles are taken dynamically, one at a time.
  subroutine compute_flux_grouped(grid, p)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, n, t, nt, ni, nj, nk, nc, ni_max, nj_max, g, m, s, r, tid, ph
    integer(kind=I4) :: s1, n1, n2, n3, n4
    logical          :: deep
    type(tile_type), allocatable :: tiles(:), tiles0(:)
    real(R8), allocatable :: fi(:,:), gj(:,:,:), hk(:,:,:,:)
    integer(kind=I4) :: jlo, jhi, klo, khi
    integer(kind=I4), allocatable :: cnt(:,:)

    nc = ncond(p)
    bad_recon = .false.
    ! Sizes of the per-group shared planes and of the thread buffers, from every group's blocks
    ni_max = 0 ; nj_max = 0 ; n1 = nc ; n2 = 1 ; n3 = 1 ; n4 = 1
    do g = 0, n_tgroups - 1
      s1 = tg_thr_hi(g) - tg_thr_lo(g) + 1
      do n = tg_first(g), tg_first(g+1) - 1
        b = tg_blocks(n)
        ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
        ni_max = max(ni_max, ni) ; nj_max = max(nj_max, nj)
        call make_tiles(ni, nj, nk, nc, s1, tiles0, deep)
        n2 = max(n2, ni) ; n3 = max(n3, nj) ; n4 = max(n4, nk, size(tiles0))
      end do
    end do
    if (allocated(tiles0)) deallocate(tiles0)
    call grow5(seam_g, n1, n2, n3, n4, n_tgroups)
    call grow5(kface_g, n1, n2, n3, n4, n_tgroups)
    allocate(cnt(0:n_tgroups-1, 3))

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(b, t, nt, ni, nj, nk, deep, tiles, g, m, s, r, tid, ph, &
    !$OMP                                  fi, gj, hk, jlo, jhi, klo, khi)
    allocate(fi(nc, 0:ni_max), gj(nc, ni_max, 0:1), hk(nc, ni_max, nj_max, 0:1))
    jlo = 0 ; jhi = 1 ; klo = 0 ; khi = 1
    tid = 0
    !$ tid = omp_get_thread_num()
    call thread_group(tid, g, m, s)
    do r = 1, tg_rounds
      b = 0 ; nt = 0 ; deep = .false. ; nk = 1
      if (r <= tg_first(g+1) - tg_first(g)) then
        b = tg_blocks(tg_first(g) + r - 1)
        ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
        call make_tiles(ni, nj, nk, nc, s, tiles, deep)
        nt = size(tiles)
      end if
      if (m == 0) cnt(g,:) = 0
      !$OMP BARRIER
      do ph = 1, 3
        if (ph == 2 .and. .not. (deep .or. nk > 1)) then
          !$OMP BARRIER
          cycle
        end if
        do
          !$OMP ATOMIC CAPTURE
          t = cnt(g, ph)
          cnt(g, ph) = cnt(g, ph) + 1
          !$OMP END ATOMIC
          t = t + 1
          if (t > nt) exit
          select case (ph)
          case (1)
            call tile_beta(grid%blk(b), p, ni, tiles(t))
          case (2)
            call tile_kfaces(grid%blk(b), p, nc, ni, nj, nk, tiles(t), t, deep, seam_g(:,:,:,:,g), kface_g(:,:,:,:,g))
          case (3)
            call tile_sweep(grid%blk(b), p, nc, ni, nj, nk, tiles, t, deep, seam_g(:,:,:,:,g), kface_g(:,:,:,:,g), &
                            fi, gj, hk, jlo, jhi, klo, khi)
          end select
        end do
        !$OMP BARRIER
      end do
      if (allocated(tiles)) deallocate(tiles)
    end do
    deallocate(fi, gj, hk)
    !$OMP END PARALLEL
    deallocate(cnt)
    if (bad_recon) then
      write(*,'(A)') ' [ERROR] [ICE::compute_flux] unphysical state at first order on an interior face'
      error stop 1
    endif

  end subroutine compute_flux_grouped


  !> grow for the per-group plane arrays.
  subroutine grow5(a, n1, n2, n3, n4, n5)
    real(R8), allocatable, intent(inout) :: a(:,:,:,:,:)
    integer, intent(in) :: n1, n2, n3, n4, n5
    if (allocated(a)) then
      if (size(a,1) < n1 .or. size(a,2) < n2 .or. size(a,3) < n3 .or. size(a,4) < n4 .or. size(a,5) < n5) deallocate(a)
    end if
    if (.not. allocated(a)) allocate(a(n1, n2, n3, n4, 0:n5-1))
  end subroutine grow5


  !> Shock-detector weight of the cells of one tile (1 without the detector).
  subroutine tile_beta(blk, p, ni, tile)
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m,     only: obj_space_scheme
    use ICE_Lib_Shock_Detector, only: SD_rho
    implicit none
    type(ICE_block_type), intent(inout) :: blk
    integer(kind=I4),     intent(in)    :: p, ni
    type(tile_type),      intent(in)    :: tile
    integer(kind=I4) :: i, j, k
    if (obj_space_scheme%SD) then
      do k = tile%k1, tile%k2
      do j = tile%j1, tile%j2
      do i = 1, ni
        blk%cond_phase(p)%beta(i,j,k) = &
          SD_rho(blk%cond_phase(p)%prim(1, i-1:i+1, j-1:j+1, k-1:k+1))
      enddo ; enddo ; enddo
    else
      do k = tile%k1, tile%k2
      do j = tile%j1, tile%j2
      do i = 1, ni
        blk%cond_phase(p)%beta(i,j,k) = 1d0
      enddo ; enddo ; enddo
    end if
  end subroutine tile_beta


  !> The k-faces one tile shares with the next: its seam (deep) or its planes' k-faces (shallow).
  subroutine tile_kfaces(blk, p, nc, ni, nj, nk, tile, t, deep, seam_a, kface_a)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_block_type), intent(inout) :: blk
    integer(kind=I4),     intent(in)    :: p, nc, ni, nj, nk, t
    type(tile_type),      intent(in)    :: tile
    logical,              intent(in)    :: deep
    real(R8),             intent(inout) :: seam_a(:,:,:,:), kface_a(:,:,:,:)
    integer(kind=I4) :: i, j, k
    if (deep) then
      k = tile%k2
      if (k < nk) then
        do j = 1, nj
        do i = 1, ni
          call face_flux(blk, p, nc, 3, i, j, k, seam_a(1:nc,i,j,t))
        enddo ; enddo
      end if
    else
      do k = tile%k1, min(tile%k2, nk-1)
      do j = tile%j1, tile%j2
      do i = 1, ni
        call face_flux(blk, p, nc, 3, i, j, k, kface_a(1:nc,i,j,k))
      enddo ; enddo ; enddo
    end if
  end subroutine tile_kfaces


  !> The gathered fluxes of the cells of tile t: the sweep of compute_flux's last loop.
  subroutine tile_sweep(blk, p, nc, ni, nj, nk, tiles, t, deep, seam, kface, fi, gj, hk, jlo, jhi, klo, khi)
    use ICE_Advanced_Types_m
    implicit none
    type(ICE_block_type), intent(inout) :: blk
    integer(kind=I4),     intent(in)    :: p, nc, ni, nj, nk, t
    type(tile_type),      intent(in)    :: tiles(:)
    logical,              intent(in)    :: deep
    real(R8),             intent(inout) :: seam(:,:,:,:), kface(:,:,:,:)
    real(R8),             intent(inout) :: fi(:,0:), gj(:,:,0:), hk(:,:,:,0:)
    integer(kind=I4),     intent(inout) :: jlo, jhi, klo, khi
    integer(kind=I4) :: i, j, k, sw
    do k = tiles(t)%k1, tiles(t)%k2

      !> Deep block: the k-faces below this plane come from the plane before it
      !> -- the previous tile's seam at the tile's first plane, carried
      !> afterwards -- and those above it are computed here, except at the
      !> tile's last plane, where they are its own seam.
      if (deep) then
        if (k == tiles(t)%k1 .and. k > 1) hk(1:nc,1:ni,1:nj,klo) = seam(1:nc,1:ni,1:nj,t-1)
        if (k == tiles(t)%k2) then
          if (k < nk) hk(1:nc,1:ni,1:nj,khi) = seam(1:nc,1:ni,1:nj,t)
        else
          do j = 1, nj
          do i = 1, ni
            call face_flux(blk, p, nc, 3, i, j, k, hk(1:nc,i,j,khi))
          enddo ; enddo
        end if
      end if

      do j = tiles(t)%j1, tiles(t)%j2

        !> The j-faces below this row: from the row before it, computed at the
        !> tile's first row and carried afterwards.
        if (j == tiles(t)%j1 .and. j > 1) then
          do i = 1, ni
            call face_flux(blk, p, nc, 2, i, j-1, k, gj(1:nc,i,jlo))
          enddo
        end if
        if (j < nj) then
          do i = 1, ni
            call face_flux(blk, p, nc, 2, i, j, k, gj(1:nc,i,jhi))
          enddo
        end if

        !> The i-faces of the line
        do i = 1, ni-1
          call face_flux(blk, p, nc, 1, i, j, k, fi(1:nc,i))
        enddo

        !> The cells of the line, each direction in the odd-then-even face order
        do i = 1, ni
          associate (r => blk%cond_phase(p)%residual(1:nc,i,j,k))
            ! direction 1: face i is r's own upper face, face i-1 its lower one
            if (mod(i,2) == 1) then
              if (i < ni) r = r - fi(1:nc,i)
              if (i > 1)  r = r + fi(1:nc,i-1)
            else
              r = r + fi(1:nc,i-1)
              if (i < ni) r = r - fi(1:nc,i)
            end if
            ! direction 2
            if (mod(j,2) == 1) then
              if (j < nj) r = r - gj(1:nc,i,jhi)
              if (j > 1)  r = r + gj(1:nc,i,jlo)
            else
              r = r + gj(1:nc,i,jlo)
              if (j < nj) r = r - gj(1:nc,i,jhi)
            end if
            ! direction 3
            if (deep) then
              if (mod(k,2) == 1) then
                if (k < nk) r = r - hk(1:nc,i,j,khi)
                if (k > 1)  r = r + hk(1:nc,i,j,klo)
              else
                r = r + hk(1:nc,i,j,klo)
                if (k < nk) r = r - hk(1:nc,i,j,khi)
              end if
            else if (nk > 1) then
              if (mod(k,2) == 1) then
                if (k < nk) r = r - kface(1:nc,i,j,k)
                if (k > 1)  r = r + kface(1:nc,i,j,k-1)
              else
                r = r + kface(1:nc,i,j,k-1)
                if (k < nk) r = r - kface(1:nc,i,j,k)
              end if
            end if
          end associate
        enddo

        ! the faces above this row are the faces below the next one
        if (j < nj) then
          sw = jlo ; jlo = jhi ; jhi = sw
        end if
      enddo

      ! the faces above this plane are the faces below the next one
      if (deep .and. k < nk) then
        sw = klo ; klo = khi ; khi = sw
      end if
    enddo
  end subroutine tile_sweep


  !> Grow a shared plane array to hold at least n1 x n2 x n3 x n4 (never shrunk).
  subroutine grow(a, n1, n2, n3, n4)
    real(R8), allocatable, intent(inout) :: a(:,:,:,:)
    integer, intent(in) :: n1, n2, n3, n4
    if (allocated(a)) then
      if (size(a,1) < n1 .or. size(a,2) < n2 .or. size(a,3) < n3 .or. size(a,4) < n4) deallocate(a)
    end if
    if (.not. allocated(a)) allocate(a(n1, n2, n3, n4))
  end subroutine grow


  !> The tiles of a block for nthreads threads. Deep (nk >= nthreads, or a 2-D block
  !> whose planes outnumber the threads in j): contiguous ranges of whole k-planes.
  !> Otherwise pieces of a plane, about two per thread, so a shallow block still
  !> spreads; a 2-D block (nk = 1) is always cut this way, along j.
  subroutine make_tiles(ni, nj, nk, nc, nthreads, tiles, deep)
    use ICE_Global_m, only: tiles_per_thread
    integer, intent(in) :: ni, nj, nk, nc, nthreads
    type(tile_type), allocatable, intent(out) :: tiles(:)
    logical, intent(out) :: deep
    integer :: t, nt, k, kk, per_plane, rows, j, jj, want

    ! Deep tiling keeps one k-range per thread: finer ranges cost 10-41 % from one
    ! socket of 20 threads to 12 x 20 (more seams, less reuse; jobs 303570/303572);
    ! and it is taken only while the threads' plane buffers together fit the budget
    ! -- beyond it the shallow tiling with dynamically scheduled pieces was faster
    ! (1 x 20 -3 %, 4 x 20 -1.5 %, 1 x 80 -18 % on the sweep, job 303578). The
    ! shallow pieces are tiles_per_thread per thread.
    want = nthreads * tiles_per_thread        ! shallow tiles asked for (ICE_TILES_PER_THREAD)
    deep = (nk >= nthreads) .and. &
           (int(nthreads, 8) * 2_8 * int(nc, 8) * int(ni, 8) * int(nj, 8) * 8_8 <= plane_buffer_budget)
    if (deep) then
      nt = min(nthreads, nk)
      allocate(tiles(nt))
      kk = 0
      do t = 1, nt
        k = (nk - kk) / (nt - t + 1)               ! planes left, shared evenly
        tiles(t) = tile_type(kk + 1, kk + k, 1, nj)
        kk = kk + k
      end do
    else
      per_plane = max(1, (2*want + nk - 1) / nk)     ! pieces per plane
      per_plane = min(per_plane, nj)
      rows = (nj + per_plane - 1) / per_plane       ! rows per piece
      nt = 0
      do k = 1, nk
        do j = 1, nj, rows
          nt = nt + 1
        end do
      end do
      allocate(tiles(nt))
      t = 0
      do k = 1, nk
        do j = 1, nj, rows
          t = t + 1
          jj = min(j + rows - 1, nj)
          tiles(t) = tile_type(k, k, j, jj)
        end do
      end do
    end if
  end subroutine make_tiles


  !> The flux through the face between cell (i,j,k) and its upper neighbour in
  !> direction d, times the area: the arithmetic of the scatter form, unchanged.
  !> The i stencil is a contiguous slab of prim and goes in as it is; the j and k
  !> stencils are gathered into four columns first, which is the copy the
  !> compiler made for the scatter form's strided sections.
  subroutine face_flux(blk, p, nc, d, i, j, k, flux)
    use ICE_Advanced_Types_m
    use ICE_Global_m,       only: mat_of, solid_of
    use ICE_Config_Types_m, only: obj_condensed
    implicit none
    type(ICE_block_type), intent(in)  :: blk
    integer(kind=I4),     intent(in)  :: p, nc, d, i, j, k
    real(kind=R8),        intent(out) :: flux(nc)
    real(kind=R8) :: prim(nc,-1:2)

    select case (d)
    case (1)
      call compute_flux_ (blk%cond_phase(p)%prim(1:nc,i-1:i+2,j,k),                            &
                          [blk%dl(i-1,j,k)%c(1), blk%dl(i,j,k)%c(1),                           &
                           blk%dl(i+1,j,k)%c(1), blk%dl(i+2,j,k)%c(1)],                        &
                          blk%dir(1)%f(i,j,k)%N, blk%dir(1)%f(i,j,k)%A, flux, nc,              &
                          blk%cond_phase(p)%beta(i,j,k), obj_condensed(mat_of(p)), solid_of(p))
    case (2)
      prim(1:nc,-1:2) = blk%cond_phase(p)%prim(1:nc,i,j-1:j+2,k)
      call compute_flux_ (prim,                                                                &
                          [blk%dl(i,j-1,k)%c(2), blk%dl(i,j,k)%c(2),                           &
                           blk%dl(i,j+1,k)%c(2), blk%dl(i,j+2,k)%c(2)],                        &
                          blk%dir(2)%f(i,j,k)%N, blk%dir(2)%f(i,j,k)%A, flux, nc,              &
                          blk%cond_phase(p)%beta(i,j,k), obj_condensed(mat_of(p)), solid_of(p))
    case (3)
      prim(1:nc,-1:2) = blk%cond_phase(p)%prim(1:nc,i,j,k-1:k+2)
      call compute_flux_ (prim,                                                                &
                          [blk%dl(i,j,k-1)%c(3), blk%dl(i,j,k)%c(3),                           &
                           blk%dl(i,j,k+1)%c(3), blk%dl(i,j,k+2)%c(3)],                        &
                          blk%dir(3)%f(i,j,k)%N, blk%dir(3)%f(i,j,k)%A, flux, nc,              &
                          blk%cond_phase(p)%beta(i,j,k), obj_condensed(mat_of(p)), solid_of(p))
    end select
  end subroutine face_flux


  subroutine compute_flux_ (prim, length, normal, area, flux, nvar, beta, mat, solid)
    use ICE_Lib_Riemann
    use ICE_Config_Types_m, only: condensed_phase_t
    implicit none
    real(kind=R8), dimension(1:nvar,-1:2), intent(in)    :: prim
    real(kind=R8), dimension(-1:2),        intent(in)    :: length
    real(kind=R8), dimension(3),           intent(in)    :: normal
    real(kind=R8),                         intent(in)    :: area
    real(kind=R8), dimension(1:nvar),      intent(out)   :: flux
    integer(kind=I4),                      intent(in)    :: nvar
    real(kind=R8),                         intent(in)    :: beta
    type(condensed_phase_t),               intent(in)    :: mat
    logical,                               intent(in)    :: solid

    real(kind=R8) :: dl0, dl1, dl2, dll, dlr
    real(kind=R8), dimension(size(prim(:,1))) :: prim_1, prim_4

    dl0 = 0.5_R8 * (length(-1)+length(0))
    dl1 = 0.5_R8 * (length(0)+length(1))
    dl2 = 0.5_R8 * (length(1)+length(2))
    dll = 0.5d0 * length(0)
    dlr = 0.5d0 * length(1)

    call state_reconstruction (prim(:,-1), prim(:,0), prim(:,1), prim(:,2), &
                               dl0, dl1, dl2, dll, dlr, prim_1, prim_4, beta, mat, solid)

    flux = riemann (prim_1, prim_4, normal, mat)
    flux = flux * area

  end subroutine compute_flux_



end module ICE_Mod_Fluxes
