module ICE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Lib_Reconstruction, only: state_reconstruction, bad_recon

  implicit none
  private
  public :: compute_flux, state_reconstruction, bad_recon

  !> The tile of one thread: the cells (i = 1..ni, j = j1..j2, k = k1..k2) of a
  !> block, a range of k-planes (the r-th of the block's nkr ranges) by a slab of
  !> j-rows. The slab is cut so that a thread's two plane buffers stay within
  !> plane_budget bytes whatever the plane size -- on a 192 x 192 plane with eight
  !> variables a whole plane is 2.4 MB per buffer, and twenty threads of them per
  !> socket overflow the L3 -- and the ranges then bring the tile count up to the
  !> thread count. Every thread builds the same list from the block's dimensions.
  type :: tile_type
    integer :: k1 = 0, k2 = 0, j1 = 0, j2 = 0, r = 0
  end type tile_type

  !> Bytes a thread's pair of plane buffers may take (two of nc x ni x jslab reals).
  integer, parameter :: plane_budget = 524288

  !> A k-face between two k-ranges is read by both: by the tile owning the cell
  !> below it and by the one owning the cell above. It is computed once, by the
  !> owner of the cell below, into a shared plane per range boundary. Grown to the
  !> largest block that needed it; freed with the program.
  real(R8), allocatable, private :: seam(:,:,:,:)

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
    integer(kind=I4) :: b, i, j, k, t, nt, nthreads, ni, nj, nk, nc, ni_max, js_max, nkr, jslab, jj
    !> tiles is each thread's inside the region; the serial pass before it uses its
    !> own tiles0: ifx 2023.1 reports "already allocated" at the first allocation of a
    !> PRIVATE copy whose original was allocated and freed before the region.
    type(tile_type), allocatable :: tiles(:), tiles0(:)
    !> Per-thread face buffers, sized once per region to the largest local block and
    !> the widest slab: the i-faces of a line, the j-faces below (jlo) and above (jhi)
    !> a row, the k-faces below (klo) and above (khi) a plane of the tile's slab.
    !> Moving up a row or a plane swaps the index: the faces above become the faces
    !> below.
    real(R8), allocatable :: fi(:,:), gj(:,:,:), hk(:,:,:,:)
    integer(kind=I4) :: jlo, jhi, klo, khi, sw

    nc = ncond(p)
    bad_recon = .false.
    !> The tiling is decided by the thread count the region will have at most, so
    !> that the shared seam planes every local block needs can be grown here, once,
    !> instead of by a SINGLE (and its barrier) per block.
    nthreads = 1
    !$ nthreads = omp_get_max_threads()
    ni_max = 0 ; js_max = 0
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
      call make_tiles(ni, nj, nk, nc, nthreads, tiles0, nkr, jslab)
      ni_max = max(ni_max, ni) ; js_max = max(js_max, jslab)
      if (nkr > 1) call grow(seam, nc, ni, nj, nkr - 1)
    end do
    if (allocated(tiles0)) deallocate(tiles0)
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(b, i, j, k, t, nt, ni, nj, nk, nkr, jslab, jj, tiles, &
    !$OMP                                  fi, gj, hk, jlo, jhi, klo, khi, sw)
    allocate(fi(nc, 0:ni_max), gj(nc, ni_max, 0:1), hk(nc, ni_max, js_max, 0:1))
    jlo = 0 ; jhi = 1 ; klo = 0 ; khi = 1
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
      call make_tiles(ni, nj, nk, nc, nthreads, tiles, nkr, jslab)
      nt = size(tiles)

      !> The shock-detector weight of every cell, before any face reads it: a face
      !> takes the weight of the cell below it, which may belong to another tile.
      if (obj_space_scheme%SD) then
        !$OMP DO SCHEDULE (STATIC)
        do t = 1, nt
          do k = tiles(t)%k1, tiles(t)%k2
          do j = tiles(t)%j1, tiles(t)%j2
          do i = 1, ni
            grid%blk(b)%cond_phase(p)%beta(i,j,k) = &
              SD_rho(grid%blk(b)%cond_phase(p)%prim(1, i-1:i+1, j-1:j+1, k-1:k+1))
          enddo ; enddo ; enddo
        enddo
        !$OMP END DO
      else
        !$OMP DO SCHEDULE (STATIC)
        do t = 1, nt
          do k = tiles(t)%k1, tiles(t)%k2
          do j = tiles(t)%j1, tiles(t)%j2
          do i = 1, ni
            grid%blk(b)%cond_phase(p)%beta(i,j,k) = 1d0
          enddo ; enddo ; enddo
        enddo
        !$OMP END DO
      end if

      !> The k-faces shared between the k-ranges, once, by the owner of the cell
      !> below each: the top faces of the last plane of every range but the top one,
      !> each tile its own slab of them.
      if (nkr > 1) then
        !$OMP DO SCHEDULE (STATIC)
        do t = 1, nt
          k = tiles(t)%k2
          if (k < nk) then
            do j = tiles(t)%j1, tiles(t)%j2
            do i = 1, ni
              call face_flux(grid%blk(b), p, nc, 3, i, j, k, seam(1:nc,i,j,tiles(t)%r))
            enddo ; enddo
          end if
        enddo
        !$OMP END DO
      end if

      !$OMP DO SCHEDULE (STATIC)
      do t = 1, nt
        do k = tiles(t)%k1, tiles(t)%k2

          !> The k-faces below this plane come from the plane before it -- the
          !> previous range's seam at the tile's first plane, carried afterwards --
          !> and those above it are computed here, except at the range's last
          !> plane, where they are its own seam.
          if (nk > 1) then
            if (k == tiles(t)%k1 .and. k > 1) then
              do j = tiles(t)%j1, tiles(t)%j2
                jj = j - tiles(t)%j1 + 1
                hk(1:nc,1:ni,jj,klo) = seam(1:nc,1:ni,j,tiles(t)%r-1)
              enddo
            end if
            if (k == tiles(t)%k2) then
              if (k < nk) then
                do j = tiles(t)%j1, tiles(t)%j2
                  jj = j - tiles(t)%j1 + 1
                  hk(1:nc,1:ni,jj,khi) = seam(1:nc,1:ni,j,tiles(t)%r)
                enddo
              end if
            else
              do j = tiles(t)%j1, tiles(t)%j2
                jj = j - tiles(t)%j1 + 1
                do i = 1, ni
                  call face_flux(grid%blk(b), p, nc, 3, i, j, k, hk(1:nc,i,jj,khi))
                enddo
              enddo
            end if
          end if

          do j = tiles(t)%j1, tiles(t)%j2
            jj = j - tiles(t)%j1 + 1

            !> The j-faces below this row: from the row before it, computed at the
            !> tile's first row and carried afterwards.
            if (j == tiles(t)%j1 .and. j > 1) then
              do i = 1, ni
                call face_flux(grid%blk(b), p, nc, 2, i, j-1, k, gj(1:nc,i,jlo))
              enddo
            end if
            if (j < nj) then
              do i = 1, ni
                call face_flux(grid%blk(b), p, nc, 2, i, j, k, gj(1:nc,i,jhi))
              enddo
            end if

            !> The i-faces of the line
            do i = 1, ni-1
              call face_flux(grid%blk(b), p, nc, 1, i, j, k, fi(1:nc,i))
            enddo

            !> The cells of the line, each direction in the odd-then-even face order
            do i = 1, ni
              associate (r => grid%blk(b)%cond_phase(p)%residual(1:nc,i,j,k))
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
                if (nk > 1) then
                  if (mod(k,2) == 1) then
                    if (k < nk) r = r - hk(1:nc,i,jj,khi)
                    if (k > 1)  r = r + hk(1:nc,i,jj,klo)
                  else
                    r = r + hk(1:nc,i,jj,klo)
                    if (k < nk) r = r - hk(1:nc,i,jj,khi)
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
          if (nk > 1 .and. k < nk) then
            sw = klo ; klo = khi ; khi = sw
          end if
        enddo
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


  !> Grow a shared plane array to hold at least n1 x n2 x n3 x n4 (never shrunk).
  subroutine grow(a, n1, n2, n3, n4)
    real(R8), allocatable, intent(inout) :: a(:,:,:,:)
    integer, intent(in) :: n1, n2, n3, n4
    if (allocated(a)) then
      if (size(a,1) < n1 .or. size(a,2) < n2 .or. size(a,3) < n3 .or. size(a,4) < n4) deallocate(a)
    end if
    if (.not. allocated(a)) allocate(a(n1, n2, n3, n4))
  end subroutine grow


  !> The tiles of a block for nthreads threads: nkr ranges of whole k-planes by
  !> slabs of jslab rows. The slab keeps a thread's two plane buffers (nc x ni x
  !> jslab reals each) within plane_budget; the ranges bring the tile count to the
  !> thread count, one plane at least each. The list is range-major, so a thread's
  !> contiguous share of it stays inside one range as far as it can.
  subroutine make_tiles(ni, nj, nk, nc, nthreads, tiles, nkr, jslab)
    integer, intent(in) :: ni, nj, nk, nc, nthreads
    type(tile_type), allocatable, intent(out) :: tiles(:)
    integer, intent(out) :: nkr, jslab
    integer :: njs, r, s, k, kk, j, jj

    jslab = max(1, min(nj, plane_budget / (2 * nc * ni * 8)))
    njs   = (nj + jslab - 1) / jslab
    nkr   = max(1, min(nk, (nthreads + njs - 1) / njs))
    allocate(tiles(nkr * njs))
    kk = 0
    do r = 1, nkr
      k = (nk - kk) / (nkr - r + 1)               ! planes left, shared evenly
      do s = 1, njs
        j  = (s - 1) * jslab + 1
        jj = min(j + jslab - 1, nj)
        tiles((r - 1) * njs + s) = tile_type(kk + 1, kk + k, j, jj, r)
      end do
      kk = kk + k
    end do
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
