module ICE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block
  use ICE_Lib_Reconstruction, only: state_reconstruction, bad_recon

  implicit none
  private
  public :: compute_flux, state_reconstruction, bad_recon

  !> The tile of one thread: the cells (i = 1..ni, j = j1..j2, k = k1..k2) of a block.
  !> Deep blocks are cut into ranges of whole k-planes; shallow ones (fewer planes
  !> than threads) into pieces of a plane, so that a thin block still spreads over
  !> the threads. Every thread builds the same list from the block's dimensions.
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
    !$ use omp_lib,             only: omp_get_num_threads, omp_get_thread_num
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4), intent(in)  :: p
    integer(kind=I4) :: b, i, j, k, t, nt, nthreads, ni, nj, nk, nc, ni_max, nj_max
    logical          :: deep
    type(tile_type), allocatable :: tiles(:)
    !> Per-thread face buffers, sized once per region to the largest local block:
    !> the i-faces of a line, the j-faces below (jlo) and above (jhi) a row, the
    !> k-faces below (klo) and above (khi) a plane (deep blocks only). Moving up a
    !> row or a plane swaps the index: the faces above become the faces below.
    real(R8), allocatable :: fi(:,:), gj(:,:,:), hk(:,:,:,:)
    integer(kind=I4) :: jlo, jhi, klo, khi, sw

    nc = ncond(p)
    bad_recon = .false.
    ni_max = 0 ; nj_max = 0
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      ni_max = max(ni_max, grid%blk(b)%dim(1)) ; nj_max = max(nj_max, grid%blk(b)%dim(2))
    end do
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(b, i, j, k, t, nt, nthreads, ni, nj, nk, deep, tiles, &
    !$OMP                                  fi, gj, hk, jlo, jhi, klo, khi, sw)
    nthreads = 1
    !$ nthreads = omp_get_num_threads()
    allocate(fi(nc, 0:ni_max), gj(nc, ni_max, 0:1), hk(nc, ni_max, nj_max, 0:1))
    jlo = 0 ; jhi = 1 ; klo = 0 ; khi = 1
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      ni = grid%blk(b)%dim(1) ; nj = grid%blk(b)%dim(2) ; nk = grid%blk(b)%dim(3)
      call make_tiles(ni, nj, nk, nthreads, tiles, deep)
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

      !> The k-faces shared between tiles, once, by the owner of the cell below
      !> each: the top faces of every tile's last plane (deep block), or every
      !> k-face of the block (shallow block).
      if (deep) then
        !$OMP SINGLE
        if (allocated(seam)) then
          if (size(seam,1) < nc .or. size(seam,2) < ni .or. size(seam,3) < nj .or. size(seam,4) < nt) &
            deallocate(seam)
        end if
        if (.not. allocated(seam)) allocate(seam(nc, ni, nj, nt))
        !$OMP END SINGLE
        !$OMP DO SCHEDULE (STATIC)
        do t = 1, nt
          k = tiles(t)%k2
          if (k < nk) then
            do j = 1, nj
            do i = 1, ni
              call face_flux(grid%blk(b), p, nc, 3, i, j, k, seam(1:nc,i,j,t))
            enddo ; enddo
          end if
        enddo
        !$OMP END DO
      else if (nk > 1) then
        !$OMP SINGLE
        if (allocated(kface)) then
          if (size(kface,1) < nc .or. size(kface,2) < ni .or. size(kface,3) < nj .or. size(kface,4) < nk) &
            deallocate(kface)
        end if
        if (.not. allocated(kface)) allocate(kface(nc, ni, nj, nk))
        !$OMP END SINGLE
        !$OMP DO SCHEDULE (STATIC)
        do t = 1, nt
          do k = tiles(t)%k1, min(tiles(t)%k2, nk-1)
          do j = tiles(t)%j1, tiles(t)%j2
          do i = 1, ni
            call face_flux(grid%blk(b), p, nc, 3, i, j, k, kface(1:nc,i,j,k))
          enddo ; enddo ; enddo
        enddo
        !$OMP END DO
      end if

      !$OMP DO SCHEDULE (STATIC)
      do t = 1, nt
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
                call face_flux(grid%blk(b), p, nc, 3, i, j, k, hk(1:nc,i,j,khi))
              enddo ; enddo
            end if
          end if

          do j = tiles(t)%j1, tiles(t)%j2

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


  !> The tiles of a block for nthreads threads. Deep (nk >= nthreads, or a 2-D block
  !> whose planes outnumber the threads in j): contiguous ranges of whole k-planes.
  !> Otherwise pieces of a plane, about two per thread, so a shallow block still
  !> spreads; a 2-D block (nk = 1) is always cut this way, along j.
  subroutine make_tiles(ni, nj, nk, nthreads, tiles, deep)
    integer, intent(in) :: ni, nj, nk, nthreads
    type(tile_type), allocatable, intent(out) :: tiles(:)
    logical, intent(out) :: deep
    integer :: t, nt, k, kk, per_plane, rows, j, jj

    deep = (nk >= nthreads)
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
      per_plane = max(1, (2*nthreads + nk - 1) / nk) ! pieces per plane
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
