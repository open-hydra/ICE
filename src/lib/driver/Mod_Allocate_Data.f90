module ICE_Mod_Allocate_Data
  implicit none
  private
  public :: allocate_data, Allocate_Block

contains

  subroutine allocate_data (grid, coupled, IOfield)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use Lib_ORION_data
    use ICE_Mod_ThreadGroups, only: block_group
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    logical, intent(in)                  :: coupled
    type(orion_data), intent(inout)      :: IOfield
    integer(kind=I4) :: ni, nj, nk, nbound
    integer(kind=I4) :: b, p

    nbound = 0
    allocate( grid%blk(1:size(IOfield%block)) )
    grid%nb = size(IOfield%block)

    do b = 1, size(IOfield%block)

      ni = IOfield%block(b)%Ni ; nj = IOfield%block(b)%Nj ; nk = IOfield%block(b)%Nk

      grid%blk(b)%dim = [ ni, nj, nk ]
      call Allocate_Block(grid%blk(b), [ni, nj, nk], b)

      !> Bound
      nbound = nbound + 2*nj*nk + 2*ni*nk + 2*nj*ni

      !> Gaseous phase
      if (coupled) then
        allocate( grid%blk(b)%gas_phase )
        allocate( grid%blk(b)%gas_phase%prim (1:5,  1:ni, 1:nj, 1:nk) )
        allocate( grid%blk(b)%gas_phase%R    (1:ni, 1:nj, 1:nk) )
        allocate( grid%blk(b)%gas_phase%gam  (1:ni, 1:nj, 1:nk) )
        allocate( grid%blk(b)%gas_phase%k    (1:ni, 1:nj, 1:nk) )
        allocate( grid%blk(b)%gas_phase%mu   (1:ni, 1:nj, 1:nk) )
        allocate( grid%blk(b)%gas_phase%dt   (1:ni, 1:nj, 1:nk) )
        if (block_group(b) >= 0) call touch_gas_by_group(grid%blk(b), nk, block_group(b))
      endif

    enddo

    !> Bound
    nbound = nbound * ngroups
    allocate( grid%bc(1:nbound) )
    allocate( grid%n_bf(1:size(IOfield%block),1:6) )

  end subroutine allocate_data


  subroutine Allocate_Block(blk, nijk, b)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_irs
    use ICE_Mod_ThreadGroups, only: block_group
    implicit none
    type(ICE_block_type), intent(inout) :: blk
    integer(kind=I4),     intent(in)    :: nijk(3)
    integer(kind=I4), optional, intent(in) :: b   !< Mesh block number: with thread groups, its group first-touches it
    integer(kind=I4) :: ni, nj, nk, p, k, g

    g = -1
    if (present(b)) g = block_group(b)

    ni = nijk(1) ; nj = nijk(2) ; nk = nijk(3)

    allocate( blk%node( -1:ni+1, -1:nj+1, -1:nk+1 ) )
    allocate( blk%M  (  0:ni+1,   0:nj+1,  0:nk+1 ) )
    allocate( blk%dl (  0:ni+1,   0:nj+1,  0:nk+1 ) )
    allocate( blk%vol(  1:ni,     1:nj,    1:nk    ) )
    allocate( blk%dir(1)%f( 0:ni, 1:nj, 1:nk ) )
    allocate( blk%dir(2)%f( 1:ni, 0:nj, 1:nk ) )
    allocate( blk%dir(3)%f( 1:ni, 1:nj, 0:nk ) )

    allocate( blk%cond_phase(1:ngroups) )
    do p = 1, ngroups
      allocate( blk%cond_phase(p)%dt       (1:ni, 1:nj, 1:nk) )
      allocate( blk%cond_phase(p)%tau      (1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc) )
      allocate( blk%cond_phase(p)%beta     (1:ni, 1:nj, 1:nk) )
      allocate( blk%cond_phase(p)%prim     (1:ncond(p), 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc) )
      allocate( blk%cond_phase(p)%prim_old (1:ncond(p), 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc) )
      allocate( blk%cond_phase(p)%source   (1:ncond(p), 1:ni, 1:nj, 1:nk) )
      allocate( blk%cond_phase(p)%residual (1:ncond(p), 1:ni, 1:nj, 1:nk) )

      !> First touch by the threads, plane by plane with the static partition the
      !> per-cell phases use (and, per socket, the flux kernel's k-ranges): a page
      !> lands on the socket of the thread that will sweep it. Written serially, a
      !> whole block sits on one socket and a rank of 80 threads reads three quarters
      !> of it across the interconnect. The values are the ones written before.
      if (g < 0) then
      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do k = 1-gc, nk+gc
        blk%cond_phase(p)%prim    (:,:,:,k) = ieee_value(1._R8, ieee_quiet_nan)
        blk%cond_phase(p)%prim_old(:,:,:,k) = ieee_value(1._R8, ieee_quiet_nan)
        blk%cond_phase(p)%tau     (:,:,k)   = ieee_value(1._R8, ieee_positive_inf)
      end do
      !$OMP END PARALLEL DO
      end if

      if (obj_irs%enabled) then
        allocate( blk%cond_phase(p)%RS1 (1:ncond(p), 0:ni+1, 0:nj+1, 0:nk+1) )
        allocate( blk%cond_phase(p)%RS2 (1:ncond(p), 0:ni+1, 0:nj+1, 0:nk+1) )
        blk%cond_phase(p)%RS1 = 0._R8
        blk%cond_phase(p)%RS2 = 0._R8
      end if
      if (g < 0) then
      !$OMP PARALLEL DO SCHEDULE(STATIC)
      do k = 1, nk
        blk%cond_phase(p)%beta    (:,:,k)   = 1d0
        blk%cond_phase(p)%source  (:,:,:,k) = 0._R8
        blk%cond_phase(p)%residual(:,:,:,k) = 0._R8
      end do
      !$OMP END PARALLEL DO
      end if
    enddo

    !> A block owned by a thread group: every array of it first written by that group's
    !> threads, so its pages land on the group's socket -- the condensed fields with the
    !> values above, the geometry and dt with zeros that set-up overwrites.
    if (g >= 0) call touch_by_group(blk, ni, nj, nk, g)

  end subroutine Allocate_Block


  !> The first write of a group-owned block, plane by plane over the group's threads.
  subroutine touch_by_group(blk, ni, nj, nk, g)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use ICE_Global_m
    use ICE_Base_Types_m,   only: ICE_vector_3D_type, ICE_tensor_3D_type, ICE_f_metrics_type
    use ICE_Advanced_Types_m
    use ICE_Mod_ThreadGroups, only: thread_group, split_range
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_block_type), intent(inout) :: blk
    integer(kind=I4),     intent(in)    :: ni, nj, nk, g
    integer(kind=I4) :: tid, gg, m, s, k, k1, k2, p, klo, khi
    real(R8) :: qnan, pinf

    qnan = ieee_value(1._R8, ieee_quiet_nan)
    pinf = ieee_value(1._R8, ieee_positive_inf)
    klo = min(-1, 1-gc) ; khi = max(nk+1, nk+gc)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid, gg, m, s, k, k1, k2, p)
    tid = 0
    !$ tid = omp_get_thread_num()
    call thread_group(tid, gg, m, s)
    if (gg == g) then
      call split_range(klo, khi, m, s, k1, k2)
      do k = k1, k2
        if (k >= -1 .and. k <= nk+1) blk%node(:,:,k) = ICE_vector_3D_type(0._R8)
        if (k >=  0 .and. k <= nk+1) then
          blk%M (:,:,k) = ICE_tensor_3D_type(0._R8)
          blk%dl(:,:,k) = ICE_vector_3D_type(0._R8)
        end if
        if (k >= 0 .and. k <= nk) blk%dir(3)%f(:,:,k) = ICE_f_metrics_type(0._R8, 0._R8)
        if (k >= 1 .and. k <= nk) then
          blk%vol(:,:,k)       = 0._R8
          blk%dir(1)%f(:,:,k)  = ICE_f_metrics_type(0._R8, 0._R8)
          blk%dir(2)%f(:,:,k)  = ICE_f_metrics_type(0._R8, 0._R8)
        end if
        do p = 1, ngroups
          if (k >= 1-gc .and. k <= nk+gc) then
            blk%cond_phase(p)%prim    (:,:,:,k) = qnan
            blk%cond_phase(p)%prim_old(:,:,:,k) = qnan
            blk%cond_phase(p)%tau     (:,:,k)   = pinf
          end if
          if (k >= 1 .and. k <= nk) then
            blk%cond_phase(p)%beta    (:,:,k)   = 1d0
            blk%cond_phase(p)%source  (:,:,:,k) = 0._R8
            blk%cond_phase(p)%residual(:,:,:,k) = 0._R8
            blk%cond_phase(p)%dt      (:,:,k)   = 0._R8
          end if
        end do
      end do
    end if
    !$OMP END PARALLEL
  end subroutine touch_by_group


  !> The gas fields of a group-owned block, first written by its group (setup_gas copies
  !> the carrier's values into them afterwards, serially, onto the pages placed here).
  subroutine touch_gas_by_group(blk, nk, g)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Mod_ThreadGroups, only: thread_group, split_range
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    type(ICE_block_type), intent(inout) :: blk
    integer(kind=I4),     intent(in)    :: nk, g
    integer(kind=I4) :: tid, gg, m, s, k, k1, k2

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid, gg, m, s, k, k1, k2)
    tid = 0
    !$ tid = omp_get_thread_num()
    call thread_group(tid, gg, m, s)
    if (gg == g) then
      call split_range(1, nk, m, s, k1, k2)
      do k = k1, k2
        blk%gas_phase%prim(:,:,:,k) = 0._R8
        blk%gas_phase%R  (:,:,k) = 0._R8
        blk%gas_phase%gam(:,:,k) = 0._R8
        blk%gas_phase%k  (:,:,k) = 0._R8
        blk%gas_phase%mu (:,:,k) = 0._R8
        blk%gas_phase%dt (:,:,k) = 0._R8
      end do
    end if
    !$OMP END PARALLEL
  end subroutine touch_gas_by_group

end module ICE_Mod_Allocate_Data