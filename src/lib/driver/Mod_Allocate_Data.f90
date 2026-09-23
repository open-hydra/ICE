module ICE_Mod_Allocate_Data
  implicit none
  private
  public :: allocate_data, Allocate_Block

contains

  subroutine allocate_data (grid, coupled, IOfield)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use Lib_ORION_data
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
      call Allocate_Block(grid%blk(b), [ni, nj, nk])

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
      endif

    enddo

    !> Bound
    nbound = nbound * ngroups
    allocate( grid%bc(1:nbound) )
    allocate( grid%n_bf(1:size(IOfield%block),1:6) )

  end subroutine allocate_data


  subroutine Allocate_Block(blk, nijk)
    use ICE_Global_m
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m, only: obj_irs
    implicit none
    type(ICE_block_type), intent(inout) :: blk
    integer(kind=I4),     intent(in)    :: nijk(3)
    integer(kind=I4) :: ni, nj, nk, p

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
      if (obj_irs%enabled) then
        allocate( blk%cond_phase(p)%RS1 (1:ncond(p), 0:ni+1, 0:nj+1, 0:nk+1) )
        allocate( blk%cond_phase(p)%RS2 (1:ncond(p), 0:ni+1, 0:nj+1, 0:nk+1) )
        blk%cond_phase(p)%RS1 = 0._R8
        blk%cond_phase(p)%RS2 = 0._R8
      end if
      blk%cond_phase(p)%beta = 1d0
    enddo

  end subroutine Allocate_Block

end module ICE_Mod_Allocate_Data