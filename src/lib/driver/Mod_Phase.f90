module ICE_Mod_Phase
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64, iostat_end
  use ICE_Global_m
  use ICE_Advanced_Types_m
  use Lib_ORION_data
  
  implicit none

contains

  subroutine setup_gas (grid, IOfield_gas, nsc, nrans)
    implicit none
    type(ICE_domain_type), intent(inout)   :: grid
    type(orion_data), intent(inout) :: IOfield_gas
    integer(kind=I4), optional, intent(in) :: nsc, nrans
    integer(kind=I4) :: b

    !> Import gaseous phase variables
    if (present(nsc)) then
      do b = 1, grid%nb
        ! Density
        grid%blk(b)%gas_phase%prim(1,1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        sum(IOfield_gas%block(b)%vars(1:nsc,:,:,:),1)
        ! Velocity (3 components)
        grid%blk(b)%gas_phase%prim(2:4,1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(nsc+1:nsc+3,:,:,:)
        ! Temperature
        grid%blk(b)%gas_phase%prim(5,1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(nsc+nrans+5,:,:,:)
        ! gamma
        grid%blk(b)%gas_phase%gam(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(nsc+nrans+6,:,:,:)
        ! R
        grid%blk(b)%gas_phase%R(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(nsc+nrans+7,:,:,:)
        ! mu
        grid%blk(b)%gas_phase%mu(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(nsc+nrans+8,:,:,:)
        ! k
        grid%blk(b)%gas_phase%k(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(nsc+nrans+9,:,:,:)
      enddo

    else

      do b = 1, grid%nb
        grid%blk(b)%gas_phase%prim(1:5,1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(1:5,:,:,:)
        grid%blk(b)%gas_phase%R(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(6,:,:,:)
        grid%blk(b)%gas_phase%gam(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(7,:,:,:)
        grid%blk(b)%gas_phase%k(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(8,:,:,:)
        grid%blk(b)%gas_phase%mu(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(9,:,:,:)
        grid%blk(b)%gas_phase%dt(1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
        IOfield_gas%block(b)%vars(10,:,:,:)
      enddo

    endif

  end subroutine setup_gas


  subroutine setup_cond (grid, IOfield_cond)
    implicit none
    type(ICE_domain_type), intent(inout)   :: grid
    type(orion_data), intent(inout) :: IOfield_cond
    integer(kind=I4) :: b, p, v, d

    !> Import condensed phase variables
    !if (restart) then
      do b = 1, grid%nb
        p = 1 ; v = 1
        do d = 1, sum(ncond)
          grid%blk(b)%cond_phase(p)%prim(v,1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
          IOfield_cond%block(b)%vars(d,:,:,:)
          v = v + 1
          if (v > ncond(p)) then
            v = 1 ; p = p +1
          endif
        enddo
      enddo

  end subroutine setup_cond

  
end module ICE_Mod_Phase

    
    


