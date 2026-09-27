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


  !> The condensed state from the IC (or the restart field), variable by variable in family order. A solidifying
  !> family's frozen and nucleated fractions are read when the file carries them, else derived from its temperature
  !> by IGLOO's injection rule.
  subroutine setup_cond (grid, IOfield_cond)
    use ICE_Config_Types_m,     only: obj_condensed
    use ICE_Lib_Solidification, only: solidPhaseAtInjection
    implicit none
    type(ICE_domain_type), intent(inout)   :: grid
    type(orion_data), intent(inout) :: IOfield_cond
    integer(kind=I4) :: b, p, v, d, nfile, i, j, k, ph
    logical          :: full

    do b = 1, grid%nb
      nfile = size(IOfield_cond%block(b)%vars, 1)
      if (any(solid_of)) then
        if (nfile /= sum(ncond) .and. nfile /= sum(nbase)) call refuse_ic(nfile)
      elseif (nfile < sum(ncond)) then
        call refuse_ic(nfile)
      endif
      full = nfile >= sum(ncond)
      d = 0
      do p = 1, ngroups
        do v = 1, merge(ncond(p), nbase(p), full)
          d = d + 1
          grid%blk(b)%cond_phase(p)%prim(v,1:grid%blk(b)%dim(1),1:grid%blk(b)%dim(2),1:grid%blk(b)%dim(3)) = &
          IOfield_cond%block(b)%vars(d,:,:,:)
        enddo
        if (solid_of(p) .and. .not. full) then
          associate (prim => grid%blk(b)%cond_phase(p)%prim, mat => obj_condensed(mat_of(p)))
          do k = 1, grid%blk(b)%dim(3) ; do j = 1, grid%blk(b)%dim(2) ; do i = 1, grid%blk(b)%dim(1)
            call solidPhaseAtInjection(prim(nbase(p)-1,i,j,k), mat%Tmelt, mat%Tnuc, ph, prim(nbase(p)+1,i,j,k))
            prim(ncond(p),i,j,k) = prim(nbase(p)+1,i,j,k)
          enddo ; enddo ; enddo
          end associate
        endif
      enddo
    enddo

  contains

    subroutine refuse_ic(n)
      integer(kind=I4), intent(in) :: n
      character(len=16) :: a, b1, b2
      write(a, '(I0)') n ; write(b1, '(I0)') sum(ncond) ; write(b2, '(I0)') sum(nbase)
      if (any(solid_of)) then
        write(*,'(A)') ' [ERROR] [ICE::setup_cond] the initial condition holds '//trim(a)//' variables per cell; expected '// &
          trim(b1)//' (with the frozen and nucleated fractions of the solidifying families) or '//trim(b2)// &
          ' (without them: both from T)'
      else
        write(*,'(A)') ' [ERROR] [ICE::setup_cond] the initial condition holds '//trim(a)//' variables per cell; expected '// &
          trim(b1)
      endif
      error stop 1
    end subroutine refuse_ic

  end subroutine setup_cond

  
end module ICE_Mod_Phase

    
    


