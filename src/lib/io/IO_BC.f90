module ICE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Global_m,         only: ICE_phase_prefix
  use ICE_Advanced_Types_m, only: ICE_domain_type, ICE_bc_type

  implicit none
  private
  public :: Setup_BC, Print_BC_Summary

  integer :: nconn = 0, nsym = 0, nio = 0, next = 0

contains


  subroutine Setup_BC(grid)
    implicit none
    type(ICE_domain_type), intent(inout), target :: grid
    type(ICE_bc_type), dimension(:), pointer     :: bc
    integer(I4), dimension(:,:), pointer         :: nbc
    integer(I4) :: i, unitfile, ios
    real(R8)    :: injtype_r

    bc  => grid%bc
    nbc => grid%n_bf
    nconn = 0; nsym = 0; nio = 0; next = 0

    open(newunit=unitfile, file='INPUT/'//trim(ICE_phase_prefix)//'bc.txt', &
         status='old', iostat=ios)
    if (ios /= 0) error stop 'BC file not found'

    nbc = 0
    do i = 1, size(bc)

      read(unitfile,*,iostat=ios) &
        bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%mat, bc(i)%p, bc(i)%type
      if (ios /= 0) write(*,*) '  Error reading BC file at entry', i

      nbc(bc(i)%b, bc(i)%f) = nbc(bc(i)%b, bc(i)%f) + 1

      select case (bc(i)%type)
        case (1)
          nconn = nconn + 1
          read(unitfile,*,iostat=ios) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
            bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
        case (3)
          nsym = nsym + 1
        case (4, 14)
          nio = nio + 1
          read(unitfile,*,iostat=ios) &
            injtype_r, bc(i)%massflux, bc(i)%velocity, &
            bc(i)%alpha, bc(i)%beta, bc(i)%temperature, bc(i)%radius
          bc(i)%injtype = nint(injtype_r)
        case (11)
          next = next + 1
      end select

      if (ios /= 0) then
        write(*,*) '  Error in BC file'
        error stop
      end if

    end do

    close(unitfile)
    call Print_BC_Summary()

  end subroutine Setup_BC


  subroutine Print_BC_Summary()
    implicit none
    write(*,*)
    write(*,'(A)') ' Boundary Conditions:'
    if (nconn > 0) write(*,'(A,T35,I0)') '   Connection',    nconn
    if (nsym  > 0) write(*,'(A,T35,I0)') '   Symmetry',      nsym
    if (nio   > 0) write(*,'(A,T35,I0)') '   Inflow',        nio
    if (next  > 0) write(*,'(A,T35,I0)') '   Extrapolation', next
    write(*,*)
  end subroutine Print_BC_Summary


end module ICE_IO_BC
