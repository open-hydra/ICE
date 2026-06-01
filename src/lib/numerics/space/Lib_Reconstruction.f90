module ICE_Lib_Reconstruction
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: state_reconstruction

contains

  subroutine state_reconstruction(prev, local, next, next2, dl0, dl1, dl2, dll, dlr, priml, primr, beta)
    use ICE_Lib_Limiters
    use ICE_Lib_Model
    implicit none
    real(R8), dimension(:), intent(in)  :: prev, local, next, next2
    real(R8), dimension(:), intent(out) :: priml, primr
    real(R8),               intent(in)  :: dl0, dl1, dl2, dll, dlr
    real(R8),               intent(in)  :: beta

    integer(I4) :: i
    real(R8)    :: slope0, slope1, slope2, slopel, sloper, limval
    logical     :: prim_status_l, prim_status_r

    prim_status_l = .false.
    prim_status_r = .false.
    limval = beta

    ! Piecewise Linear Reconstruction
    if (prev(1)<=1d-6 .or. local(1)<=1d-6 .or. next(1)<=1e-6) then
      priml = local
      primr = next
    else
      do while (.not. prim_status_l .or. .not. prim_status_r)
        do i = 1, size(prev)
          slope0   = (local(i)-prev(i)) / dl0
          slope1   = (next(i)-local(i)) / dl1
          slope2   = (next2(i)-next(i)) / dl2
          slopel   = limiter(slope1, slope0) * limval
          sloper   = limiter(slope2, slope1) * limval
          priml(i) = local(i) + slopel*dll
          primr(i) =  next(i) - sloper*dlr
        enddo
        prim_status_l = check_prim(priml)
        prim_status_r = check_prim(primr)
        limval = 0.5_R8*limval
      end do
    endif

  end subroutine state_reconstruction

end module ICE_Lib_Reconstruction
