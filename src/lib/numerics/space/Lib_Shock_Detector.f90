module ICE_Lib_Shock_Detector
  use iso_fortran_env, only: R8 => real64

  implicit none
  private
  public :: SD_rho

  real(R8), parameter, private :: DELTA = 20.d0

contains

  pure function SD_rho(rho_matrix) result(beta)
    implicit none
    real(R8), intent(in) :: rho_matrix(0:2, 0:2, 0:2)
    real(R8)             :: beta
    real(R8) :: kp1, kp2, kp3, sw, phi
    real(R8) :: r, r1, r2, r3, r4, r5, r6

    r  = rho_matrix(1,1,1)
    r1 = rho_matrix(0,1,1) ; r2 = rho_matrix(2,1,1)
    r3 = rho_matrix(1,0,1) ; r4 = rho_matrix(1,2,1)
    r5 = rho_matrix(1,1,0) ; r6 = rho_matrix(1,1,2)

    ! Jameson-type sensor on density
    kp1 = abs( (r2 - 2d0*r + r1) / (r2 + 2d0*r + r1 + 1d-30) )
    kp2 = abs( (r4 - 2d0*r + r3) / (r4 + 2d0*r + r3 + 1d-30) )
    kp3 = abs( (r6 - 2d0*r + r5) / (r6 + 2d0*r + r5 + 1d-30) )
    sw = max(kp1, kp2, kp3)

    if (sw < 1d0/DELTA) then
      phi = max(sw/DELTA, 0d0)
      beta = 1d0 - tanh(10d0 * phi*phi*phi)
    else
      beta = 0d0
    end if

  end function SD_rho

end module ICE_Lib_Shock_Detector
