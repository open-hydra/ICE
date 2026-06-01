module ICE_Lib_RK
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: RK_stage

  real(R8), parameter, private :: RKC(3,3) = reshape([ 1._R8,    0._R8,      0._R8,  &
                                                        1._R8, 0.50_R8,      0._R8,  &
                                                        1._R8, 0.25_R8, 2._R8/3._R8 ],&
                                                        shape(RKC), order=[2,1])

contains

  pure function RK_stage(irk, n_rk, cons, cons_old, residual) result(newcons)
    implicit none
    integer(I4),            intent(in) :: irk, n_rk
    real(R8), dimension(:), intent(in) :: cons, cons_old, residual
    real(R8), dimension(size(cons))    :: newcons, residual_

    residual_ = (cons - cons_old + residual) * RKC(n_rk, irk)
    newcons   = cons_old + residual_

  end function RK_stage

end module ICE_Lib_RK
