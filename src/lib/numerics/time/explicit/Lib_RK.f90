module ICE_Lib_RK
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: RK_stage
  public :: RK_accumulate

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


  !> The running sum of a quantity of the stages, weighted as RK_stage weighs their residuals: since
  !> U_irk - U_0 = c_irk (U_(irk-1) - U_0 + residual_irk), after the last stage it holds the stages' values
  !> with the scheme's weights (1; 1/2, 1/2; 1/6, 1/6, 2/3), what the update integrated over the step. The first stage
  !> starts the sum, as every scheme of RKC has c_1 = 1.
  elemental function RK_accumulate(irk, n_rk, acc, stage) result(newacc)
    implicit none
    integer(I4), intent(in) :: irk, n_rk
    real(R8),    intent(in) :: acc, stage
    real(R8)                :: newacc

    if (irk == 1) then
      newacc = stage
    else
      newacc = (acc + stage) * RKC(n_rk, irk)
    end if

  end function RK_accumulate

end module ICE_Lib_RK
