!> The test of a real value for NaN and infinities by its bits.
module ICE_Finite_m
  use iso_fortran_env, only: R8 => real64

  implicit none
  private
  public :: Not_Finite, Any_Not_Finite

contains

  !> .true. for NaN and infinities, the values with all the exponent bits set. The test reads the bits:
  !> -ffast-math (-ffinite-math-only) lets gfortran fold isnan, ieee_is_nan and ieee_is_finite to constants, and
  !> isnan misses an infinity in every build. Reading the bits is an operation on integers, which raises no
  !> floating-point exception, also in a build that traps the invalid operations.
  elemental logical function Not_Finite(x)
    use iso_fortran_env, only: int64
    implicit none
    real(R8), intent(in) :: x
    integer(int64), parameter :: expo = int(z'7FF', int64)

    Not_Finite = iand(shiftr(transfer(x, 0_int64), 52), expo) == expo

  end function Not_Finite

  !> .true. when any value of x is NaN or infinite, in one call: the per-cell and per-face state checks use it, as
  !> Not_Finite applied to an array from another file costs one call per value (it is not inlined there).
  pure logical function Any_Not_Finite(x)
    implicit none
    real(R8), intent(in) :: x(:)

    Any_Not_Finite = any(Not_Finite(x))

  end function Any_Not_Finite

end module ICE_Finite_m
