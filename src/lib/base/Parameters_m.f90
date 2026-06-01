module ICE_Parameters_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

  integer, parameter  :: hlen = 1024
  integer, parameter  :: llen = 256
  integer, parameter  :: clen = 16
  character(len=clen) :: codename = 'ICE'

  real(R8), parameter :: pi      = 3.14159265358979323846_R8
  real(R8), parameter :: sigma_SB = 5.67d-8

  integer, dimension(6,3), parameter :: guide = reshape( [ 1, 0, 0, &
                                                           -1, 0, 0, &
                                                            0, 1, 0, &
                                                            0,-1, 0, &
                                                            0, 0, 1, &
                                                            0, 0,-1  ] , shape(guide), order=[2,1] )

end module ICE_Parameters_m
