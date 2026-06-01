module ICE_Procedures_m
  use ICE_Wrap_Setup
  use ICE_Wrap_Solve
  use ICE_Wrap_Postprocess

  implicit none

  type :: ICE_type

  contains

    procedure, nopass  :: Setup => ICE_setup
    procedure, nopass  :: Solve => ICE_solve
    procedure, nopass  :: Postprocess => ICE_postprocess

  end type ICE_type


end module ICE_Procedures_m