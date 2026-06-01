module ICE_Base_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none


  type :: ICE_tensor_3D_type
    real(R8) :: c(3,3)       !> Metric tensor components.
  end type ICE_tensor_3D_type


  type :: ICE_vector_3D_type
    real(R8) :: c(3)         !> Average cell length components.
  end type ICE_vector_3D_type


  type :: ICE_tensor_3D_R3_type
    real(R8) :: c(3,3,3)     !> Metric tensor components.
  end type ICE_tensor_3D_R3_type


  type :: ICE_f_metrics_type
    real(R8) :: N(3)         !> Unit normal vector.
    real(R8) :: A            !> Interface area.
  end type ICE_f_metrics_type


  type :: ICE_d_metrics_type
    type(ICE_f_metrics_type), allocatable :: f(:,:,:)
  end type ICE_d_metrics_type

  
end module ICE_Base_Types_m
