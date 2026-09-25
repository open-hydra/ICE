module ICE_Lib_Model
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Lib_MK
  use ICE_Lib_IG
  use ICE_Lib_AG
  use ICE_Lib_Drag
  use ICE_Lib_Heat
  use ICE_Config_Types_m, only: obj_time_scheme, condensed_phase_t
  implicit none
  private
  public :: assign_prim_2_cons
  public :: assign_cons_2_prim
  public :: assign_check_prim
  public :: assign_pressure_make
  public :: assign_sound_make
  public :: assign_wavespeed_make
  public :: assign_flux_make
  public :: assign_source_make
  public :: assign_all

  !> Concrete "model" procedure pointing to one of the function realizations
  procedure(prim_2_cons_if), pointer, public   :: prim_2_cons
  procedure(cons_2_prim_if), pointer, public   :: cons_2_prim
  procedure(check_prim_if), pointer, public    :: check_prim
  procedure(pressure_make_if), pointer, public :: pressure_make
  procedure(sound_make_if), pointer, public    :: sound_make
  procedure(wavespeed_make_if), pointer, public :: wavespeed_make
  procedure(flux_make_if), pointer, public     :: flux_make
  procedure(source_make_if), pointer, public   :: source_make

  !> Abstract interface relative to the "model" procedure
  abstract interface
  function prim_2_cons_if(prim, mat) result(cons)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    import :: condensed_phase_t
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: cons(size(prim))
  end function prim_2_cons_if

  function cons_2_prim_if(cons, mat) result(prim)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    import :: condensed_phase_t
    implicit none
    real(kind=R8), intent(in) :: cons(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: prim(size(cons))
  end function cons_2_prim_if

  function check_prim_if(prim) result(prim_status)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    logical                   :: prim_status
  end function check_prim_if

  function pressure_make_if(prim)  result(pressure)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: pressure
  end function pressure_make_if

  function sound_make_if(prim)  result(sound)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: sound
  end function sound_make_if

  function wavespeed_make_if(prim, normal) result(speed)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    real(kind=R8)             :: speed
  end function wavespeed_make_if

  function flux_make_if(prim, normal, mat) result(flux)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    import :: condensed_phase_t
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim))
  end function flux_make_if

  function source_make_if(prim,force,mat) result(source)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    import :: condensed_phase_t
    implicit none
    real(kind=R8), intent(in) :: prim(:), force(6)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: source(size(prim))
  end function source_make_if
  end interface


contains


  subroutine assign_prim_2_cons(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      prim_2_cons => prim_2_cons_MK
    case ('IG')
      prim_2_cons => prim_2_cons_IG
    case ('AG')
      prim_2_cons => prim_2_cons_AG
    end select

  end subroutine assign_prim_2_cons

  subroutine assign_cons_2_prim(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      cons_2_prim => cons_2_prim_MK
    case ('IG')
      cons_2_prim => cons_2_prim_IG
    case ('AG')
      cons_2_prim => cons_2_prim_AG
    end select

  end subroutine assign_cons_2_prim

  subroutine assign_check_prim(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      check_prim => check_prim_MK
    case ('IG')
      check_prim => check_prim_IG
    case ('AG')
      check_prim => check_prim_AG
    end select

  end subroutine assign_check_prim

  subroutine assign_pressure_make(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      pressure_make => pressure_make_MK
    case ('IG')
      pressure_make => pressure_make_IG
    case ('AG')
      pressure_make => pressure_make_AG
    end select

  end subroutine assign_pressure_make

  subroutine assign_sound_make(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      sound_make => sound_make_MK
    case ('IG')
      sound_make => sound_make_IG
    case ('AG')
      sound_make => sound_make_AG
    end select

  end subroutine assign_sound_make

  subroutine assign_wavespeed_make(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      wavespeed_make => wavespeed_make_MK
    case ('IG')
      wavespeed_make => wavespeed_make_IG
    case ('AG')
      wavespeed_make => wavespeed_make_AG
    end select

  end subroutine assign_wavespeed_make

  subroutine assign_flux_make(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      flux_make => flux_make_MK
    case ('IG')
      flux_make => flux_make_IG
    case ('AG')
      flux_make => flux_make_AG
    end select

  end subroutine assign_flux_make

  subroutine assign_source_make(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    select case (trim(obj_time_scheme%model(phase)))
    case ('MK')
      source_make => source_make_MK
    case ('IG')
      source_make => source_make_IG
    case ('AG')
      source_make => source_make_AG
    end select

  end subroutine assign_source_make

  subroutine assign_all(phase)
    implicit none
    integer(kind=I4), intent(in) :: phase

    call assign_prim_2_cons (phase)
    call assign_cons_2_prim (phase)
    call assign_check_prim (phase)
    call assign_pressure_make (phase)
    call assign_sound_make (phase)
    call assign_wavespeed_make (phase)
    call assign_flux_make (phase)
    call assign_source_make (phase)

  end subroutine assign_all


end module ICE_Lib_Model

  

  