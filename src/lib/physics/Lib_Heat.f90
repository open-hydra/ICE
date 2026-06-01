module ICE_Lib_Heat
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: assign_heat

  !> Concrete "drag" procedure pointing to one of the function realizations
  procedure(heat_if), pointer, public :: heat

  !> Abstract interface relative to the "limiter" procedure
  abstract interface
  pure function heat_if(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu
  end function heat_if
  end interface

  real(kind=R8), parameter :: toll = 1.e-20_R8
  
contains

  subroutine assign_heat(heat_word)
    implicit none
    character(len=*), intent(in) :: heat_word

    select case (heat_word)
    case ('Stokes')
      heat => heat_Stokes
    case ('JAXA1')
      heat => heat_JAXA_1
    case ('JAXA2')
      heat => heat_JAXA_2
    case ('JAXA3')
      heat => heat_JAXA_3
    case ('Chang')
      heat => heat_Chang
    case ('Ranz-Marshall')
      heat => heat_Ranz_Marshall
    case ('Kavanau-Drake')
      heat => heat_Kavanau_Drake
    case default
      write(*,*)
      write(*,*)
      write(*,*) "Wrong heat input ---> "//heat_word
      write(*,*) "Choose one of the following :"
      write(*,*) "- Stokes "
      write(*,*) "- JAXA1 "
      write(*,*) "- JAXA2 "
      write(*,*) "- JAXA3 "
      write(*,*) "- Chang "
      write(*,*) "- Ranz-Marshall "
      write(*,*) "- Kavanau-Drake "
      write(*,*)
      stop
    end select

  end subroutine assign_heat

  !> Stokes constant heat model
  pure function heat_Stokes(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = 2.0_R8

  end function heat_Stokes

  !> Jaxa Report 1 Heat Model
  pure function heat_JAXA_1(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = 2.5_R8*Re**0.15_R8 + 0.04_R8*Re

  end function heat_JAXA_1

  !> Jaxa Report 2 Heat Model
  pure function heat_JAXA_2(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.37_R8*Re**0.6_R8*Pr**(1._R8/3._R8)

  end function heat_JAXA_2

  !> Jaxa Report 3 Heat Model
  pure function heat_JAXA_3(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu  = ((2._R8+0.645_R8*Re**0.5_R8*Pr**(1._R8/3._R8))**(-1._I4) + 3.42_R8*Ma/(Re*Pr))**(-1._I4)

  end function heat_JAXA_3

  !> Chang Heat Model
  pure function heat_Chang(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.459_R8*Re**0.55_R8*Pr**(1._R8/3._R8)

  end function heat_Chang

  !> Ranz Marshall Heat Model
  pure function heat_Ranz_Marshall(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.6_R8*Re**0.5_R8*Pr**(1._R8/3._R8)
    
  end function heat_Ranz_Marshall

  !> Kavanau Drake Heat Model
  pure function heat_Kavanau_Drake(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.459_R8*Re**0.55_R8*Pr**0.33_R8
    Nu = Nu / (1._R8 + 3.42_R8*Ma/(Re*Pr+toll)*Nu)

  end function heat_Kavanau_Drake

end module ICE_Lib_Heat
