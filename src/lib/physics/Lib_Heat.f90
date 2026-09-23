!>@brief Nusselt-number correlations Nu(Re, Pr, Ma) selected by keyword.
!> Same shape as Lib_Drag: the correlations are pure functions and the choice
!> travels as an integer, so nothing mutable is shared between threads.
module ICE_Lib_Heat
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private

  public :: assign_heat
  public :: heat

  real(kind=R8), parameter :: toll = 1.e-20_R8

contains

  !> Maps the [ICE-Scheme] heat keyword to the selector used by heat().
  subroutine assign_heat(heat_word, heatSelect)
    implicit none
    character(len=*), intent(in)  :: heat_word
    integer(kind=I4), intent(out) :: heatSelect

    select case (heat_word)
    case ('Stokes')
      heatSelect = 1
    case ('JAXA1')
      heatSelect = 2
    case ('JAXA2')
      heatSelect = 3
    case ('JAXA3')
      heatSelect = 4
    case ('Chang')
      heatSelect = 5
    case ('Ranz-Marshall')
      heatSelect = 6
    case ('Kavanau-Drake')
      heatSelect = 7
    case default
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
      error stop 'ICE: unknown heat model'
    end select

  end subroutine assign_heat


  !> Dispatches to the selected Nusselt correlation.
  pure function heat(Re,Pr,Ma,heatSelect) result(Nu)
    implicit none
    integer(kind=I4), intent(in) :: heatSelect
    real(kind=R8),    intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    select case (heatSelect)
    case (1) !> Stokes constant heat model
      Nu = 2._R8
    case (2)
      Nu = heat_JAXA_1(Re)
    case (3)
      Nu = heat_JAXA_2(Re,Pr)
    case (4)
      Nu = heat_JAXA_3(Re,Pr,Ma)
    case (5)
      Nu = heat_Chang(Re,Pr)
    case (6)
      Nu = heat_Ranz_Marshall(Re,Pr)
    case (7)
      Nu = heat_Kavanau_Drake(Re,Pr,Ma)
    end select

  end function heat


  !> JAXA report 1 heat model: tends to 0, not to 2, as the slip vanishes.
  pure function heat_JAXA_1(Re) result(Nu)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(kind=R8) :: Nu

    Nu = 2.5_R8*Re**0.15_R8 + 0.04_R8*Re

  end function heat_JAXA_1


  !> JAXA report 2 heat model
  pure function heat_JAXA_2(Re,Pr) result(Nu)
    implicit none
    real(kind=R8), intent(in) :: Re, Pr
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.37_R8*Re**0.6_R8*Pr**(1._R8/3._R8)

  end function heat_JAXA_2


  !> JAXA report 3 heat model: rarefaction correction in Ma/(Re Pr).
  pure function heat_JAXA_3(Re,Pr,Ma) result(Nu)
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = ((2._R8+0.645_R8*Re**0.5_R8*Pr**(1._R8/3._R8))**(-1._I4) + 3.42_R8*Ma/(Re*Pr+toll))**(-1._I4)

  end function heat_JAXA_3


  !> Chang heat model
  pure function heat_Chang(Re,Pr) result(Nu)
    implicit none
    real(kind=R8), intent(in) :: Re, Pr
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.459_R8*Re**0.55_R8*Pr**(1._R8/3._R8)

  end function heat_Chang


  !> Ranz-Marshall heat model
  pure function heat_Ranz_Marshall(Re,Pr) result(Nu)
    implicit none
    real(kind=R8), intent(in) :: Re, Pr
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.6_R8*Re**0.5_R8*Pr**(1._R8/3._R8)

  end function heat_Ranz_Marshall


  !> Kavanau-Drake heat model: the Chang plateau with a rarefaction denominator.
  pure function heat_Kavanau_Drake(Re,Pr,Ma) result(Nu)
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.459_R8*Re**0.55_R8*Pr**0.33_R8
    Nu = Nu / (1._R8 + 3.42_R8*Ma/(Re*Pr+toll)*Nu)

  end function heat_Kavanau_Drake

end module ICE_Lib_Heat
