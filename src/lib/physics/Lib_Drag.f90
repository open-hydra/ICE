module ICE_Lib_Drag
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: assign_drag

  !> Concrete "drag" procedure pointing to one of the function realizations
  procedure(drag_if), pointer, public :: drag

  !> Abstract interface relative to the "drag" procedure
  abstract interface
  function drag_if(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd
  end function drag_if
  end interface

  real(kind=R8), parameter :: toll = 1.e-20_R8
  
contains

  subroutine assign_drag(drag_word)
    implicit none
    character(len=*), intent(in) :: drag_word

    select case (drag_word)
    case ('Newton')
      drag => drag_Newton
    case ('Stokes')
      drag => drag_Stokes
    case ('Schlichting')
      drag => drag_Schlichting
    case ('Schiller-Naumann')
      drag => drag_Schiller_Naumann
    case ('Chang')
      drag => drag_Chang
    case ('Wen-Yu')
      drag => drag_Wen_Yu
    case ('Putnam')
      drag => drag_Putnam
    case ('Clift-Gauvin')
      drag => drag_Clift_Gauvin
    case ('Morsi-Alexander')
      drag => drag_Morsi_Alexander
    case ('Carlson-Hoglund')
      drag => drag_Carlson_Hoglund
    case ('Henderson')
      drag => drag_Henderson
    case ('Crowe')
      drag => drag_Crowe
    case ('Hermsen')
      drag => drag_Hermsen
    case default
      write(*,*)
      write(*,*)
      write(*,*) "Wrong drag input ---> "//drag_word
      write(*,*) "Choose one of the following :"
      write(*,*) "- Newton "
      write(*,*) "- Stokes "
      write(*,*) "- Schlichting "
      write(*,*) "- Schiller-Naumann "
      write(*,*) "- Chang "
      write(*,*) "- Wen-Yu "
      write(*,*) "- Putnam "
      write(*,*) "- Clift-Gauvin "
      write(*,*) "- Morsi-Alexander "
      write(*,*) "- Carlson-Hoglund "
      write(*,*) "- Henderson "
      write(*,*) "- Crowe "
      write(*,*) "- Hermsen "
      write(*,*)
      stop
    end select

  end subroutine assign_drag

  !> Newton Drag Model
  function drag_Newton(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    Cd = 0.45_R8
    
  end function drag_Newton

  !> Stokes Drag Model
  function drag_Stokes(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll)
    
  end function drag_Stokes

  !> Schlichting Drag Model
  function drag_Schlichting(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll) * (1._R8+3._R8*Re/16._R8)
    
  end function drag_Schlichting

  !> Schiller Naumann Drag Model
  function drag_Schiller_Naumann(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll) * (1._R8+0.15_R8*Re**0.687_R8)
    
  end function drag_Schiller_Naumann

  !> Chang Drag Model
  function drag_Chang(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll) * (1._R8+0.15_R8*Re**0.687_R8) + 0.42_R8 / (1._R8+42500_R8*(Re+toll)**(-1.16_R8))
    
  end function drag_Chang

  !> Wen Yu Drag Model
  function drag_Wen_Yu(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd
    
    if (Re <= 1000) then
      Cd = 24._R8/(Re+toll) * (1._R8+0.15_R8*Re**0.687_R8)
    else
      Cd = 0.43_R8
    endif

  end function drag_Wen_Yu

  !> Putnam Drag Model
  function drag_Putnam(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    if (Re < 1000) then
      Cd = 24._R8/(Re+toll) * (1._R8+Re**(2._R8/3._R8)/6._R8)
    else
      Cd = 0.4392_R8
    endif
    
  end function drag_Putnam

  !> Clift Gauvin Drag Model
  function drag_Clift_Gauvin(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll) * (1._R8+0.15_R8*Re**0.687_R8+0.0175_R8*Re/(1._R8+4.25_R8*1e4*Re**(-1.16_R8)))
    
  end function drag_Clift_Gauvin  

  !> Morsi Alexander Drag Model
  function drag_Morsi_Alexander(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re,Ma,G,Tr
    real(R8) :: Cd, a1, a2, a3

    if (Re <= 0.1_R8) then
      a1 = 0_R8; a2 = 24_R8; a3 = 0_R8;
    elseif (Re <= 1.0_R8 .and. Re > 0.1_R8) then
      a1 = 3.69_R8; a2 = 22.73_R8; a3 = 0.0903_R8;
    elseif (Re <= 10_R8 .and. Re > 1.0_R8) then
      a1 = 1.222_R8; a2 = 29.1667_R8; a3 = -3.8889_R8;
    elseif (Re <= 100_R8 .and. Re > 10.0_R8) then
      a1 = 0.6167_R8; a2 = 46.5_R8; a3 = -116.67_R8;
    elseif (Re <= 1000_R8 .and. Re > 100.0_R8) then
      a1 = 0.3644_R8; a2 = 98.33_R8; a3 = -2778_R8;
    elseif (Re <= 5000_R8 .and. Re > 1000_R8) then
      a1 = 0.357_R8; a2 = 148.62_R8; a3 = -4.75e+4_R8;
    elseif (Re <= 10000_R8 .and. Re > 5000_R8) then
      a1 = 0.46_R8; a2 = -490.546_R8; a3 = 57.87e+4_R8;
    elseif (Re >= 10000_R8) then
      a1 = 0.5191_R8; a2 = -1662.5_R8; a3 = 5.4167e+6_R8;
    end if
    Cd = a1+a2/(Re+toll)+a3/(Re*Re+toll)

  end function drag_Morsi_Alexander

  !> Carlson Hoglund Drag Model
  function drag_Carlson_Hoglund(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd, Cd0

    Cd0 = drag_Wen_Yu(Re,Ma,G,Tr)
    Cd  = Cd0 * (1._R8 + exp(-0.427_R8/(Ma**4.63_R8+toll) - 3._R8/(Re**0.88_R8+toll))) / &
                  (1._R8 + Ma/(Re+toll)*(3.82_R8+1.28_R8*exp(-1.25_R8*Re/(Ma+toll))))    
    
  end function drag_Carlson_Hoglund

  !> Henderson Drag Model
  !> Subsonic Flow
  function drag_Henderson_1(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8) :: Cd

    Cd = 24._R8 * (Re + Ma*(0.5_R8*G)**0.5_R8*(4.33_R8 + (3.65_R8-1.53_R8*Tr)/(1._R8+0.353_R8*Tr)*exp(-0.247_R8*Re/(Ma*(0.5_R8*G) &
          **0.5_R8+toll)))+toll)**(-1._R8) + exp(-0.5_R8*Ma/(Re**0.5_R8+toll))*((4.5_R8+0.38_R8*(0.03_R8*Re+0.48_R8*Re**0.5_R8)) /&
          (1._R8+0.03_R8*Re+0.48_R8*Re**0.5_R8) + 0.1_R8*Ma**2._R8 + 0.2_R8*Ma**8._R8) + 0.6_R8*Ma*(0.5_R8*G)**0.5_R8*(1._R8-exp  &
          (-Ma/(Re+toll)))

  end function drag_Henderson_1

  !> High Supersonic Flow
  function drag_Henderson_2(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8) :: Cd

    Cd = (0.9_R8 + 0.34_R8/(Ma**2._R8+toll) + 1.86_R8*(Ma/(Re+toll))**0.5_R8*(2._R8 + 2._R8/(Ma**2*0.5_R8*G+toll) + (1.058_R8/ &
          (Ma*(0.5_R8*G)**0.5_R8+toll))*Tr**0.5_R8) - 1_R8/(Ma**4._R8*G**2._R8*0.25_R8+toll)) / (1._R8 + 1.86_R8*(Ma/(Re+toll))&
          **0.5_R8)

  end function drag_Henderson_2

  function drag_Henderson(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd, Ma1, Ma2, Cd1, Cd2

    if (Ma<=1) then
      Cd = drag_Henderson_1(Re,Ma,G,Tr)
    elseif (Ma>=1.75) then
      Cd = drag_Henderson_2(Re,Ma,G,Tr)
    else
      Ma1 = 1.00
      Ma2 = 1.75
      Cd1 = drag_Henderson_1(Re,Ma1,G,Tr)
      Cd2 = drag_Henderson_2(Re,Ma2,G,Tr)
      Cd = Cd1 + 0.75_R8*(Ma-1._R8)*(Cd2-Cd1)
    endif

  end function drag_Henderson

  function drag_Crowe(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8) :: Cd, Cd0, gfun, hfun

    gfun = 1.25_R8 * (1._R8+dtanh(0.77_R8*dlog10(Re)-1.92_R8))
    gfun = 10**gfun
    hfun = 2.3_R8 + 1.7_R8*Tr**0.5_R8 - 2.3_R8*dtanh(1.17_R8*dlog10(Ma))

    Cd0 = drag_Wen_Yu(Re,Ma,G,Tr)
    Cd = 2._R8 + (Cd0-2._R8)*exp(-3.07_R8*G**0.5_R8*Ma/(Re+toll)*gfun) + hfun/(Ma*G**0.5_R8+toll)*exp(-Re/(2*Ma+toll))

  end function drag_Crowe

  function drag_Hermsen(Re,Ma,G,Tr) result(Cd)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8) :: Cd, Cd0, gfun, hfun

    gfun = (1._R8 + Re*(12.278_R8+0.548_R8*Re)) / (1._R8+11.278_R8*Re)
    hfun = 5.6_R8/(1._R8+Ma) + 1.7_R8*Tr**0.5_R8

    Cd0 = drag_Wen_Yu(Re,Ma,G,Tr)
    Cd = 2._R8 + (Cd0-2._R8)*exp(-3.07_R8*G**0.5_R8*Ma/(Re+toll)*gfun) + hfun/(Ma*G**0.5_R8+toll)*exp(-Re/(2*Ma+toll))
    
  end function drag_Hermsen

end module ICE_Lib_Drag
