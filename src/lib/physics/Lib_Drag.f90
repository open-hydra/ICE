!>@brief Drag-coefficient correlations Cd(Re, Ma, gamma, Tr) selected by keyword.
!> Ported from IGLOO's Lib_Drag. The correlations are pure functions and the choice
!> travels as an integer, so nothing mutable is shared between threads.
!>
!> The Chang correlation is deliberately absent: `24/Re (1 + 0.15 Re^0.687) +
!> 0.42/(1 + 42500 Re^-1.16)` is the Clift-Gauvin correlation written differently
!> (24 * 0.0175 = 0.42), so the two were the same function under two names.
module ICE_Lib_Drag
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private

  public :: assign_drag
  public :: drag

contains

  !> Maps the [ICE-Scheme] drag keyword to the selector used by drag().
  subroutine assign_drag(drag_word, dragSelect)
    implicit none
    character(len=*), intent(in)  :: drag_word
    integer(kind=I4), intent(out) :: dragSelect

    select case (drag_word)
    case ('Newton')
      dragSelect = 1
    case ('Stokes')
      dragSelect = 2
    case ('Schlichting')
      dragSelect = 3
    case ('Schiller-Naumann')
      dragSelect = 4
    case ('Wen-Yu')
      dragSelect = 5
    case ('Putnam')
      dragSelect = 6
    case ('Clift-Gauvin')
      dragSelect = 7
    case ('Morsi-Alexander')
      dragSelect = 8
    case ('Carlson-Hoglund')
      dragSelect = 9
    case ('Henderson')
      dragSelect = 10
    case ('Crowe')
      dragSelect = 11
    case ('Hermsen')
      dragSelect = 12
    case default
      write(*,*)
      if (drag_word == 'none') then
        write(*,'(A)') ' [ERROR] [ICE::assign_drag] drag is not set: [ICE-Physics] drag is required for a coupled run'
      else
        write(*,*) "Wrong drag input ---> "//drag_word
      endif
      if (drag_word == 'Chang') &
        write(*,*) "Chang is the Clift-Gauvin correlation written differently: use Clift-Gauvin."
      write(*,*) "Choose one of the following :"
      write(*,*) "- Newton "
      write(*,*) "- Stokes "
      write(*,*) "- Schlichting "
      write(*,*) "- Schiller-Naumann "
      write(*,*) "- Wen-Yu "
      write(*,*) "- Putnam "
      write(*,*) "- Clift-Gauvin "
      write(*,*) "- Morsi-Alexander "
      write(*,*) "- Carlson-Hoglund "
      write(*,*) "- Henderson "
      write(*,*) "- Crowe "
      write(*,*) "- Hermsen "
      write(*,*)
      if (drag_word == 'none') error stop 1
      error stop 'ICE: unknown drag model'
    end select

  end subroutine assign_drag


  !> Dispatches to the selected drag correlation.
  pure function drag(Re,Ma,G,Tr,dragSelect) result(Cd)
    implicit none
    integer(kind=I4), intent(in) :: dragSelect
    real(kind=R8),    intent(in) :: Re, Ma, G, Tr
    real(kind=R8) :: Cd

    select case (dragSelect)
    case (1) !> Newton Drag Model
      Cd = 0.45_R8
    case (2)
      Cd = drag_Stokes(Re)
    case (3)
      Cd = drag_Schlichting(Re)
    case (4)
      Cd = drag_Schiller_Naumann(Re)
    case (5)
      Cd = drag_Wen_Yu(Re)
    case (6)
      Cd = drag_Putnam(Re)
    case (7)
      Cd = drag_Clift_Gauvin(Re)
    case (8)
      Cd = drag_Morsi_Alexander(Re)
    case (9)
      Cd = drag_Carlson_Hoglund(Re,Ma)
    case (10)
      Cd = drag_Henderson(Re,Ma,G,Tr)
    case (11)
      Cd = drag_Crowe(Re,Ma,G,Tr)
    case (12)
      Cd = drag_Hermsen(Re,Ma,G,Tr)
    end select

  end function drag


  !> Stokes Drag Model
  pure function drag_Stokes(Re) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(R8), parameter :: toll=1e-20
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll)

  end function drag_Stokes


  !> Schlichting Drag Model
  pure function drag_Schlichting(Re) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(R8), parameter :: toll=1e-20
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll) * (1._R8+3._R8*Re/16._R8)

  end function drag_Schlichting


  !> Schiller Naumann Drag Model
  pure function drag_Schiller_Naumann(Re) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(R8), parameter :: toll=1e-20
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll) * (1._R8+0.15_R8*Re**0.687_R8)

  end function drag_Schiller_Naumann


  !> Wen Yu Drag Model
  pure function drag_Wen_Yu(Re) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(R8), parameter :: toll=1e-20
    real(kind=R8) :: Cd

    if (Re <= 1000) then
      Cd = 24._R8/(Re+toll) * (1._R8+0.15_R8*Re**0.687_R8)
    else
      Cd = 0.43_R8   ! plateau per Shimada2006 eq.16
    endif

  end function drag_Wen_Yu


  !> Putnam Drag Model
  pure function drag_Putnam(Re) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(R8), parameter :: toll=1e-20
    real(kind=R8) :: Cd

    if (Re < 1000) then
      Cd = 24._R8/(Re+toll) * (1._R8+Re**(2._R8/3._R8)/6._R8)
    else
      Cd = 0.4392_R8 ! plateau per Shimada2006 eq.17
    endif

  end function drag_Putnam


  !> Clift Gauvin Drag Model
  pure function drag_Clift_Gauvin(Re) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(R8), parameter :: toll=1e-20
    real(kind=R8) :: Cd

    Cd = 24._R8/(Re+toll) * (1._R8+0.15_R8*Re**0.687_R8+0.0175_R8*Re/(1._R8+4.25_R8*1e4*Re**(-1.16_R8)))

  end function drag_Clift_Gauvin


  !> Carlson Hoglund Drag Model
  pure function drag_Carlson_Hoglund(Re,Ma) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re, Ma
    real(R8), parameter :: toll=1e-20
    real(kind=R8) :: Cd, Cd0

    Cd0 = drag_Wen_Yu(Re)
    Cd  = Cd0 * (1._R8 + exp(-0.427_R8/(Ma**4.63_R8+toll) - 3._R8/(Re**0.88_R8+toll))) / &
                  (1._R8 + Ma/(Re+toll)*(3.82_R8+1.28_R8*exp(-1.25_R8*Re/(Ma+toll))))

  end function drag_Carlson_Hoglund


  !> Henderson drag model, subsonic branch (Ma <= 1).
  pure function drag_Henderson_1(Re,Ma,G,Tr) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8), parameter :: toll=1e-20
    real(R8) :: Cd

    Cd = 24._R8 * (Re + Ma*(0.5_R8*G)**0.5_R8*(4.33_R8 + (3.65_R8-1.53_R8*Tr)/(1._R8+0.353_R8*Tr)*exp(-0.247_R8*Re/(Ma*(0.5_R8*G) &
          **0.5_R8+toll)))+toll)**(-1._R8) + exp(-0.5_R8*Ma/(Re**0.5_R8+toll))*((4.5_R8+0.38_R8*(0.03_R8*Re+0.48_R8*Re**0.5_R8)) /&
          (1._R8+0.03_R8*Re+0.48_R8*Re**0.5_R8) + 0.1_R8*Ma**2._R8 + 0.2_R8*Ma**8._R8) + 0.6_R8*Ma*(0.5_R8*G)**0.5_R8*(1._R8-exp  &
          (-Ma/(Re+toll)))

  end function drag_Henderson_1


  !> Henderson drag model, supersonic branch (Ma >= 1.75).
  pure function drag_Henderson_2(Re,Ma,G,Tr) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8), parameter :: toll=1e-20
    real(R8) :: Cd

    Cd = (0.9_R8 + 0.34_R8/(Ma**2._R8+toll) + 1.86_R8*(Ma/(Re+toll))**0.5_R8*(2._R8 + 2._R8/(Ma**2*0.5_R8*G+toll) + (1.058_R8/ &
          (Ma*(0.5_R8*G)**0.5_R8+toll))*Tr**0.5_R8) - 1_R8/(Ma**4._R8*G**2._R8*0.25_R8+toll)) / (1._R8 + 1.86_R8*(Ma/(Re+toll))&
          **0.5_R8)

  end function drag_Henderson_2


  !> Henderson drag model: subsonic and supersonic branches, linear bridge for 1 < Ma < 1.75.
  pure function drag_Henderson(Re,Ma,G,Tr) result(Cd)
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
      Cd = Cd1 + 4._R8/3._R8*(Ma-1._R8)*(Cd2-Cd1)  ! 4/3 = 1/(1.75-1); Shimada eq.22 prints 3/4 (typo)
    endif

  end function drag_Henderson


  !> Crowe drag model (Shimada2006 eq.23).
  pure function drag_Crowe(Re,Ma,G,Tr) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8), parameter :: toll=1e-20
    real(R8) :: Cd, Cd0, gfun, hfun

    gfun = 1.25_R8 * (1._R8+dtanh(0.77_R8*dlog10(Re+toll)-1.92_R8))
    gfun = 10**gfun
    hfun = 2.3_R8 + 1.7_R8*Tr**0.5_R8 - 2.3_R8*dtanh(1.17_R8*dlog10(Ma+toll))

    Cd0 = drag_Wen_Yu(Re)
    Cd = 2._R8 + (Cd0-2._R8)*exp(-3.07_R8*G**0.5_R8*Ma*gfun/(Re+toll)) + hfun/(Ma*G**0.5_R8+toll)*exp(-Re/(2*Ma+toll))

  end function drag_Crowe


  !> Hermsen drag model (Shimada2006 eq.24).
  pure function drag_Hermsen(Re,Ma,G,Tr) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re, Ma, G, Tr
    real(R8), parameter :: toll=1e-20
    real(R8) :: Cd, Cd0, gfun, hfun

    gfun = (1._R8 + Re*(12.278_R8+0.548_R8*Re)) / (1._R8+11.278_R8*Re)
    hfun = 5.6_R8/(1._R8+Ma) + 1.7_R8*Tr**0.5_R8

    Cd0 = drag_Wen_Yu(Re)
    Cd = 2._R8 + (Cd0-2._R8)*exp(-3.07_R8*G**0.5_R8*Ma*gfun/(Re+toll)) + hfun/(Ma*G**0.5_R8+toll)*exp(-Re/(2*Ma+toll))

  end function drag_Hermsen


  !> Morsi-Alexander drag model: a1 + a2/Re + a3/Re^2, piecewise over eight Re ranges.
  pure function drag_Morsi_Alexander(Re) result(Cd)
    implicit none
    real(kind=R8), intent(in) :: Re
    real(R8), parameter :: toll=1e-20
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

end module ICE_Lib_Drag
