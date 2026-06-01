module ICE_Lib_Limiters
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: assign_limiter

  !> Concrete "limiter" procedure pointing to one of the function realizations
  procedure(limiter_if), pointer, public :: limiter

  !> Abstract interface relative to the "limiter" procedure
  abstract interface
  pure function limiter_if(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
  end function limiter_if
  end interface
  
contains

  subroutine assign_limiter(limiter_word)
    implicit none
    character(len=*), intent(in) :: limiter_word

    select case (limiter_word)
    case ('IORD')
      limiter => limiter_IORD
    case ('MINMOD')
      limiter => limiter_MINMOD
    case ('VANALBADA')
      limiter => limiter_VANALBADA
    case ('VANLEER')
      limiter => limiter_VANLEER
    case ('OSPRE')
      limiter => limiter_OSPRE
    case ('UMIST')
      limiter => limiter_UMIST
    case ('OSHER')
      limiter => limiter_OSHER
    case ('SWEBY')
      limiter => limiter_SWEBY
    case ('MC')
      limiter => limiter_MC
    case ('KOREN')
      limiter => limiter_KOREN
    case ('SUPERBEE')
      limiter => limiter_SUPERBEE
    case default
     write(*,*)
     write(*,*)
     write(*,*) "Wrong limiter input ---> "//limiter_word
     write(*,*) "Choose one of the following :"
     write(*,*) "- IORD "
     write(*,*) "- MINMOD "
     write(*,*) "- VANALBADA "
     write(*,*) "- VANLEER "
     write(*,*) "- OSPRE "
     write(*,*) "- UMIST "
     write(*,*) "- OSHER "
     write(*,*) "- SWEBY "
     write(*,*) "- MC "
     write(*,*) "- KOREN "
     write(*,*) "- SUPERBEE "
     write(*,*)
     stop
    end select

  end subroutine assign_limiter

  !> First order limiter
  pure function limiter_IORD(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    result = 0._R8

  end function limiter_IORD

  !> Minmod flux limiter
  pure function limiter_MINMOD(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([1d0,ru])])
    result = phi*x

  end function limiter_MINMOD

  !> Van Albada flux limiter
  pure function limiter_VANALBADA(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = (ru*ru+ru)/(ru*ru+1d0)
    result = phi*x

  end function limiter_VANALBADA

  !> Van Leer flux limiter
  pure function limiter_VANLEER(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = (ru+abs(ru))/(1d0+abs(ru))
    result = phi*x

  end function limiter_VANLEER

  !> Ospre flux limiter
  pure function limiter_OSPRE(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = (1.5d0*(ru*ru+ru))/(ru*ru+ru+1d0)
    result = phi*x

  end function limiter_OSPRE

  !> Umist flux limiter
  pure function limiter_UMIST(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([2d0*ru,(0.25d0+0.75d0*ru),(0.75d0+0.25d0*ru),2d0])])
    result = phi*x

  end function limiter_UMIST
  
  !> Osher flux limiter
  pure function limiter_OSHER(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([ru,1.5d0])])
    result = phi*x

  end function limiter_OSHER

  !> Sweby flux limiter
  pure function limiter_SWEBY(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([1.5d0*ru,1d0]),minval([ru,1.5d0])])
    result = phi*x

  end function limiter_SWEBY

  !> Monotonized central flux limiter
  pure function limiter_MC(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([2d0*ru,0.5d0*(1d0+ru),2d0])])
    result = phi*x

  end function limiter_MC

  !> Koren flux limiter
  pure function limiter_KOREN(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([2d0*ru,minval([(1d0+2d0*ru)/3d0,2d0])])])
    result = phi*x

  end function limiter_KOREN

  !> Superbee flux limiter
  pure function limiter_SUPERBEE(x,y) result(result)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: x, y
    real(kind=R8) :: ru
    real(kind=R8) :: phi
    real(kind=R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,min(2d0*ru,1d0),min(ru,2d0)])
    result = phi*x

  end function limiter_SUPERBEE

end module ICE_Lib_Limiters
