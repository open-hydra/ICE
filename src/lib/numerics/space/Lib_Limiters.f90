module ICE_Lib_Limiters
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: assign_limiter

  !> The limiter assigned at set-up (its word; 'IORD' for first order) and its forward-Euler step bound (step_bound)
  character(len=16), public, protected :: limiter_name  = 'IORD'
  real(R8),          public, protected :: limiter_bound = 1._R8

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

  !> `IORD` returns a zero slope and is how the solver expresses a first-order
  !  reconstruction. It is not offered as an input: `space-reconstruction` selects it.
  subroutine assign_limiter(limiter_word)
    implicit none
    character(len=*), intent(in) :: limiter_word

    select case (limiter_word)
    case ('IORD')
      limiter => limiter_IORD
    case ('minmod')
      limiter => limiter_MINMOD
    case ('vanalbada')
      limiter => limiter_VANALBADA
    case ('vanleer')
      limiter => limiter_VANLEER
    case ('ospre')
      limiter => limiter_OSPRE
    case ('umist')
      limiter => limiter_UMIST
    case ('osher')
      limiter => limiter_OSHER
    case ('sweby')
      limiter => limiter_SWEBY
    case ('mc')
      limiter => limiter_MC
    case ('koren')
      limiter => limiter_KOREN
    case ('superbee')
      limiter => limiter_SUPERBEE
    case default
     write(*,*)
     write(*,*)
     write(*,*) "Wrong limiter input ---> "//limiter_word
     write(*,*) "Choose one of the following :"
     write(*,*) "- minmod "
     write(*,*) "- vanalbada "
     write(*,*) "- vanleer "
     write(*,*) "- ospre "
     write(*,*) "- umist "
     write(*,*) "- osher "
     write(*,*) "- sweby "
     write(*,*) "- mc "
     write(*,*) "- koren "
     write(*,*) "- superbee "
     write(*,*)
     write(*,*) "For first order, set space-reconstruction = first-order."
     write(*,*)
     error stop 'ICE: unknown flux limiter'
    end select

    limiter_name  = limiter_word
    limiter_bound = step_bound(limiter_word)

  end subroutine assign_limiter

  !> The largest CFL number at which forward Euler keeps the MUSCL reconstruction TVD for 1-D scalar advection:
  !> Harten's sufficient condition in incremental form, over every slope ratio r (backward over forward slope, as
  !> Lib_Reconstruction passes them) and both signs of the speed, nu = 1/(1 + A/2 - Pmin/2) for a > 0 and
  !> 1/(1 + Pmax/2 - B/2) for a < 0, with A, B the largest and smallest phi(r)/r and Pmax, Pmin those of phi(r) of
  !> the formulas below (plan-bucket/records/pre-existing-defects/d7_limiter_bounds.py). Literal constants, so that
  !> cfl = 0.5 with van Leer is not above its bound by a rounding. First order: 1 (not checked).
  pure function step_bound(limiter_word) result(nu)
    implicit none
    character(len=*), intent(in) :: limiter_word
    real(R8) :: nu

    select case (limiter_word)
    case ('vanleer', 'mc', 'superbee', 'koren', 'umist', 'ospre')
      nu = 0.5_R8
    case ('vanalbada')
      nu = 0.585786437626905_R8       ! 2 - sqrt(2), correctly rounded
    case ('sweby', 'osher')
      nu = 0.5714285714285714_R8      ! 4/7
    case ('minmod')
      nu = 0.6666666666666666_R8      ! 2/3
    case default
      nu = 1._R8                      ! IORD
    end select

  end function step_bound

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
