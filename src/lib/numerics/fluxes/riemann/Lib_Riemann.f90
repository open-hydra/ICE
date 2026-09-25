module ICE_Lib_Riemann
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Lib_Model
  use ICE_Config_Types_m, only: condensed_phase_t
  use ICE_Lib_Properties, only: mat_cp_const, mat_e
  implicit none
  private
  public :: assign_riemann

  !> Concrete "drag" procedure pointing to one of the function realizations
  procedure(riemann_if), pointer, public :: riemann

  !> Abstract interface relative to the "limiter" procedure
  abstract interface
  function riemann_if(prim_1,prim_4,normal,mat) result(flux)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    import :: condensed_phase_t
    implicit none
    real(kind=R8), intent(in) :: prim_1(:)
    real(kind=R8), intent(in) :: prim_4(:)
    real(kind=R8), intent(in) :: normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim_1))
  end function riemann_if
  end interface

contains

subroutine assign_riemann(riemann_word)
  implicit none
  character(len=*), intent(in) :: riemann_word

  select case (riemann_word)
  case ('Saurel')
    riemann => riemann_Saurel
  case ('Rusanov')
    riemann => riemann_Rusanov
  case ('HLLE')
    riemann => riemann_HLLE
  case default
      write(*,*)
      write(*,*)
      write(*,*) "Wrong riemann input ---> "//riemann_word
      write(*,*) "Choose one of the following :"
      write(*,*) "- Saurel "
      write(*,*) "- Rusanov "
      write(*,*) "- HLLE "
      write(*,*)
      error stop 'ICE: unknown Riemann solver'
  end select

end subroutine assign_riemann


  !> Saurel Riemann Solver - EX CPM
  function riemann_Saurel(prim_1,prim_4,normal,mat) result(flux)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: prim_1(:)
    real(kind=R8), intent(in) :: prim_4(:)
    real(kind=R8), intent(in) :: normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim_1))

    real(kind=R8) :: ut_1, vt_1, wt_1, uu_1, vv_1, ww_1
    real(kind=R8) :: ut_4, vt_4, wt_4, uu_4, vv_4, ww_4
    real(kind=R8) :: veln_1, veln_4, veln
    real(kind=R8) :: flux_1(size(prim_1)), flux_4(size(prim_4))

    veln_1 = prim_1(2)*normal(1)+prim_1(3)*normal(2)+prim_1(4)*normal(3)
    veln_4 = prim_4(2)*normal(1)+prim_4(3)*normal(2)+prim_4(4)*normal(3)
    veln = 0.5_R8*(veln_1+veln_4)

    if (veln > 0._R8) then
      ut_1 = (prim_1(2)-veln_1*normal(1)); vt_1 = (prim_1(3)-veln_1*normal(2)); wt_1 = (prim_1(4)-veln_1*normal(3))
      uu_1 = veln_1*normal(1)+ut_1; vv_1 = veln_1*normal(2)+vt_1; ww_1 = veln_1*normal(3)+wt_1
      flux(1) = prim_1(1)*veln
      flux(2) = uu_1*flux(1)
      flux(3) = vv_1*flux(1)
      flux(4) = ww_1*flux(1)
      if (mat%cs_varies) then
        flux(5) = flux(1)*(0.5_R8*(uu_1*uu_1+vv_1*vv_1+ww_1*ww_1)+mat_e(mat, prim_1(5)))
      else
        flux(5) = flux(1)*(0.5_R8*(uu_1*uu_1+vv_1*vv_1+ww_1*ww_1)+mat_cp_const(mat)*prim_1(5))
      endif
      flux(6) = prim_1(6)*veln
    elseif (veln < 0._R8) then
      ut_4 = (prim_4(2)-veln_4*normal(1)); vt_4 = (prim_4(3)-veln_4*normal(2)); wt_4 = (prim_4(4)-veln_4*normal(3))
      uu_4 = veln_4*normal(1)+ut_4; vv_4 = veln_4*normal(2)+vt_4; ww_4 = veln_4*normal(3)+wt_4
      flux(1) = prim_4(1)*veln
      flux(2) = uu_4*flux(1)
      flux(3) = vv_4*flux(1)
      flux(4) = ww_4*flux(1)
      if (mat%cs_varies) then
        flux(5) = flux(1)*(0.5_R8*(uu_4*uu_4+vv_4*vv_4+ww_4*ww_4)+mat_e(mat, prim_4(5)))
      else
        flux(5) = flux(1)*(0.5_R8*(uu_4*uu_4+vv_4*vv_4+ww_4*ww_4)+mat_cp_const(mat)*prim_4(5))
      endif
      flux(6) = prim_4(6)*veln
    else 
      flux = 0._R8
    endif

  end function riemann_Saurel


  !> Rusanov Riemann Solver
  function riemann_Rusanov(prim_1,prim_4,normal,mat) result(flux)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: prim_1(:)
    real(kind=R8), intent(in) :: prim_4(:)
    real(kind=R8), intent(in) :: normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim_1))

    real(kind=R8) :: veln_1, veln_4
    real(kind=R8) :: sound_1, sound_4

    real(kind=R8) :: A
    
    real(kind=R8) :: cons_1(size(prim_1)), cons_4(size(prim_4))
    real(kind=R8) :: flux_1(size(prim_1)), flux_4(size(prim_4))

    veln_1 = prim_1(2)*normal(1)+prim_1(3)*normal(2)+prim_1(4)*normal(3)
    veln_4 = prim_4(2)*normal(1)+prim_4(3)*normal(2)+prim_4(4)*normal(3)

    sound_1 = wavespeed_make(prim_1,normal)
    sound_4 = wavespeed_make(prim_4,normal)
    
    cons_1 = prim_2_cons(prim_1, mat)
    cons_4 = prim_2_cons(prim_4, mat)
    
    flux_1 = flux_make(prim_1,normal,mat)
    flux_4 = flux_make(prim_4,normal,mat)

    A = MAX(ABS(veln_1-sound_1),ABS(veln_4-sound_4),ABS(veln_1+sound_1),ABS(veln_4+sound_4))

    flux = 0.5 * (flux_1 + flux_4 - A*(cons_4-cons_1))
     
  end function riemann_Rusanov


  !> HLLE Riemann Solver
  function riemann_HLLE(prim_1,prim_4,normal,mat) result(flux)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: prim_1(:)
    real(kind=R8), intent(in) :: prim_4(:)
    real(kind=R8), intent(in) :: normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim_1))

    real(kind=R8) :: veln_1, veln_4    
    real(kind=R8) :: sound_1, sound_4

    real(kind=R8) :: u_ROE, v_ROE, w_ROE, veln_ROE, sound_ROE
    real(kind=R8) :: S1, S4

    real(kind=R8) :: cons_1(size(prim_1)), cons_4(size(prim_1))
    real(kind=R8) :: flux_1(size(prim_1)), flux_4(size(prim_1))

    veln_1 = prim_1(2)*normal(1)+prim_1(3)*normal(2)+prim_1(4)*normal(3)
    veln_4 = prim_4(2)*normal(1)+prim_4(3)*normal(2)+prim_4(4)*normal(3)

    sound_1 = wavespeed_make(prim_1,normal)
    sound_4 = wavespeed_make(prim_4,normal)
   
    u_ROE = (sqrt(prim_1(1))*prim_1(2)+sqrt(prim_4(1))*prim_4(2)) / (sqrt(prim_1(1))+sqrt(prim_4(1)))
    v_ROE = (sqrt(prim_1(1))*prim_1(3)+sqrt(prim_4(1))*prim_4(3)) / (sqrt(prim_1(1))+sqrt(prim_4(1)))
    w_ROE = (sqrt(prim_1(1))*prim_1(4)+sqrt(prim_4(1))*prim_4(4)) / (sqrt(prim_1(1))+sqrt(prim_4(1)))

    veln_ROE = u_ROE*normal(1)+v_ROE*normal(2)+w_ROE*normal(3)

    sound_ROE = (sqrt(prim_1(1))*sound_1+sqrt(prim_4(1))*sound_4) / (sqrt(prim_1(1))+sqrt(prim_4(1)))

    S1 = min(0.d0,veln_1-sound_1,veln_ROE-sound_ROE)
    S4 = max(0.d0,veln_4+sound_4,veln_ROE+sound_ROE)

    ! S1 = min(0.d0,veln_1-sound_1,veln_4-sound_4)
    ! S4 = max(0.d0,veln_1+sound_1,veln_4+sound_4)

    select case(minloc([-S1,S1*S4,S4],dim=1))
    case(1)
      flux = flux_make(prim_1,normal,mat)
    case(2)
      flux_1 = flux_make(prim_1,normal,mat)
      flux_4 = flux_make(prim_4,normal,mat)
      cons_1 = prim_2_cons(prim_1, mat)
      cons_4 = prim_2_cons(prim_4, mat)
      flux   = (S4*flux_1 - S1*flux_4 + S1*S4*(cons_4-cons_1)) / (S4-S1)
    case(3)
      flux = flux_make(prim_4,normal,mat)
    endselect

  end function riemann_HLLE


end module ICE_Lib_Riemann
