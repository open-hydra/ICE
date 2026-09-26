!> Monokinetic closure MK functions
module ICE_Lib_MK
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Global_m
  use ICE_Config_Types_m, only: condensed_phase_t
  use ICE_Lib_Properties, only: mat_cp_const, mat_e, mat_T_from_rhoe
  implicit none

contains

  function prim_2_cons_MK(prim, mat) result(cons)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: cons(size(prim))

    real(kind=R8) :: norm2V

    norm2V  = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

    cons(1) = prim(1)

    cons(2) = prim(1)*prim(2)
    cons(3) = prim(1)*prim(3)
    cons(4) = prim(1)*prim(4)

    if (mat%cs_varies) then
      cons(5) = prim(1)*mat_e(mat, prim(5)) + 0.5_R8*prim(1)*norm2V
    else
      cons(5) = prim(1)*mat_cp_const(mat)*prim(5) + 0.5_R8*prim(1)*norm2V
    endif

    cons(6) = prim(6)
  
  end function prim_2_cons_MK


  function cons_2_prim_MK(cons, mat) result(prim)
    implicit none
    real(kind=R8), intent(in) :: cons(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: prim(size(cons))
    
    real(kind=R8) :: norm2V

    real(kind=R8), parameter  :: eps = 1e-25

    prim(1) = cons(1)

    prim(2) = cons(2) / (prim(1)+eps)
    prim(3) = cons(3) / (prim(1)+eps)
    prim(4) = cons(4) / (prim(1)+eps)

    norm2V  = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

    if (mat%cs_varies) then
      prim(5) = mat_T_from_rhoe(mat, cons(5) - 0.5_R8*prim(1)*norm2V, prim(1), eps)
    else
      prim(5) = (cons(5) - 0.5_R8*prim(1)*norm2V) / (mat_cp_const(mat)*prim(1)+eps)
    endif

    prim(6) = cons(6)

    if (prim(1)<eps) then
      prim = eps
    endif

    if (prim(5)<eps) then
      prim(5) = eps
    endif

    if (prim(6)<eps) then
      prim(6) = eps
    endif

  end function cons_2_prim_MK


  function check_prim_MK(prim) result(prim_status)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    logical                   :: prim_status

    !> Check density
    if (prim(1) < 0._R8 .or. any(isnan(prim))) then
      prim_status = .false.
    else
      prim_status = .true.
    endif

  end function check_prim_MK


  function pressure_make_MK(prim)  result(pressure)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: pressure

    pressure = 0._R8

  end function pressure_make_MK
  

  function sound_make_MK(prim)  result(sound)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: sound

    sound = 0._R8

  end function sound_make_MK


  !> A pressureless cloud has no signal speed of its own
  function wavespeed_make_MK(prim,normal) result(speed)
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    real(kind=R8)             :: speed

    speed = 0._R8

  end function wavespeed_make_MK


  function flux_make_MK(prim,normal,mat) result(flux)
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim))

    real(kind=R8) :: norm2V
    real(kind=R8) :: un
    
    norm2V  = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)
    
    un      = prim(2)*normal(1) + prim(3)*normal(2) + prim(4)*normal(3)

    flux(1) = prim(1)*un

    flux(2) = flux(1)*prim(2)
    flux(3) = flux(1)*prim(3)
    flux(4) = flux(1)*prim(4)

    if (mat%cs_varies) then
      flux(5) = flux(1)*(0.5_R8*norm2V + mat_e(mat, prim(5)))
    else
      flux(5) = flux(1)*(0.5_R8*norm2V + mat_cp_const(mat)*prim(5))
    endif

    flux(6) = prim(6)*un

  end function flux_make_MK


  function source_make_MK(prim,force,mat) result(source)
    implicit none
    real(kind=R8), intent(in) :: prim(:), force(6)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: source(size(prim))

    real(kind=R8) :: norm2V

    real(kind=R8), parameter  :: eps = 1e-25

    norm2V    = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

    source(1) = - force(1)

    source(2) = - force(1)*prim(2) + prim(1)*force(2)/(force(6)+eps)
    source(3) = - force(1)*prim(3) + prim(1)*force(3)/(force(6)+eps)
    source(4) = - force(1)*prim(4) + prim(1)*force(4)/(force(6)+eps)

    if (mat%cs_varies) then
      source(5) = - force(1)*(0.5_R8*norm2V + mat_e(mat, prim(5))) - force(1)*mat%lv_al + force(5) + &
                    prim(1)*(force(2)*prim(2)+force(3)*prim(3)+force(4)*prim(4))/(force(6)+eps)
    else
      source(5) = - force(1)*(0.5_R8*norm2V + mat_cp_const(mat)*prim(5)) - force(1)*mat%lv_al + force(5) + &
                    prim(1)*(force(2)*prim(2)+force(3)*prim(3)+force(4)*prim(4))/(force(6)+eps)
    endif

    source(6) = 0._R8

  end function source_make_MK


end module ICE_Lib_MK
