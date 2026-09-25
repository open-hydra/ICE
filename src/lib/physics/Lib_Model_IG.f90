!> Isotropic Gaussian closure IG functions
module ICE_Lib_IG
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Global_m
  use ICE_Config_Types_m, only: condensed_phase_t
  use ICE_Lib_Properties, only: mat_cp_const, mat_e, mat_T_from_rhoe
  implicit none

contains

  function prim_2_cons_IG(prim, mat) result(cons)
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

    cons(5) = prim(1)*norm2V + 3._R8*prim(5)
    
    if (mat%cs_varies) then
      cons(6) = prim(1)*mat_e(mat, prim(6)) + 0.5_R8*cons(5)
    else
      cons(6) = prim(1)*mat_cp_const(mat)*prim(6) + 0.5_R8*cons(5)
    endif
    
    cons(7) = prim(7)
  
  end function prim_2_cons_IG


  function cons_2_prim_IG(cons, mat) result(prim)
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

    prim(5) = 1._R8/3._R8 * (cons(5) - prim(1)*norm2V)
    
    if (mat%cs_varies) then
      prim(6) = mat_T_from_rhoe(mat, cons(6) - 0.5_R8*cons(5), prim(1), eps)
    else
      prim(6) = (cons(6) - 0.5_R8*cons(5)) / (mat_cp_const(mat)*prim(1)+eps)
    endif
    
    prim(7) = cons(7)

    if (prim(1)<eps) then
      prim = eps
    endif

    if (prim(5)<eps) then
      prim(5) = eps
    endif

    if (prim(6)<eps) then
      prim(6) = eps
    endif

    if (prim(7)<eps) then
      prim(7) = eps
    endif

  end function cons_2_prim_IG


  function check_prim_IG(prim) result(prim_status)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    logical                   :: prim_status

    !> Check density and pseudo pressure
    if (prim(1) < 0._R8 .or. prim(5) < 0._R8 .or. any(isnan(prim))) then
      prim_status = .false.
    else
      prim_status = .true.
    endif

  end function check_prim_IG


  function pressure_make_IG(prim)  result(pressure)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: pressure

    pressure = (prim(5) + prim(5) + prim(5)) / 3._R8

  end function pressure_make_IG
  

  function sound_make_IG(prim)  result(sound)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: sound

    real(kind=R8), parameter  :: eps = 1e-25

    sound = sqrt( 3._R8*pressure_make_IG(prim)/(prim(1)+eps) )

  end function sound_make_IG


  !> Fastest signal speed across a face of normal n: the isotropic sqrt(3P/rho)
  function wavespeed_make_IG(prim,normal) result(speed)
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    real(kind=R8)             :: speed

    real(kind=R8), parameter  :: eps = 1e-25

    speed = sqrt( 3._R8*pressure_make_IG(prim)/(prim(1)+eps) )

  end function wavespeed_make_IG


  function flux_make_IG(prim,normal,mat) result(flux)
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim))

    real(kind=R8) :: norm2V
    real(kind=R8) :: un
    
    norm2V  = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

    un      = prim(2)*normal(1) + prim(3)*normal(2) + prim(4)*normal(3)

    flux(1) = prim(1)*un

    flux(2) = flux(1)*prim(2) + prim(5)*normal(1)
    flux(3) = flux(1)*prim(3) + prim(5)*normal(2)
    flux(4) = flux(1)*prim(4) + prim(5)*normal(3)

    flux(5) = flux(1)*norm2V + 5._R8*prim(5)*un
   
    if (mat%cs_varies) then
      flux(6) = flux(1)*(0.5_R8*norm2V + mat_e(mat, prim(6))) + 2.5_R8*prim(5)*un
    else
      flux(6) = flux(1)*(0.5_R8*norm2V + mat_cp_const(mat)*prim(6)) + 2.5_R8*prim(5)*un
    endif
   
    flux(7) = prim(7)*un

  end function flux_make_IG


  function source_make_IG(prim,force,mat) result(source)
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
 
    source(5) = - force(1)*(norm2V + 3._R8*prim(5)/prim(1)) + &
                  2._R8*prim(1)*(force(2)*prim(2)+force(3)*prim(3)+force(4)*prim(4))/(force(6)+eps) - &
                  6._R8*prim(5)/(force(6)+eps)

    if (mat%cs_varies) then
      source(6) = - force(1)*mat_e(mat, prim(6)) - force(1)*mat%lv_al + force(5) + 0.5_R8*source(5)
    else
      source(6) = - force(1)*mat_cp_const(mat)*prim(6) - force(1)*mat%lv_al + force(5) + 0.5_R8*source(5)
    endif
    
    source(7) = 0._R8

  end function source_make_IG


end module ICE_Lib_IG
