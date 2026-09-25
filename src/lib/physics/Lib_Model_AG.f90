!> Anisotropic Gaussian closure AG functions   
module ICE_Lib_AG
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Global_m
  use ICE_Config_Types_m, only: obj_condensed
  use ICE_Load_Table,     only: get_cs_al
  implicit none

contains

  function prim_2_cons_AG(prim) result(cons)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: cons(size(prim))

    real(kind=R8) :: norm2V

    norm2V   = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

    cons(1)  = prim(1)

    cons(2)  = prim(1)*prim(2)
    cons(3)  = prim(1)*prim(3)
    cons(4)  = prim(1)*prim(4)

    cons(5)  = prim(1)*prim(2)*prim(2) + prim(5)
    cons(6)  = prim(1)*prim(2)*prim(3) + prim(6)
    cons(7)  = prim(1)*prim(2)*prim(4) + prim(7)
    cons(8)  = prim(1)*prim(3)*prim(3) + prim(8)
    cons(9)  = prim(1)*prim(3)*prim(4) + prim(9)
    cons(10) = prim(1)*prim(4)*prim(4) + prim(10)

    cons(11) = prim(1)*get_cs_al(prim(11))*prim(11) + 0.5_R8*(cons(5)+cons(8)+cons(10))

    cons(12) = prim(12)
  
  end function prim_2_cons_AG


  function cons_2_prim_AG(cons) result(prim)
    implicit none
    real(kind=R8), intent(in) :: cons(:)
    real(kind=R8)             :: prim(size(cons))

    real(kind=R8) :: norm2V

    real(kind=R8), parameter  :: eps = 1e-25

    prim(1) = cons(1)

    prim(2) = cons(2) / (prim(1)+eps)
    prim(3) = cons(3) / (prim(1)+eps)
    prim(4) = cons(4) / (prim(1)+eps)

    norm2V  = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

    prim(5)  = cons(5)  - prim(1)*prim(2)*prim(2)
    prim(6)  = cons(6)  - prim(1)*prim(2)*prim(3)
    prim(7)  = cons(7)  - prim(1)*prim(2)*prim(4)
    prim(8)  = cons(8)  - prim(1)*prim(3)*prim(3)
    prim(9)  = cons(9)  - prim(1)*prim(3)*prim(4)
    prim(10) = cons(10) - prim(1)*prim(4)*prim(4)

    prim(11) = (cons(11) - 0.5_R8*(cons(5)+cons(8)+cons(10))) / (obj_condensed%cs_al*prim(1)+eps)   ! initial estimate
    prim(11) = (cons(11) - 0.5_R8*(cons(5)+cons(8)+cons(10))) / (get_cs_al(prim(11))*prim(1)+eps)  ! table correction

    prim(12) = cons(12)

    if (prim(1)<eps) then
      prim = eps
      prim(6) = 0._R8
      prim(7) = 0._R8
      prim(9) = 0._R8
    endif

    if (prim(5)<eps) then
      prim(5) = eps
    endif

    if (abs(prim(6))<eps) then
      prim(6) = 0._R8
    endif

    if (abs(prim(7))<eps) then
      prim(7) = 0._R8
    endif

    if (prim(8)<eps) then
      prim(8) = eps
    endif

    if (abs(prim(9))<eps) then
      prim(9) = 0._R8
    endif

    if (prim(10)<eps) then
      prim(10) = eps
    endif

    if (prim(11)<eps) then
      prim(11) = eps
    endif

    if (prim(12)<eps) then
      prim(12) = eps
    endif

  end function cons_2_prim_AG


  function check_prim_AG(prim) result(prim_status)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    logical                   :: prim_status

    !> Check density and pseudo pressure
    if (prim(1) < 0._R8 .or. prim(5) < 0._R8 .or. prim(8) < 0._R8 .or. prim(10) < 0._R8 .or. any(isnan(prim))) then
      prim_status = .false.
    else
      prim_status = .true.
    endif

  end function check_prim_AG


  function pressure_make_AG(prim) result(pressure)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: pressure

    pressure = (prim(5) + prim(8) + prim(10)) / 3._R8

  end function pressure_make_AG
  

  function sound_make_AG(prim) result(sound)
    implicit none
    real(kind=R8), intent(in) :: prim(:)
    real(kind=R8)             :: sound

    real(kind=R8), parameter  :: eps = 1e-25

    sound = sqrt( 3._R8*pressure_make_AG(prim)/(prim(1)+eps) )

  end function sound_make_AG


  !> Fastest signal speed across a face of normal n: sqrt(3 P_nn/rho), P_nn = n.P.n
  function wavespeed_make_AG(prim,normal) result(speed)
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    real(kind=R8)             :: speed
    real(kind=R8)             :: pnn
    real(kind=R8), parameter  :: eps = 1e-25

    pnn = normal(1)*(prim(5)*normal(1) + prim(6)*normal(2) + prim(7)*normal(3)) + &
          normal(2)*(prim(6)*normal(1) + prim(8)*normal(2) + prim(9)*normal(3)) + &
          normal(3)*(prim(7)*normal(1) + prim(9)*normal(2) + prim(10)*normal(3))
    speed = sqrt( 3._R8*max(pnn,0._R8)/(prim(1)+eps) )

  end function wavespeed_make_AG


  !> Reflect the symmetric tensor (xx,xy,xz,yy,yz,zz) across the plane of unit normal n: H P H, H = I - 2 n n^T
  pure function mirror_tensor_AG(p, n) result(pm)
    implicit none
    real(kind=R8), intent(in) :: p(6), n(3)
    real(kind=R8)             :: pm(6)
    real(kind=R8)             :: t(3), s

    t(1) = p(1)*n(1) + p(2)*n(2) + p(3)*n(3)
    t(2) = p(2)*n(1) + p(4)*n(2) + p(5)*n(3)
    t(3) = p(3)*n(1) + p(5)*n(2) + p(6)*n(3)
    s    = n(1)*t(1) + n(2)*t(2) + n(3)*t(3)

    pm(1) = p(1) - 4._R8*n(1)*t(1) + 4._R8*s*n(1)*n(1)
    pm(2) = p(2) - 2._R8*(n(1)*t(2) + t(1)*n(2)) + 4._R8*s*n(1)*n(2)
    pm(3) = p(3) - 2._R8*(n(1)*t(3) + t(1)*n(3)) + 4._R8*s*n(1)*n(3)
    pm(4) = p(4) - 4._R8*n(2)*t(2) + 4._R8*s*n(2)*n(2)
    pm(5) = p(5) - 2._R8*(n(2)*t(3) + t(2)*n(3)) + 4._R8*s*n(2)*n(3)
    pm(6) = p(6) - 4._R8*n(3)*t(3) + 4._R8*s*n(3)*n(3)

  end function mirror_tensor_AG


  function flux_make_AG(prim,normal) result(flux)
    implicit none
    real(kind=R8), intent(in) :: prim(:), normal(3)
    real(kind=R8)             :: flux(size(prim))

    real(kind=R8) :: norm2V
    real(kind=R8) :: un
    real(kind=R8) :: p1n, p2n, p3n

    real(kind=R8), parameter  :: eps = 1e-25

    norm2V   = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

    un       = prim(2)*normal(1) + prim(3)*normal(2) + prim(4)*normal(3)

    p1n      = prim(5)*normal(1) + prim(6)*normal(2) + prim(7)*normal(3)
    p2n      = prim(6)*normal(1) + prim(8)*normal(2) + prim(9)*normal(3)
    p3n      = prim(7)*normal(1) + prim(9)*normal(2) + prim(10)*normal(3)

    flux(1)  = prim(1)*un

    flux(2)  = flux(1)*prim(2) + p1n
    flux(3)  = flux(1)*prim(3) + p2n
    flux(4)  = flux(1)*prim(4) + p3n

    flux(5)  = flux(1)*prim(2)*prim(2) + un*prim(5) + prim(2)*p1n + prim(2)*p1n
    flux(6)  = flux(1)*prim(2)*prim(3) + un*prim(6) + prim(2)*p2n + prim(3)*p1n
    flux(7)  = flux(1)*prim(2)*prim(4) + un*prim(7) + prim(2)*p3n + prim(4)*p1n
    flux(8)  = flux(1)*prim(3)*prim(3) + un*prim(8) + prim(3)*p2n + prim(3)*p2n
    flux(9)  = flux(1)*prim(3)*prim(4) + un*prim(9) + prim(3)*p3n + prim(4)*p2n
    flux(10) = flux(1)*prim(4)*prim(4) + un*prim(10) + prim(4)*p3n + prim(4)*p3n

    flux(11) = flux(1)*(0.5_R8*norm2V + 0.5_R8*(prim(5)+prim(8)+prim(10))/(prim(1)+eps) + get_cs_al(prim(11))*prim(11)) &
             + prim(2)*p1n + prim(3)*p2n + prim(4)*p3n

    flux(12) = prim(12)*un

  end function flux_make_AG


  function source_make_AG(prim,force) result(source)
    implicit none
    real(kind=R8), intent(in) :: prim(:), force(6)
    real(kind=R8)             :: source(size(prim))

    real(kind=R8), parameter  :: eps = 1e-25

    source(1) = - force(1)

    source(2) = - force(1)*prim(2) + prim(1)*force(2)/(force(6)+eps)
    source(3) = - force(1)*prim(3) + prim(1)*force(3)/(force(6)+eps)
    source(4) = - force(1)*prim(4) + prim(1)*force(4)/(force(6)+eps)

    source(5) = - force(1)*(prim(2)*prim(2)+prim(5)/prim(1)) + &
                  prim(1)*(force(2)*prim(2)+force(2)*prim(2))/(force(6)+eps) - &
                  2._R8*prim(5)/(force(6)+eps)

    source(6) = - force(1)*(prim(2)*prim(3)+prim(6)/prim(1)) + &
                  prim(1)*(force(3)*prim(2)+force(2)*prim(3))/(force(6)+eps) - &
                  2._R8*prim(6)/(force(6)+eps)

    source(7) = - force(1)*(prim(2)*prim(4)+prim(7)/prim(1)) + &
                  prim(1)*(force(4)*prim(2)+force(2)*prim(4))/(force(6)+eps) - &
                  2._R8*prim(7)/(force(6)+eps)

    source(8) = - force(1)*(prim(3)*prim(3)+prim(8)/prim(1)) + &
                  prim(1)*(force(3)*prim(3)+force(3)*prim(3))/(force(6)+eps) - &
                  2._R8*prim(8)/(force(6)+eps)
              
    source(9) = - force(1)*(prim(3)*prim(4)+prim(9)/prim(1)) + &
                  prim(1)*(force(4)*prim(3)+force(3)*prim(4))/(force(6)+eps) - &
                  2._R8*prim(9)/(force(6)+eps)

    source(10) = - force(1)*(prim(4)*prim(4)+prim(10)/prim(1)) + &
                   prim(1)*(force(4)*prim(4)+force(4)*prim(4))/(force(6)+eps) - &
                   2._R8*prim(10)/(force(6)+eps)

    source(11) = - force(1)*get_cs_al(prim(11))*prim(11) - force(1)*obj_condensed%lv_al + force(5) + 0.5_R8*source(5) + 0.5_R8*source(8) + 0.5_R8*source(10)

    source(12) = 0._R8

  end function source_make_AG


end module ICE_Lib_AG
