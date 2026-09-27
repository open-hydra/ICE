!>@brief The closures of a solidifying family: the base closure (MK, IG, AG) on its own slots 1:nb, then the frozen
!> fraction f (nb+1) and the nucleated fraction chi (nb+2), both carried with the mass. The energy holds the latent part
!> -f*L(T); after every conversion the phase follows from the energy and chi (ICE_Lib_Solidification::solid_state).
module ICE_Lib_Solid
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Config_Types_m,     only: condensed_phase_t
  use ICE_Lib_Properties,     only: mat_cp_const
  use ICE_Lib_Solidification, only: solid_de, solid_state
  use ICE_Lib_MK, only: prim_2_cons_MK, cons_2_prim_MK, flux_make_MK, source_make_MK
  use ICE_Lib_IG, only: prim_2_cons_IG, cons_2_prim_IG, flux_make_IG, source_make_IG
  use ICE_Lib_AG, only: prim_2_cons_AG, cons_2_prim_AG, flux_make_AG, source_make_AG
  implicit none
  private
  public :: prim_2_cons_MK_S, cons_2_prim_MK_S, flux_make_MK_S, source_make_MK_S
  public :: prim_2_cons_IG_S, cons_2_prim_IG_S, flux_make_IG_S, source_make_IG_S
  public :: prim_2_cons_AG_S, cons_2_prim_AG_S, flux_make_AG_S, source_make_AG_S

  !> The bases' own guard, so T_liq is their temperature bit for bit
  real(kind=R8), parameter :: eps = 1e-25

contains

  function prim_2_cons_MK_S(prim, mat) result(cons)
    real(kind=R8), intent(in) :: prim(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: cons(size(prim))
    cons(1:6) = prim_2_cons_MK(prim(1:6), mat)
    call latent_cons(cons, prim, 6, mat)
  end function prim_2_cons_MK_S

  function cons_2_prim_MK_S(cons, mat) result(prim)
    real(kind=R8), intent(in) :: cons(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: prim(size(cons))
    real(kind=R8) :: norm2V, Tliq
    prim(1:6) = cons_2_prim_MK(cons(1:6), mat)
    norm2V = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)
    Tliq   = (cons(5) - 0.5_R8*prim(1)*norm2V) / (mat_cp_const(mat)*prim(1)+eps)
    call phase_prims(prim, cons, 6, Tliq, mat)
  end function cons_2_prim_MK_S

  function flux_make_MK_S(prim, normal, mat) result(flux)
    real(kind=R8), intent(in) :: prim(:), normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim))
    flux(1:6) = flux_make_MK(prim(1:6), normal, mat)
    call latent_flux(flux, prim, 6, mat)
  end function flux_make_MK_S

  function source_make_MK_S(prim, force, mat) result(source)
    real(kind=R8), intent(in) :: prim(:), force(6)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: source(size(prim))
    source(1:6) = source_make_MK(prim(1:6), force, mat)
    call latent_source(source, prim, force, 6, mat)
  end function source_make_MK_S


  function prim_2_cons_IG_S(prim, mat) result(cons)
    real(kind=R8), intent(in) :: prim(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: cons(size(prim))
    cons(1:7) = prim_2_cons_IG(prim(1:7), mat)
    call latent_cons(cons, prim, 7, mat)
  end function prim_2_cons_IG_S

  function cons_2_prim_IG_S(cons, mat) result(prim)
    real(kind=R8), intent(in) :: cons(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: prim(size(cons))
    real(kind=R8) :: Tliq
    prim(1:7) = cons_2_prim_IG(cons(1:7), mat)
    Tliq = (cons(6) - 0.5_R8*cons(5)) / (mat_cp_const(mat)*prim(1)+eps)
    call phase_prims(prim, cons, 7, Tliq, mat)
  end function cons_2_prim_IG_S

  function flux_make_IG_S(prim, normal, mat) result(flux)
    real(kind=R8), intent(in) :: prim(:), normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim))
    flux(1:7) = flux_make_IG(prim(1:7), normal, mat)
    call latent_flux(flux, prim, 7, mat)
  end function flux_make_IG_S

  function source_make_IG_S(prim, force, mat) result(source)
    real(kind=R8), intent(in) :: prim(:), force(6)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: source(size(prim))
    source(1:7) = source_make_IG(prim(1:7), force, mat)
    call latent_source(source, prim, force, 7, mat)
  end function source_make_IG_S


  function prim_2_cons_AG_S(prim, mat) result(cons)
    real(kind=R8), intent(in) :: prim(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: cons(size(prim))
    cons(1:12) = prim_2_cons_AG(prim(1:12), mat)
    call latent_cons(cons, prim, 12, mat)
  end function prim_2_cons_AG_S

  function cons_2_prim_AG_S(cons, mat) result(prim)
    real(kind=R8), intent(in) :: cons(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: prim(size(cons))
    real(kind=R8) :: Tliq
    prim(1:12) = cons_2_prim_AG(cons(1:12), mat)
    Tliq = (cons(11) - 0.5_R8*(cons(5)+cons(8)+cons(10))) / (mat_cp_const(mat)*prim(1)+eps)
    call phase_prims(prim, cons, 12, Tliq, mat)
  end function cons_2_prim_AG_S

  function flux_make_AG_S(prim, normal, mat) result(flux)
    real(kind=R8), intent(in) :: prim(:), normal(3)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: flux(size(prim))
    flux(1:12) = flux_make_AG(prim(1:12), normal, mat)
    call latent_flux(flux, prim, 12, mat)
  end function flux_make_AG_S

  function source_make_AG_S(prim, force, mat) result(source)
    real(kind=R8), intent(in) :: prim(:), force(6)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)             :: source(size(prim))
    source(1:12) = source_make_AG(prim(1:12), force, mat)
    call latent_source(source, prim, force, 12, mat)
  end function source_make_AG_S


  !> rho*f and rho*chi; the energy slot nb-1 gains rho times the latent part (nothing when f = 0).
  pure subroutine latent_cons(cons, prim, nb, mat)
    real(kind=R8), intent(inout) :: cons(:)
    real(kind=R8), intent(in)    :: prim(:)
    integer,       intent(in)    :: nb
    type(condensed_phase_t), intent(in) :: mat
    cons(nb+1) = prim(1)*prim(nb+1)
    cons(nb+2) = prim(1)*prim(nb+2)
    if (prim(nb+1) /= 0._R8) cons(nb-1) = cons(nb-1) + prim(1)*solid_de(prim(nb-1), prim(nb+1), &
                                          mat_cp_const(mat), mat%cpSol, mat%hFus, mat%Tmelt)
  end subroutine latent_cons

  !> f and chi move with the mass flux; the energy flux gains the latent part (nothing when f = 0).
  pure subroutine latent_flux(flux, prim, nb, mat)
    real(kind=R8), intent(inout) :: flux(:)
    real(kind=R8), intent(in)    :: prim(:)
    integer,       intent(in)    :: nb
    type(condensed_phase_t), intent(in) :: mat
    flux(nb+1) = flux(1)*prim(nb+1)
    flux(nb+2) = flux(1)*prim(nb+2)
    if (prim(nb+1) /= 0._R8) flux(nb-1) = flux(nb-1) + flux(1)*solid_de(prim(nb-1), prim(nb+1), &
                                          mat_cp_const(mat), mat%cpSol, mat%hFus, mat%Tmelt)
  end subroutine latent_flux

  !> Evaporated mass leaves with its f, chi and latent energy (zero here: evaporation is refused with solidification).
  pure subroutine latent_source(source, prim, force, nb, mat)
    real(kind=R8), intent(inout) :: source(:)
    real(kind=R8), intent(in)    :: prim(:), force(6)
    integer,       intent(in)    :: nb
    type(condensed_phase_t), intent(in) :: mat
    source(nb+1) = - force(1)*prim(nb+1)
    source(nb+2) = - force(1)*prim(nb+2)
    if (prim(nb+1) /= 0._R8) source(nb-1) = source(nb-1) - force(1)*solid_de(prim(nb-1), prim(nb+1), &
                                            mat_cp_const(mat), mat%cpSol, mat%hFus, mat%Tmelt)
  end subroutine latent_source

  !> T, f and chi of the converted cell from T_liq = e/c_l (unfloored) and the transported chi. A liquid content keeps
  !  the base's primitives; the base's vacuum branch carries no phase.
  pure subroutine phase_prims(prim, cons, nb, Tliq, mat)
    real(kind=R8), intent(inout) :: prim(:)
    real(kind=R8), intent(in)    :: cons(:), Tliq
    integer,       intent(in)    :: nb
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8) :: chi, T, f, chiOut
    prim(nb+1:nb+2) = 0._R8
    if (cons(1) < eps) return
    chi = min(max(cons(nb+2)/(prim(1)+eps), 0._R8), 1._R8)
    call solid_state(Tliq, chi, mat_cp_const(mat), mat%cpSol, mat%Tmelt, mat%Tnuc, mat%hFus, T, f, chiOut)
    prim(nb+2) = chiOut
    if (f > 0._R8) then
      prim(nb-1) = max(T, eps)
      prim(nb+1) = f
    endif
  end subroutine phase_prims

end module ICE_Lib_Solid
