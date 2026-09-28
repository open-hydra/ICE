!>@brief Solidification with supercooling, recalescence and melting. solidPhaseAtInjection, nucleationJump and hSolid
!> are IGLOO's model-6 functions (same names, same algebra); solid_de and solid_state give the latent part of the energy
!> and the state of a cell from its energy and its nucleated fraction.
module ICE_Lib_Solidification
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: solidPhaseAtInjection, nucleationJump, hSolid, solid_de, solid_state

  !> Phases: liquid, undercooled liquid, freezing plateau at T-melt, solid.
  integer, parameter, public :: phLiquid=0, phUndercooled=1, phPlateau=2, phSolid=3

contains

  !> Phase and frozen fraction at injection: solid at or below T-nuc, undercooled below T-melt, else liquid.
  pure subroutine solidPhaseAtInjection(Tp, Tmelt, Tnuc, phase, f)
    real(R8), intent(in)  :: Tp, Tmelt, Tnuc
    integer,  intent(out) :: phase
    real(R8), intent(out) :: f

    if (Tp <= Tnuc) then
      phase = phSolid;       f = 1._R8
    elseif (Tp < Tmelt) then
      phase = phUndercooled; f = 0._R8
    else
      phase = phLiquid;      f = 0._R8
    endif

  end subroutine solidPhaseAtInjection


  !> Adiabatic recalescence from Tev: back to T-melt with f0 = cpl*(Tmelt-Tev)/hFus when f0 < 1,
  !  else frozen whole at the temperature that conserves hSolid.
  pure subroutine nucleationJump(Tev, cpl, cpSol, Tmelt, hFus, Tafter, fAfter, phaseAfter)
    real(R8), intent(in)  :: Tev, cpl, cpSol, Tmelt, hFus
    real(R8), intent(out) :: Tafter, fAfter
    integer,  intent(out) :: phaseAfter
    real(R8) :: f0

    f0 = cpl*(Tmelt - Tev)/hFus
    if (f0 < 1._R8) then
      Tafter = Tmelt;                                    fAfter = f0;    phaseAfter = phPlateau
    else
      Tafter = Tmelt - (cpl*(Tmelt - Tev) - hFus)/cpSol; fAfter = 1._R8; phaseAfter = phSolid
    endif

  end subroutine nucleationJump


  !> Specific enthalpy with the latent heat of fusion; hOff is the table datum.
  pure function hSolid(T, f, phase, cpl, cpSol, hFus, Tmelt, hOff) result(h)
    real(R8), intent(in) :: T, f, cpl, cpSol, hFus, Tmelt, hOff
    integer,  intent(in) :: phase
    real(R8) :: h

    select case (phase)
    case (phPlateau); h = cpl*Tmelt + hOff - f*hFus
    case (phSolid);   h = cpl*Tmelt + hOff - hFus - cpSol*(Tmelt - T)
    case default;     h = cpl*T + hOff
    end select

  end function hSolid


  !> Latent part of the energy: e(T,f) - cpl*T = -f*L(T), L(T) = hFus - (cpl - cpSol)*(Tmelt - T).
  pure function solid_de(T, f, cpl, cpSol, hFus, Tmelt) result(de)
    real(R8), intent(in) :: T, f, cpl, cpSol, hFus, Tmelt
    real(R8) :: de

    de = -f*(hFus - (cpl - cpSol)*(Tmelt - T))

  end function solid_de


  !> Temperature, frozen fraction and nucleated fraction of a content of liquid-branch temperature Tliq = e/cpl and
  !  transported nucleated fraction chi: liquid from Tmelt (chi cleared), nucleated up to Tnuc (chi set), in between
  !  nucleated iff chi >= 1/2 (chi kept); a nucleated content is nucleationJump(Tliq), the lever rule below Tmelt.
  pure subroutine solid_state(Tliq, chi, cpl, cpSol, Tmelt, Tnuc, hFus, T, f, chiOut)
    real(R8), intent(in)  :: Tliq, chi, cpl, cpSol, Tmelt, Tnuc, hFus
    real(R8), intent(out) :: T, f, chiOut
    integer :: phase

    if (Tliq >= Tmelt) then
      T = Tliq; f = 0._R8; chiOut = 0._R8
    elseif (Tliq <= Tnuc) then
      call nucleationJump(Tliq, cpl, cpSol, Tmelt, hFus, T, f, phase); chiOut = 1._R8
    elseif (chi >= 0.5_R8) then
      call nucleationJump(Tliq, cpl, cpSol, Tmelt, hFus, T, f, phase); chiOut = chi
    else
      T = Tliq; f = 0._R8; chiOut = chi
    endif

  end subroutine solid_state

end module ICE_Lib_Solidification
