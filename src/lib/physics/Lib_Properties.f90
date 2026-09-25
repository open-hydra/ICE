!>@brief Condensed-material properties as functions of temperature.
!> A constant property returns its value with no arithmetic; a tabulated one is linear
!> between the integer-kelvin nodes and takes the end value outside them. The energy
!> e = h - hOff is extended linearly beyond both ends (through e(0) = 0 below), and T(e)
!> inverts it by bisection over the nodes. The tables must be allocated (use_table).
module ICE_Lib_Properties
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use ICE_Config_Types_m, only: condensed_phase_t
  implicit none
  private
  public :: mat_rho, mat_cp, mat_e, mat_T_from_e

contains

  pure function mat_rho(mat, T) result(rho)
    type(condensed_phase_t), intent(in) :: mat
    real(R8),                intent(in) :: T
    real(R8) :: rho
    if (.not. mat%use_table) then
      rho = mat%rho_al
    elseif (.not. mat%rho_varies) then
      rho = mat%rho_tab(mat%T_min)
    else
      rho = table_value(mat%rho_tab, mat%T_min, mat%T_max, T)
    endif
  end function mat_rho


  pure function mat_cp(mat, T) result(cp)
    type(condensed_phase_t), intent(in) :: mat
    real(R8),                intent(in) :: T
    real(R8) :: cp
    if (.not. mat%use_table) then
      cp = mat%cs_al
    elseif (.not. mat%cs_varies) then
      cp = mat%cs_tab(mat%T_min)
    else
      cp = table_value(mat%cs_tab, mat%T_min, mat%T_max, T)
    endif
  end function mat_cp


  !> Energy e = h - hOff [J/kg]: the table between the nodes, the line from e(0) = 0 to
  !  eTab(Tmin) below Tmin (the first segment when Tmin = 0), the last-segment slope above Tmax.
  pure function mat_e(mat, T) result(e)
    type(condensed_phase_t), intent(in) :: mat
    real(R8),                intent(in) :: T
    real(R8) :: e
    integer  :: i
    associate (tab => mat%e_tab, lo => mat%T_min, hi => mat%T_max)
      if (T /= T) then
        e = T
      elseif (T < real(lo, R8)) then
        if (lo > 0) then
          e = tab(lo)*(T/real(lo, R8))
        else
          e = (tab(lo+1) - tab(lo))*T
        endif
      elseif (T >= real(hi, R8)) then
        e = tab(hi) + (tab(hi) - tab(hi-1))*(T - real(hi, R8))
      else
        i = int(T)
        e = tab(i) + (tab(i+1) - tab(i))*(T - real(i, R8))
      endif
    end associate
  end function mat_e


  !> Inverse of mat_e: the same three pieces, the bracketing node found by bisection.
  pure function mat_T_from_e(mat, e) result(T)
    type(condensed_phase_t), intent(in) :: mat
    real(R8),                intent(in) :: e
    real(R8) :: T
    integer  :: i, j, k
    associate (tab => mat%e_tab, lo => mat%T_min, hi => mat%T_max)
      if (e /= e) then
        T = e
      elseif (e < tab(lo)) then
        if (lo > 0) then
          T = real(lo, R8)*(e/tab(lo))
        else
          T = e/(tab(lo+1) - tab(lo))
        endif
      elseif (e >= tab(hi)) then
        T = real(hi, R8) + (e - tab(hi))/(tab(hi) - tab(hi-1))
      else
        i = lo
        j = hi
        do while (j - i > 1)
          k = (i + j)/2
          if (tab(k) <= e) then
            i = k
          else
            j = k
          endif
        enddo
        T = real(i, R8) + (e - tab(i))/(tab(i+1) - tab(i))
      endif
    end associate
  end function mat_T_from_e


  !> Linear between the nodes, the end value outside; a NaN temperature returns NaN.
  pure function table_value(tab, lo, hi, T) result(v)
    integer,  intent(in) :: lo, hi
    real(R8), intent(in) :: tab(lo:hi), T
    real(R8) :: v
    real(R8) :: Tc
    integer  :: i
    if (T /= T) then
      v = T
      return
    endif
    Tc = min(max(T, real(lo, R8)), real(hi, R8))
    if (Tc >= real(hi, R8)) then
      v = tab(hi)
      return
    endif
    i  = int(Tc)
    v  = tab(i) + (tab(i+1) - tab(i))*(Tc - real(i, R8))
  end function table_value

end module ICE_Lib_Properties
