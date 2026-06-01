module ICE_Mod_Space
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Space_Scheme

contains

  subroutine Setup_Space_Scheme()
    use ICE_Config_Types_m, only: obj_space_scheme
    implicit none

    obj_space_scheme%SD = (index(obj_space_scheme%space_reconstruction, 'SD') > 0)

    if (index(obj_space_scheme%space_reconstruction, 'MUSCL') > 0) then
      if (trim(obj_space_scheme%flux_limiter) == 'none') then
        obj_space_scheme%flux_limiter = 'VANLEER'
        write(*, '(A)') ' [WARNING] MUSCL specified without flux-limiter. Van Leer by default.'
      end if
    end if

  end subroutine Setup_Space_Scheme

end module ICE_Mod_Space
