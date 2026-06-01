module ICE_Read_Physics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m
  use ICE_Config_Types_m, only: obj_condensed
  use ICE_Input_Registry

  implicit none
  private
  public :: Register_Physics

contains

  subroutine Register_Physics()
    implicit none
    character(len=:), allocatable :: section

    obj_condensed%warning_message = 'none'
    obj_condensed%error_message   = 'none'
    obj_condensed%description     = 'none'

    !! ------------------------------------------------------
    !! Condensed Phase Material Properties ------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Physics'

    call reg%add(section, 'rho',   obj_condensed%rho_al, &
                 '2700.0', 'Aluminium density [kg/m^3]',         '> 0',  .false.)
    call reg%add(section, 'cs',    obj_condensed%cs_al,  &
                 '1598.0', 'Aluminium specific heat [J/(kg K)]', '> 0',  .false.)
    call reg%add(section, 'lv',    obj_condensed%lv_al,  &
                 '1.08e7', 'Aluminium latent heat [J/kg]',       '> 0',  .false.)
    call reg%add(section, 'q',     obj_condensed%q_al,   &
                 '9.53e6', 'Aluminium combustion energy [J/kg]', '> 0',  .false.)
    call reg%add(section, 'emiss', obj_condensed%emiss,  &
                 '1.0',    'Particle surface emissivity [-]',    '>= 0', .false.)

  end subroutine Register_Physics

end module ICE_Read_Physics
