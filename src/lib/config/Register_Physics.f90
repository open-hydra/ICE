module ICE_Read_Physics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m
  use ICE_Config_Types_m, only: obj_condensed, obj_time_scheme
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

    ! Interphase exchange models -------------------------
    call reg%add(section, 'drag', obj_time_scheme%drag, 'none', &
                 'Drag model, global for all families', &
                 'Newton, Stokes, Schlichting, Schiller-Naumann, Wen-Yu, Putnam, '// &
                 'Clift-Gauvin, Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, '// &
                 'Hermsen, none', .false.)
    call reg%add(section, 'heat-transfer', obj_time_scheme%heat, 'none', &
                 'Convective heat transfer model, global for all families', &
                 'Stokes, JAXA1, JAXA2, JAXA3, Chang, Ranz-Marshall, Kavanau-Drake, none', &
                 .false.)

    ! Condensed-material properties ----------------------
    call reg%add(section, 'density', obj_condensed%rho_al, '2700.0', &
                 'Condensed-material density [kg/m^3], used when no property table is given', &
                 '> 0',  .false.)
    call reg%add(section, 'specific-heat', obj_condensed%cs_al, '1598.0', &
                 'Condensed-material specific heat [J/(kg K)], used when no property table is given', &
                 '> 0',  .false.)
    call reg%add(section, 'latent-heat', obj_condensed%lv_al, '1.08e7', &
                 'Latent heat of vaporisation [J/kg]; only acts through the mass-transfer term', &
                 '> 0',  .false.)
    call reg%add(section, 'combustion-energy', obj_condensed%q_al, '9.53e6', &
                 'Heat of combustion [J/kg]; only acts through the mass-transfer term', &
                 '> 0',  .false.)
    call reg%add(section, 'emissivity', obj_condensed%emiss, '1.0', &
                 'Particle surface emissivity; 0 switches radiative exchange off', &
                 '>= 0', .false.)

  end subroutine Register_Physics

end module ICE_Read_Physics
