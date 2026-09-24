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
                 'Drag model, global for all families; required for a coupled run (none = not set)', &
                 'Newton, Stokes, Schlichting, Schiller-Naumann, Wen-Yu, Putnam, '// &
                 'Clift-Gauvin, Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, '// &
                 'Hermsen, none', .false.)
    call reg%add(section, 'heat-transfer', obj_time_scheme%heat, 'none', &
                 'Convective heat transfer model, global for all families; required for a coupled run '// &
                 '(none = not set)', &
                 'Stokes, JAXA1, JAXA2, JAXA3, Chang, Ranz-Marshall, Kavanau-Drake, none', &
                 .false.)
    call reg%add(section, 'evaporation', obj_time_scheme%evaporation, 'none', &
                 'Evaporation model, global for all families', &
                 'd2-law, CEM, CEM-B, ASM, TC, none', .false.)
    call reg%add(section, 'evaporation-interface', obj_time_scheme%interface_model, 'VLE', &
                 'Vapour-liquid interface: VLE equilibrium, or LK Langmuir-Knudsen '// &
                 'non-equilibrium; ignored when evaporation is none', &
                 'VLE, LK', .false.)
    call reg%add(section, 'evaporation-blowing', obj_time_scheme%blowing, 'none', &
                 'Stefan-blowing reduction of the convective heat; LK applies '// &
                 'Miller-Harstad-Bellan f2. Ignored under ASM and TC, which carry '// &
                 'their own gas-side heat', &
                 'LK, none', .false.)

    ! Condensed-material properties ----------------------
    call reg%add(section, 'density', obj_condensed%rho_al, '2700.0', &
                 'Condensed-material density [kg/m^3], used when no property table is given', &
                 '> 0',  .false.)
    call reg%add(section, 'specific-heat', obj_condensed%cs_al, '1598.0', &
                 'Condensed-material specific heat [J/(kg K)], used when no property table is given', &
                 '> 0',  .false.)
    call reg%add(section, 'latent-heat', obj_condensed%lv_al, '1.08e7', &
                 'Latent heat of vaporisation [J/kg]; the evaporation models use it '// &
                 'both as the energy sink and as the anchor of the saturation curve', &
                 '> 0',  .false.)
    call reg%add(section, 'emissivity', obj_condensed%emiss, '1.0', &
                 'Particle surface emissivity; 0 switches radiative exchange off', &
                 '>= 0', .false.)

    ! Vapour properties (only read when an evaporation model is selected) --
    call reg%add(section, 'vapour-molar-mass', obj_condensed%Mv, '26.98', &
                 'Molar mass of the vapour [kg/kmol]', '> 0', .false.)
    call reg%add(section, 'boiling-temperature', obj_condensed%Tboil, '2792.0', &
                 'Boiling temperature at 1 atm [K], the anchor of the '// &
                 'Clausius-Clapeyron saturation pressure', '> 0', .false.)
    call reg%add(section, 'vapour-specific-heat', obj_condensed%cpv, '0.0', &
                 'Specific heat of the vapour [J/(kg K)]; 0 falls back to the gas cp', &
                 '>= 0', .false.)
    call reg%add(section, 'lewis-number', obj_condensed%Le, '1.0', &
                 'Lewis number of the vapour in the gas, Le = k/(rho cp D)', '> 0', .false.)
    call reg%add(section, 'vapour-mass-fraction', obj_condensed%Yinf, '0.0', &
                 'Vapour mass fraction in the far-field gas; evaporation stops once '// &
                 'the surface value falls to it', '>= 0', .false.)
    call reg%add(section, 'evaporation-coefficient', obj_condensed%alphaE, '1.0', &
                 'Evaporation (accommodation) coefficient of the Langmuir-Knudsen '// &
                 'interface; unused under VLE', '> 0', .false.)

  end subroutine Register_Physics

end module ICE_Read_Physics
