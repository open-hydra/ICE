module ICE_Read_Physics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m
  use ICE_Config_Types_m, only: ini_condensed, obj_time_scheme
  use ICE_Input_Registry

  implicit none
  private
  public :: Register_Physics

contains

  subroutine Register_Physics()
    implicit none
    character(len=:), allocatable :: section

    ini_condensed%warning_message = 'none'
    ini_condensed%error_message   = 'none'
    ini_condensed%description     = 'none'

    !! ------------------------------------------------------
    !! Condensed Phase Material Properties ------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Physics'

    ! Interphase exchange models -------------------------
    call reg%add(section, 'drag', obj_time_scheme%drag, 'none', &
                 'Drag model, global for all families; required for a coupled run (none = not set, '// &
                 'NoDrag = no momentum exchange)', &
                 'Newton, Stokes, Schlichting, Schiller-Naumann, Wen-Yu, Putnam, '// &
                 'Clift-Gauvin, Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, '// &
                 'Hermsen, NoDrag, none', .false.)
    call reg%add(section, 'heat-transfer', obj_time_scheme%heat, 'none', &
                 'Convective heat transfer model, global for all families; required for a coupled run '// &
                 '(none = not set, NoHeat = no convective exchange; Chang stops with a pointer to JAXA3, '// &
                 'which is its formula)', &
                 'Stokes, JAXA1, JAXA2, JAXA3, JAXA4, Chang, Ranz-Marshall, Kavanau-Drake, NoHeat, none', &
                 .false.)
    call reg%add(section, 'evaporation', obj_time_scheme%evaporation, 'none', &
                 'Evaporation model, the default of every material; a phase-file token overrides it', &
                 'd2-law, CEM, CEM-B, ASM, TC, none', .false.)
    call reg%add(section, 'evaporation-interface', obj_time_scheme%interface_model, 'VLE', &
                 'Vapour-liquid interface: VLE equilibrium, or LK Langmuir-Knudsen '// &
                 'non-equilibrium; ignored when evaporation is none; a phase-file token overrides it', &
                 'VLE, LK', .false.)
    call reg%add(section, 'evaporation-blowing', obj_time_scheme%blowing, 'none', &
                 'Stefan-blowing reduction of the convective heat; LK applies '// &
                 'Miller-Harstad-Bellan f2. Ignored under ASM and TC, which carry '// &
                 'their own gas-side heat', &
                 'LK, none', .false.)

    ! Condensed-material properties ----------------------
    call reg%add(section, 'density', ini_condensed%rho_al, '2700.0', &
                 'Condensed-material density [kg/m^3], used when no property table is given; with a table '// &
                 'it may be left out, and if given it must equal the table''s constant density'//'; one value per material', &
                 '> 0',  .false., per_material=.true.)
    call reg%add(section, 'specific-heat', ini_condensed%cs_al, '1598.0', &
                 'Condensed-material specific heat [J/(kg K)], used when no property table is given; with a '// &
                 'table it may be left out, and if given it must equal the table''s constant cp'//'; one value per material', &
                 '> 0',  .false., per_material=.true.)
    call reg%add(section, 'latent-heat', ini_condensed%lv_al, '1.08e7', &
                 'Latent heat of vaporisation [J/kg]; the evaporation models use it '// &
                 'both as the energy sink and as the anchor of the Clausius-Clapeyron saturation curve, '// &
                 'which a Psat column in the property table replaces'//'; one value per material', &
                 '> 0',  .false., per_material=.true.)
    call reg%add(section, 'emissivity', ini_condensed%emiss, '1.0', &
                 'Particle surface emissivity; 0 switches radiative exchange off'//'; one value per material', &
                 '>= 0', .false., per_material=.true.)

    ! Vapour properties (only read when an evaporation model is selected) --
    call reg%add(section, 'vapour-molar-mass', ini_condensed%Mv, '26.98', &
                 'Molar mass of the vapour [kg/kmol]'//'; one value per material', '> 0', .false., per_material=.true.)
    call reg%add(section, 'boiling-temperature', ini_condensed%Tboil, '2792.0', &
                 'Boiling temperature at 1 atm [K], the anchor of the '// &
                 'Clausius-Clapeyron saturation pressure; with a Psat column in the property table '// &
                 'it must lie in [Tmin, Tmax-1] of the table, where the column must give 0.5 to 2 atm'//'; one value per material', &
                 '> 0', .false., per_material=.true.)
    call reg%add(section, 'vapour-specific-heat', ini_condensed%cpv, '0.0', &
                 'Specific heat of the vapour [J/(kg K)]; 0 falls back to the gas cp'//'; one value per material', &
                 '>= 0', .false., per_material=.true.)
    call reg%add(section, 'lewis-number', ini_condensed%Le, '1.0', &
                 'Lewis number of the vapour in the gas, Le = k/(rho cp D)'//'; one value per material', '> 0', .false., per_material=.true.)
    call reg%add(section, 'vapour-mass-fraction', ini_condensed%Yinf, '0.0', &
                 'Vapour mass fraction in the far-field gas; evaporation stops once '// &
                 'the surface value falls to it'//'; one value per material', '>= 0', .false., per_material=.true.)
    call reg%add(section, 'evaporation-coefficient', ini_condensed%alphaE, '1.0', &
                 'Evaporation (accommodation) coefficient of the Langmuir-Knudsen '// &
                 'interface; unused under VLE; the default of alpha-e'//'; one value per material', '> 0', .false., per_material=.true.)

  end subroutine Register_Physics

end module ICE_Read_Physics
