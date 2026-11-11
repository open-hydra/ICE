# CMake project file for 'ICE'

cmake_minimum_required(VERSION 3.23)

set(solver ICE)
project (${solver} LANGUAGES Fortran)
include(CTest)

# The core-cycle counter (src/lib/diagnostic/ice_cycles.c) is the only C in the tree. It is optional, because a toolchain that cannot link a C program must still be able to build ICE -- the timers then report wall-clock only. 
# Configure with -DICE_CYCLES=OFF to skip it deliberately.
OPTION(ICE_CYCLES "Count core-cycles per iteration (needs a working C compiler)" ON)
set(ICE_HAVE_CYCLES OFF)
IF(ICE_CYCLES)
  include(CheckLanguage)
  check_language(C)
  IF(CMAKE_C_COMPILER)
    enable_language(C)
    set(ICE_HAVE_CYCLES ON)
  ELSE()
    message(STATUS "No working C compiler: core-cycle counting disabled, "
                   "timers will report wall-clock only")
  ENDIF()
ENDIF()

# Libraries requirements
# - ORION. 
# Required by a range of releases with semver_dependency of cmake/SemVerDeps.cmake. 
# Standalone, the release is checked out in lib/ORION, out of the build folder, and kept in deps.lock (both git ignored); -DSEMVER_UPDATE=ALL selects the newest one again.
set(ICE_ORION_VERSION 2.2)
set(ICE_ORION_URL https://github.com/MarcoGrossi92/ORION.git)
# - FiNeR.
# Directly pinned at one release (bump the gitlink and the tag together).
set(ICE_FINER_TAG v2.3.0)
# FiNeR's own dependencies (cmake/FiNeRDependencies.cmake). 
# FREE fetches them as FiNeR's author does (`fobis fetch`)
# BACKUP checks out the releases below, which ICE builds and passes its tests with.
# When an upstream change breaks ICE, set BACKUP here, or configure one build with -DICE_FINER_DEPS=BACKUP (kept in its cache until -UICE_FINER_DEPS).
if (NOT DEFINED ICE_FINER_DEPS)
  set(ICE_FINER_DEPS "FREE")
endif()
set(ICE_FINER_BACKUP_PINS
  "PENF=v2.0.1"
  "StringiFor=v1.8.0"
  "FACE=v1.1.17"
  "FLAP=v2.5.3"
  "BeFoR64=v1.2.4")

# If ICE is built as a standalone project, the following code is applied.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # Search path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  OPTION(USE_MPI    "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization"          OFF)

  # Binary Tecplot (.szplt) through TecIO, which ORION builds as a static C++ library. The archive brings no C++ runtime, and CMake adds that runtime to a Fortran link only where CXX is enabled.
  OPTION(USE_TECIO  "Binary Tecplot I/O through TecIO (needs a C++ compiler)" OFF)
  if (USE_TECIO)
    enable_language(CXX)
  endif()

  # Interprocedural optimisation across the modules.
  # Off by default, because it moves the last bits of the fields and the references were made without it.
  OPTION(ICE_ENABLE_IPO "Interprocedural optimisation (moves round-off: re-baseline)" OFF)
  if (ICE_ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ice_ipo_ok OUTPUT ice_ipo_msg LANGUAGES Fortran)
    if (ice_ipo_ok)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    else()
      message(WARNING "ICE_ENABLE_IPO: not supported by this toolchain: ${ice_ipo_msg}")
    endif()
  endif()

  # A test build only: the flux kernel's thread-group rounds NaN-fill their shared
  # k-face planes and delay one writer, so a missing barrier between the k-faces
  # and the sweep stops the run (test Equiv3DGroupsRace, label race). Slower, and
  # never for production.
  OPTION(ICE_RACE_PROBE "Race probe of the grouped flux rounds (test build only)" OFF)

  # This INCLUDE statement executes code that sets the compile flags for DEBUG, RELEASE, and TESTING.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
  # Locate and set parallelization libraries.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # Libraries building
  # An external clone given by a relative path is taken from the source folder
  foreach(_name ORION FINER)
    if (FETCHCONTENT_SOURCE_DIR_${_name})
      get_filename_component(FETCHCONTENT_SOURCE_DIR_${_name} "${FETCHCONTENT_SOURCE_DIR_${_name}}" ABSOLUTE
                             BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    endif()
  endforeach()
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
  semver_dependency(ORION GIT_REPOSITORY ${ICE_ORION_URL} VERSION ${ICE_ORION_VERSION} SOURCE_DIR lib/ORION)
  # FiNeR and its dependencies, the checks of both libraries, and FetchContent_MakeAvailable
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/AddLibraries.cmake)

else()

  # Part of a larger project, which gives ICE its ORION and FiNeR. 
  # When that project selected ORION with SemVerDeps, its release is checked against ICE's requirement; nothing is added here.
  get_property(_ice_orion_selected GLOBAL PROPERTY SEMVER_ORION_SOURCE SET)
  if (_ice_orion_selected)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
    semver_dependency(ORION GIT_REPOSITORY ${ICE_ORION_URL} VERSION ${ICE_ORION_VERSION})
  endif()

endif()

# Build library
add_subdirectory(src/lib)

# Build application (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)

  # Test suite: ctest --test-dir <build> (see test/CMakeLists.txt)
  enable_testing()
  add_subdirectory(test)
endif()
