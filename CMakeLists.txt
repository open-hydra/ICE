# CMake project file for 'ICE'

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver ICE)

project (${solver} LANGUAGES Fortran)

include(CTest)

# The core-cycle counter (src/lib/diagnostic/ice_cycles.c) is the only C in the
# tree: perf_event_open has no Fortran binding. It is optional, because a
# toolchain that cannot link a C program must still be able to build ICE -- the
# timers then report wall-clock only. Configure with -DICE_CYCLES=OFF to skip
# it deliberately.
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

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # Search path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  OPTION(USE_MPI    "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization"          OFF)

  # Interprocedural optimisation across the modules. The flux path reaches the
  # reconstruction, the limiter and the Riemann solver through procedures of
  # other files, which -O3 alone never inlines. Off by default, because it moves
  # the last bits of the fields and the references were made without it;
  # measured 6-10 % off the flux region at one thread (ifx 2023.1, 2026-09-30).
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

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING. You should review this file and make sure the flags
  # are to your liking.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
  # Locate and set parallelization libraries. There are some CMake peculiarities
  # taken care of here, such as the fact that the FindOpenMP routine doesn't know
  # about Fortran.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # -------- Path fallbacks (standalone builds) --------
  if (NOT DEFINED ORION_PATH)
    set(ORION_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/")
  endif()
  if (NOT DEFINED FINER_PATH)
    set(FINER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/")
  endif()

  # -------- Helper: initialize a git submodule from the correct git root --------
  macro(ensure_submodule_initialized SUBMOD_PATH)
    if (NOT EXISTS "${SUBMOD_PATH}/CMakeLists.txt")
      message(STATUS "Submodule at ${SUBMOD_PATH} not initialized. Attempting git submodule update --init...")
      get_filename_component(_sm_parent "${SUBMOD_PATH}" DIRECTORY)
      execute_process(
        COMMAND git rev-parse --show-toplevel
        WORKING_DIRECTORY "${_sm_parent}"
        OUTPUT_VARIABLE _sm_git_root
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _sm_git_result)
      if(_sm_git_result EQUAL 0)
        file(RELATIVE_PATH _sm_rel "${_sm_git_root}" "${SUBMOD_PATH}")
        execute_process(
          COMMAND git submodule update --init "${_sm_rel}"
          WORKING_DIRECTORY "${_sm_git_root}"
          RESULT_VARIABLE _sm_init_result)
        if(NOT _sm_init_result EQUAL 0)
          message(FATAL_ERROR "Failed to initialize submodule: ${SUBMOD_PATH}")
        endif()
      else()
        message(FATAL_ERROR "Could not find git root for submodule: ${SUBMOD_PATH}")
      endif()
    endif()
  endmacro()

  # ORION
  if (NOT TARGET ORION)
    ensure_submodule_initialized("${ORION_PATH}")
    set(orionD "${ORION_PATH}" "${CMAKE_BINARY_DIR}/lib/ORION")
    add_subdirectory(${orionD})
  endif()

  # FiNeR
  if (NOT TARGET FiNeR)
    ensure_submodule_initialized("${FINER_PATH}")
    set(finerD "${FINER_PATH}" "${CMAKE_BINARY_DIR}/lib/third_party/FiNeR")
    add_subdirectory(${finerD})
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
