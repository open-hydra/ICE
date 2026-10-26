file(GLOB_RECURSE ICE_SRC
    CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.f90
)

# The cycle counter is the only C in the tree; it needs perf_event_open, which
# has no Fortran binding. It compiles to a no-op away from Linux, and is left
# out entirely when no C compiler was found (see the top-level CMakeLists).
set(ICE_SRC_C)
IF(ICE_HAVE_CYCLES)
    file(GLOB_RECURSE ICE_SRC_C
        CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c
    )
ENDIF()

add_library(ICEL ${ICE_SRC} ${ICE_SRC_C})
add_library(ICEL::ICEL ALIAS ICEL)
target_link_libraries(ICEL PRIVATE FiNeR::FiNeR ORION)

target_include_directories(ICEL INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(ICE_HAVE_CYCLES)
    target_compile_definitions(ICEL PUBLIC ICE_HAVE_CYCLES)
ENDIF()

IF(USE_OPENMP)
    target_link_libraries(ICEL PUBLIC OpenMP::OpenMP_Fortran)
ENDIF(USE_OPENMP)

IF(ICE_RACE_PROBE)
    target_compile_definitions(ICEL PRIVATE ICE_RACE_PROBE)
ENDIF()

IF(USE_MPI)
    target_compile_definitions(ICEL PUBLIC USE_MPI)
    target_link_libraries(ICEL PUBLIC MPI::MPI_Fortran)
    # An mpi module made for a compiler that cannot ignore the type of a choice buffer declares no
    # interface for the routines that take one (MPI 4.1, Section 19.1.3), as the gfortran modules of
    # Intel MPI do: gfortran >= 10 then stops on calls of one such routine with buffers of different
    # types or ranks in one file. Mod_MPI.f90 makes such calls (MPI_BCAST of a LOGICAL and of an
    # INTEGER, MPI_ALLREDUCE of a scalar and of an array); -fallow-argument-mismatch turns those errors
    # into warnings for that file only, and the other files keep the full check.
    IF(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
        set_source_files_properties(parallel/Mod_MPI.f90 PROPERTIES COMPILE_OPTIONS -fallow-argument-mismatch)
    ENDIF()
ENDIF(USE_MPI)
