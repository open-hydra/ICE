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

IF(USE_MPI)
    target_compile_definitions(ICEL PUBLIC USE_MPI)
    target_link_libraries(ICEL PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)
