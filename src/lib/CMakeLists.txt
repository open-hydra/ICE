file(GLOB_RECURSE ICE_SRC
    CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.f90
)

add_library(ICEL ${ICE_SRC})
add_library(ICEL::ICEL ALIAS ICEL)
target_link_libraries(ICEL PRIVATE FiNeR::FiNeR ORION)

target_include_directories(ICEL INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    target_link_libraries(ICEL PUBLIC OpenMP::OpenMP_Fortran)
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(ICEL PUBLIC USE_MPI)
    target_link_libraries(ICEL PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)
