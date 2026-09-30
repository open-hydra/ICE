function(add_fortran_executable target source)

  add_executable(${target} ${source})
  set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${target} PRIVATE ICEL FiNeR::FiNeR)

  if (USE_TECIO)
    target_link_libraries(${target} PRIVATE -pthread)
  endif()

  target_include_directories(${target} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

  if (USE_OPENMP)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  if (USE_MPI)
    target_compile_definitions(${target} PRIVATE USE_MPI)
    target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
  endif()

endfunction()

# Where the executables go. A cache entry, so that two configurations of one
# tree (serial and MPI, release and strict floating point) can each keep their
# own binary: with the plain variable every build wrote the same bin/ICE.
set(ICE_BIN_DIR "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Directory that receives the ICE executables")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ICE_BIN_DIR})

add_fortran_executable(${solver} "main.f90")
add_fortran_executable(DocGen "docgen.f90")