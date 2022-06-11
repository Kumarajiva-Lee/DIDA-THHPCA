#cmake_minimum_required(VERSION 3.0)
#
#project(flogger LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  add_definitions(-DINTEL)
endif ()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")

#get_directory_property(parent_dir PARENT_DIRECTORY)
#if (EXISTS ${PROJECT_SOURCE_DIR}/LIB/container/CMakeLists.txt AND NOT parent_dir)
#    add_subdirectory(LIB/container)
#endif ()
#if (EXISTS ${PROJECT_SOURCE_DIR}/LIB/string/CMakeLists.txt AND NOT parent_dir)
#    add_subdirectory(LIB/string)
#endif ()

set(source_files
  src/face_mod.F90
  src/log_mod.F90
  src/flogger.F90
  src/c_code.c
)

add_library(flogger ${source_files})
target_link_libraries(flogger coupler fortran_container fortran_string  ${MPI_Fortran_LIBRARIES})

if (flogger_unit)
  add_executable(test_flogger src/test_flogger.F90)
  target_link_libraries(test_flogger flogger)
endif ()
