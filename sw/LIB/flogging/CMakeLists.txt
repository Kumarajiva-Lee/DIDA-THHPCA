#cmake_minimum_required(VERSION 2.8)

project(flogging LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  add_definitions(-DINTEL)
endif ()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -O0 -cpp -Iinclude ")

include_directories(include)

set(source_files
  src/flogging.F90
  src/vt100.F90
)

add_library(flogging ${source_files})
target_link_libraries(flogging coupler)

if (RUN_TEST)
  add_executable(test_flogging test/test.f90)
  target_link_libraries(test_flogger flogger)
endif ()
