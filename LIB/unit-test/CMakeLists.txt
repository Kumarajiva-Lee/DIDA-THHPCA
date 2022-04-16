cmake_minimum_required(VERSION 3.0)

set(CMAKE_MACOSX_RPATH ON)

project(FortranUnitTest LANGUAGES Fortran)

set(version 0.0.1)
 
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

set(sources
  src/assert_mod.F90
  src/test_common_mod.F90
  src/test_case_mod.F90
  src/test_suite_mod.F90
  src/unit_test.F90
)

add_library(fortran_unit_test STATIC ${sources})
