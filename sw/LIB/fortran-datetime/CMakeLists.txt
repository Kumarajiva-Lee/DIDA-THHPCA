cmake_minimum_required(VERSION 3.0)
set (CMAKE_MACOSX_RPATH ON)
project(fortran-datetime LANGUAGES Fortran)
set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")

set (version 0.0.2)
 
set (sources
  src/timedelta_mod.F90
  src/datetime_mod.F90
  src/datetime.F90
)

add_library (fortran_datetime STATIC ${sources})
