cmake_minimum_required(VERSION 3.0)

project(letkf)

#find_package(MPI REQUIRED)
#include_directories(${MPI_Fortran_INCLUDE_PATH})
#link_directories("/software/intel/2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64")

#enable_language(Fortran)
#set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_Fortran_COMPILER mpiifort)
set(CMAKE_Fortran_FLAGS "-g -traceback -lp0 -cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

set(namelist_file
  src/utils/namelist_mod.F90)

set(pro_info_file
  src/pro_info.f90
)

set(source_files
  src/letkf.F90
)

add_library(namelist ${namelist_file})

add_library(pro_info ${pro_info_file})

add_library(letkf ${source_files})
target_link_libraries(letkf pro_info namelist ${MPI_Fortran_LIBRARIES})

#add_executable(test.exe unit/main.f90)
#target_link_libraries(test.exe letkf ${MPI_Fortran_LIBRARIES})
