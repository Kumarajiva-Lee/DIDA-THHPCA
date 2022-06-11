cmake_minimum_required(VERSION 3.0)

project(da)

#find_package(MPI REQUIRED)
#include_directories(${MPI_Fortran_INCLUDE_PATH})
#link_directories("/software/intel/2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64")

#enable_language(Fortran)
#set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_Fortran_COMPILER mpiifort)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

set(utils_files
  pro_info.f90
  da_namelist_mod.F90
  splitchar.f90
  da_output.f90
)

set(letkf_files
  ../LETKF/src/letkf.f90
)

set (dnn_files
  ../EnDNN/src/module_da.F90
  ../EnDNN/src/module_da_config.F90
  ../EnDNN/src/module_da_coupler.F90
  ../EnDNN/src/module_da_parallel.F90
)

add_library(utils ${utils_files})

#add_library(letkf ${letkf_files} ${utils_files})

add_library(da ${dnn_files} ${letkf_files})
target_link_libraries(da coupler redis utils ${MPI_Fortran_LIBRARIES})
