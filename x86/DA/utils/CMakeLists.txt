cmake_minimum_required(VERSION 3.0)

project(da)

#find_package(MPI REQUIRED)
#include_directories("../rpe/modules")
#link_directories("../rpe/lib")

#enable_language(Fortran)
#set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_Fortran_COMPILER mpiifort)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fp-model precise")

set(utils_files
    pro_info.F90
    da_namelist_mod.F90
    splitchar.F90
    var_info.F90
    vert_coord_mod.F90
    randn_mod.F90
    da_output.F90)

set(letkf_core_file
  ../LETKF/src/letkf_core.F90
  ../LETKF/src/letkf_core_mpal.F90)

set(letkf_files
  ../LETKF/src/letkf.F90
  ../LETKF/src/read_obs_mod.F90
#  ../LETKF/src/obsmaker_mod.F90
)

set (dnn_files
  ../EnDNN/src/module_da.F90
  ../EnDNN/src/module_da_config.F90
  ../EnDNN/src/module_da_coupler.F90
  ../EnDNN/src/module_da_parallel.F90
)

add_library(utils ${utils_files})
target_link_libraries(utils atm coupler redis fiona flogger fortran_string ${MPI_Fortran_LIBRARIES} )

add_library(letkf_core ${letkf_core_file})
target_link_libraries(letkf_core flogger)

add_library(letkf ${letkf_files})
#add_library(letkf ${letkf_files} ${utils_files})
target_link_libraries(letkf letkf_core atm utils atm coupler redis ${MPI_Fortran_LIBRARIES} )
add_library(da ${dnn_files})
#target_link_libraries(letkf atm utils coupler redis ${MPI_Fortran_LIBRARIES} )
target_link_libraries(da letkf letkf_core utils atm flogger coupler redis mpal rpe cadnaF ${MPI_Fortran_LIBRARIES})
