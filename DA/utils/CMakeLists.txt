cmake_minimum_required(VERSION 3.0)

project(da)

#find_package(MPI REQUIRED)
#include_directories("../rpe/modules")
#link_directories("../rpe/lib")

#enable_language(Fortran)
#set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_Fortran_COMPILER mpiifort)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fcray-pointer")

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
  ../EnDNN/src/module_da_parallel.F90
  ../EnDNN/src/module_da_coupler.F90
  ../CNN/c_model/main.c
  ../CNN/c_model/init_data.c
  ../CNN/c_model/predict.c
)

add_library(utils ${utils_files})
if(_is_sw_platform)
  target_link_libraries(utils atm coupler redis fiona flogger fortran_string netcdff netcdf hdf5_hl hdf5 z ${MPI_Fortran_LIBRARIES} )
  #target_link_libraries(utils atm coupler redis fiona flogger fortran_string ${MPI_Fortran_LIBRARIES} )
else()
  target_link_libraries(utils atm coupler redis fiona flogger fortran_string ${MPI_Fortran_LIBRARIES} )
endif()

add_library(letkf_core ${letkf_core_file})
if (_is_sw_platform)
  if (_is_sw_master)
      #target_link_libraries(letkf_core lapack refblas tmglib m)
    target_link_libraries(letkf_core letkf_half swhgemm mpal rpe lapack swblas m_slave tmglib)
  else()
    #target_link_libraries(letkf_core letkf_half swhgemm mpal rpe xMath m_slave)
    target_link_libraries(letkf_core letkf_half swhgemm mpal rpe swlapack swblas m_slave)
  endif()
else()
  target_link_libraries(letkf_core mpal rpe)
endif()

add_library(letkf ${letkf_files})
if(_is_sw_platform)
  if (_is_sw_master)
    target_link_libraries(letkf letkf_core atm utils atm flogger coupler redis mpal rpe m ${MPI_Fortran_LIBRARIES})
  else()
    target_link_libraries(letkf letkf_core atm utils atm flogger coupler redis mpal rpe m_slave ${MPI_Fortran_LIBRARIES})
  endif()
else()
  #add_library(letkf ${letkf_files} ${utils_files})
  target_link_libraries(letkf letkf_core atm utils atm flogger coupler redis mpal rpe ${MPI_Fortran_LIBRARIES})
endif()

add_library(da ${dnn_files})
#target_link_libraries(letkf atm utils coupler redis ${MPI_Fortran_LIBRARIES} )
if(_is_sw_platform)
    target_link_libraries(da letkf utils atm coupler redis m_slave darknet swdnn cblas_LINUX_SW blas_LINUX_SW swblas2 ${MPI_Fortran_LIBRARIES})
    #target_link_libraries(da letkf utils atm coupler redis lapack refblas tmglib ${MPI_Fortran_LIBRARIES})
else()
    target_link_libraries(da letkf utils atm coupler redis ${MPI_Fortran_LIBRARIES})
endif()
